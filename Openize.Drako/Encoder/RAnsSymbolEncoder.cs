using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Encoder
{

    /// <summary>
    /// A helper class for encoding symbols using the rANS algorithm (see ans.h).
    /// The class can be used to initialize and encode probability table needed by
    /// rANS, and to perform encoding of symbols into the provided EncoderBuffer.
    /// </summary>
    class RAnsSymbolEncoder : RAnsBitCodec
    {

        private int maxSymbols;
        private int ransPrecisionBits;
        private int ransPrecision;

        private RansSym[] probabilityTable;
        /// <summary>
        /// The number of symbols in the input alphabet.
        /// </summary>
        uint numSymbols;
        /// <summary>
        /// Expected number of bits that is needed to encode the input.
        /// </summary>
        ulong numExpectedBits;

        RAnsEncoder ans;
        /// <summary>
        /// Initial offset of the encoder buffer before any ans data was encoded.
        /// </summary>
        ulong bufferOffset;

        public RAnsSymbolEncoder(int maxSymbolBitLength, ulong[] frequencies, EncoderBuffer buffer)
        {

            maxSymbols = 1 << maxSymbolBitLength;
            ransPrecisionBits = ComputeRAnsPrecisionFromMaxSymbolBitLength(maxSymbolBitLength);
            ransPrecision = 1 << ransPrecisionBits;

            ans = new RAnsEncoder(ransPrecisionBits);

            // Compute the total of the input frequencies.
            ulong totalFreq = 0;
            uint maxValidSymbol = 0;
            for (uint i = 0; i < frequencies.Length; ++i)
            {
                totalFreq += frequencies[i];
                if (frequencies[i] > 0)
                    maxValidSymbol = i;
            }
            uint numSymbols = maxValidSymbol + 1;
            this.numSymbols = numSymbols;
            probabilityTable = new RansSym[numSymbols];
            double totalFreqD = totalFreq;
            double ransPrecisionD = ransPrecision;
            // Compute probabilities by rescaling the normalized frequencies into interval
            // [1, ransPrecision - 1]. The total probability needs to be equal to
            // ransPrecision.
            int totalRansProb = 0;
            for (int i = 0; i < numSymbols; ++i)
            {
                ulong freq = frequencies[i];

                // Normalized probability.
                double prob = freq / totalFreqD;

                // RAns probability in range of [1, ransPrecision - 1].
                uint ransProb = (uint) (prob * ransPrecisionD + 0.5f);
                if (ransProb == 0 && freq > 0)
                    ransProb = 1;
                probabilityTable[i].prob = ransProb;
                totalRansProb += (int) ransProb;
            }
            // Because of rounding errors, the total precision may not be exactly accurate
            // and we may need to adjust the entries a little bit.
            if (totalRansProb != ransPrecision)
            {
                int[] sortedProbabilities = new int[numSymbols];
                for (int i = 0; i < numSymbols; ++i)
                {
                    sortedProbabilities[i] = i;
                }
                Array.Sort(sortedProbabilities, delegate(int a, int b)
                {
                    return probabilityTable[a].prob.CompareTo(probabilityTable[b].prob);
                });
                if (totalRansProb < ransPrecision)
                {
                    // This happens rather infrequently, just add the extra needed precision
                    // to the most frequent symbol.
                    probabilityTable[sortedProbabilities.Length - 1].prob +=
                        (uint) (ransPrecision - totalRansProb);
                }
                else
                {
                    // We have over-allocated the precision, which is quite common.
                    // Rescale the probabilities of all symbols.
                    int error = totalRansProb - ransPrecision;
                    while (error > 0)
                    {
                        double actTotalProbD = (double) totalRansProb;
                        double actRelErrorD = ransPrecisionD / actTotalProbD;
                        for (int j = (int) (numSymbols - 1); j > 0; --j)
                        {
                            int symbolId = sortedProbabilities[j];
                            if (probabilityTable[symbolId].prob <= 1)
                            {
                                if (j == numSymbols - 1)
                                    return; // Most frequent symbol would be empty.
                                break;
                            }
                            int newProb = (int) Math.Floor(actRelErrorD * probabilityTable[symbolId].prob);
                            int fix = (int) probabilityTable[symbolId].prob - newProb;
                            if (fix == 0u)
                                fix = 1;
                            if (fix >= (int) probabilityTable[symbolId].prob)
                                fix = (int) probabilityTable[symbolId].prob - 1;
                            if (fix > error)
                                fix = error;
                            probabilityTable[symbolId].prob -= (uint) fix;
                            totalRansProb -= fix;
                            error -= fix;
                            if (totalRansProb == ransPrecision)
                                break;
                        }
                    }
                }
            }

            // Compute the cumulative probability (cdf).
            uint totalProb = 0;
            for (int i = 0; i < numSymbols; ++i)
            {
                probabilityTable[i].cumProb = totalProb;
                totalProb += probabilityTable[i].prob;
            }
            if (totalProb != ransPrecision)
                throw new Exception("Failed to initialize RAns symbol encoder");

            // Estimate the number of bits needed to encode the input.
            // From Shannon entropy the total number of bits N is:
            //   N = -sum{i : allSymbols}(F(i) * log2(P(i)))
            // where P(i) is the normalized probability of symbol i and F(i) is the
            // symbol's frequency in the input data.
            double numBits = 0;
            for (int i = 0; i < numSymbols; ++i)
            {
                if (probabilityTable[i].prob == 0)
                    continue;
                double normProb = probabilityTable[i].prob / ransPrecisionD;
                numBits += frequencies[i] * Math.Log(normProb, 2);
            }
            numExpectedBits = (ulong) Math.Ceiling(-numBits);
            EncodeTable(buffer);
        }

        private void EncodeTable(EncoderBuffer buffer)
        {
            Encoding.EncodeVarint((uint) numSymbols, buffer);
            // Use varint encoding for the probabilities (first two bits represent the
            // number of bytes used - 1).
            for (uint i = 0; i < numSymbols; ++i)
            {
                uint prob = probabilityTable[i].prob;
                int numExtraBytes = 0;
                if (prob >= (1 << 6))
                {
                    numExtraBytes++;
                    if (prob >= (1 << 14))
                    {
                        numExtraBytes++;
                        if (prob >= (1 << 22))
                        {
                            numExtraBytes++;
                        }
                    }
                }

                if (prob == 0)
                {

                    // When the probability of the symbol is 0, set the first two bits to 1
                    // (unique identifier) and use the remaining 6 bits to store the offset
                    // to the next symbol with non-zero probability.
                    uint offset = 0;
                    for (; offset < (1 << 6) - 1; ++offset)
                    {
                        // Note: we don't have to check whether the next symbol id is larger
                        // than num_symbols_ because we know that the last symbol always has
                        // non-zero probability.
                        uint next_prob = probabilityTable[i + offset + 1].prob;
                        if (next_prob > 0)
                        {
                            break;
                        }
                    }

                    buffer.Encode((byte)((offset << 2) | 3));
                    i += offset;
                }
                else
                {
                    // Encode the first byte (including the number of extra bytes).
                    buffer.Encode((byte) ((prob << 2) | (uint) (numExtraBytes & 3)));
                    // Encode the extra bytes.
                    for (int b = 0; b < numExtraBytes; ++b)
                    {
                        buffer.Encode((byte) (prob >> (8 * (b + 1) - 2)));
                    }
                }
            }
        }

        public void StartEncoding(EncoderBuffer buffer)
        {
            // Allocate extra storage just in case.
            ulong requiredBits = 2 * numExpectedBits + 32;

            bufferOffset = (ulong) buffer.Bytes;
            ulong requiredBytes = (requiredBits + 7) / 8;
            buffer.Resize((int) bufferOffset + (int) requiredBytes + 8);
            var data = buffer.Data;
            // Offset the encoding by sizeof(bufferOffset). We will use this memory to
            // store the number of encoded bytes.
            ans.Reset(new BytePointer(data, (int) bufferOffset));
        }

        public void EncodeSymbol(int symbol)
        {
            ans.Write(probabilityTable[symbol]);
        }

        public void EndEncoding(EncoderBuffer buffer)
        {

            int src = (int) bufferOffset;
            // char *const src = const_cast<char *>(buffer->data()) + buffer_offset_; 
            // TODO(fgalligan): Look into changing this to uint32_t as write_end()
            // returns an int.
            //const uint64_t bytes_written = static_cast<uint64_t>(ans_.write_end());
            ulong bytes_written = (ulong) ans.writeEnd();

            EncoderBuffer var_size_buffer = new EncoderBuffer();
            Encoding.EncodeVarint(bytes_written, var_size_buffer);
            int size_len = var_size_buffer.Bytes;
            int dst = src + size_len;
            Array.Copy(buffer.Data, src, buffer.Data, dst, (int) bytes_written);

            // Store the size of the encoded data.
            //memcpy(src, var_size_buffer.data(), size_len);
            Array.Copy(var_size_buffer.Data, 0, buffer.Data, src, size_len);

            // Resize the buffer to match the number of encoded bytes.
            buffer.Resize((int) bufferOffset + (int) bytes_written + size_len);
        }
    }
}
