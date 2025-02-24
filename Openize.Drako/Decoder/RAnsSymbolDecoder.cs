using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Decoder
{
    class RAnsSymbolDecoder : RAnsBitCodec
    {
        private int maxSymbolBitLength;
        private int maxSymbols;
        private int ransPrecision;
        private uint[] probabilityTable;
        private int numSymbols;
        private RAnsDecoder ans;

        public RAnsSymbolDecoder(int maxSymbolBitLength)
        {
            this.maxSymbolBitLength = maxSymbolBitLength;
            this.maxSymbols = 1 << maxSymbolBitLength;
            int ransPrecisionBits = ComputeRAnsPrecisionFromMaxSymbolBitLength(maxSymbolBitLength);
            this.ransPrecision = 1 << ransPrecisionBits;
            ans= new RAnsDecoder(ransPrecisionBits);
        }

        public void Create(DecoderBuffer buffer)
        {
            if (buffer.BitstreamVersion == 0)
                throw DracoUtils.Failed();
            // Decode the number of alphabet symbols.
            if (buffer.BitstreamVersion < 20)
            {
                numSymbols = buffer.DecodeI32();
            }
            else
            {
                uint n = Decoding.DecodeVarintU32(buffer);
                numSymbols = (int)n;
            }
            probabilityTable = new uint[numSymbols];
            if (numSymbols == 0)
                return;
            // Decode the table.
            for (int i = 0; i < numSymbols; ++i)
            {
                uint prob = 0;
                // Decode the first byte and extract the number of extra bytes we need to
                // get.
                byte byteProb = buffer.DecodeU8();
                int token = byteProb & 3;
                if (token == 3)
                {
                    var offset = byteProb >> 2;
                    if (i + offset >= numSymbols)
                        throw DracoUtils.Failed();
                    // Set zero probability for all symbols in the specified range.
                    for (int j = 0; j < offset + 1; ++j)
                    {
                        probabilityTable[i + j] = 0;
                    }

                    i += offset;
                }
                else
                {
                    int extraBytes = byteProb & 3;
                    prob = (uint) (byteProb >> 2);
                    for (int b = 0; b < extraBytes; ++b)
                    {
                        byte eb = buffer.DecodeU8();
                        // Shift 8 bits for each extra byte and subtract 2 for the two first bits.
                        prob |= (uint) (eb) << (8 * (b + 1) - 2);
                    }
                }
                probabilityTable[i] = prob;
            }
            if (!ans.BuildLookupTable(probabilityTable, numSymbols))
                throw DracoUtils.Failed();
        }

        public void StartDecoding( DecoderBuffer buffer)
        {
            long bytesEncoded;
            // Decode the number of bytes encoded by the encoder.
            if (buffer.BitstreamVersion < 20)
            {
                bytesEncoded = buffer.DecodeI64();
            }
            else
            {
                ulong n = Decoding.DecodeVarintU64(buffer);
                bytesEncoded = (long)n;
            }

            if (bytesEncoded > buffer.RemainingSize)
                throw DracoUtils.Failed();
            BytePointer dataHead = buffer.Pointer + buffer.DecodedSize;
            // Advance the buffer past the rANS data.
            buffer.Advance((int)bytesEncoded);
            if (ans.readInit(dataHead, (int)bytesEncoded) != 0)
                throw DracoUtils.Failed();
        }

        public int NumSymbols
        {
            get { return numSymbols; }
        }

        public int DecodeSymbol()
        {
            return ans.Read();
        }

        public void EndDecoding()
        {
            ans.ReadEnd();
        }
    }
}
