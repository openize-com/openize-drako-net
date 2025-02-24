using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Encoder
{
    
    class Encoding
    {
        const int kMaxTagSymbolBitLength = 32;
        const int kMaxRawEncodingBitLength = 18;
        private const int Tagged = 0;
        private const int Raw = 1;

        public static void ConvertSignedIntsToSymbols(Span<int> input, int num, Span<int> output)
        {
            // Convert the quantized values into a format more suitable for entropy
            // encoding.
            // Put the sign bit into LSB pos and shift the rest one bit left.
            for (int i = 0; i < num; ++i)
            {
                int val = input[i];
                bool isNegative = (val < 0);
                if (isNegative)
                    val = -val - 1; // Map -1 to 0, -2 to -1, etc..
                val <<= 1;
                if (isNegative)
                    val |= 1;
                output[i] = val;
            }
        }

        /// <summary>
        /// Computes bit lengths of the input values. If numComponents > 1, the values
        /// are processed in "numComponents" sized chunks and the bit length is always
        /// computed for the largest value from the chunk.
        /// </summary>
        public static int[] ComputeBitLengths(Span<int> symbols,
            int numComponents,
            out int outMaxValue)
        {
            var outBitLengths = new int[symbols.Length / numComponents];
            int p = 0;
            outMaxValue = 0;
            // Maximum integer value across all components.
            for (int i = 0; i < symbols.Length; i += numComponents)
            {
                // Get the maximum value for a given entry across all attribute components.
                int maxComponentValue = symbols[i];
                for (int j = 1; j < numComponents; ++j)
                {
                    if (maxComponentValue < symbols[i + j])
                        maxComponentValue = symbols[i + j];
                }
                int valueMsbPos = 0;
                if (maxComponentValue > 0)
                {
                    valueMsbPos = DracoUtils.MostSignificantBit((uint) maxComponentValue);
                }
                if (maxComponentValue > outMaxValue)
                {
                    outMaxValue = maxComponentValue;
                }
                outBitLengths[p++] = valueMsbPos + 1;
            }
            return outBitLengths;
        }

        static long ComputeShannonEntropy(Span<int> symbols, int num_symbols,
            int max_value, out int out_num_unique_symbols)
        {
            // First find frequency of all unique symbols in the input array.
            int num_unique_symbols = 0;
            int[] symbol_frequencies = new int[max_value + 1];
            for (int i = 0; i < num_symbols; ++i)
            {
                ++symbol_frequencies[symbols[i]];
            }

            double total_bits = 0;
            double num_symbols_d = num_symbols;
            double log2 = Math.Log(2);
            for (int i = 0; i < max_value + 1; ++i)
            {
                if (symbol_frequencies[i] > 0)
                {
                    ++num_unique_symbols;
                    // Compute Shannon entropy for the symbol.
                    // We don't want to use std::log2 here for Android build.
                    total_bits +=
                        symbol_frequencies[i] * Math.Log(1.0 * symbol_frequencies[i] / num_symbols_d) / log2;
                }
            }

            out_num_unique_symbols = num_unique_symbols;
            // Entropy is always negative.
            return (long)-total_bits;
        }

        // Compute approximate frequency table size needed for storing the provided
        // symbols.
        static long ApproximateRAnsFrequencyTableBits(int max_value, int num_unique_symbols)
        {
            // Approximate number of bits for storing zero frequency entries using the
            // run length encoding (with max length of 64).
            long table_zero_frequency_bits =
                8 * (num_unique_symbols + (max_value - num_unique_symbols) / 64);
            return 8 * num_unique_symbols + table_zero_frequency_bits;
        }

        static long ApproximateTaggedSchemeBits(Span<int> bit_lengths, int num_components)
        {
            // Compute the total bit length used by all values (the length of data encode
            // after tags).
            ulong total_bit_length = 0;
            for (int i = 0; i < bit_lengths.Length; ++i)
            {
                total_bit_length += (ulong)bit_lengths[i];
            }

            // Compute the number of entropy bits for tags.
            int num_unique_symbols;
            long tag_bits = ComputeShannonEntropy(bit_lengths, bit_lengths.Length, 32, out num_unique_symbols);
            long tag_table_bits = ApproximateRAnsFrequencyTableBits(num_unique_symbols, num_unique_symbols);
            return tag_bits + tag_table_bits + (long)total_bit_length * num_components;
        }

        static long ApproximateRawSchemeBits(Span<int> symbols,
            int num_symbols, uint max_value,
            out int out_num_unique_symbols)
        {
            int num_unique_symbols;
            long data_bits = ComputeShannonEntropy(symbols, num_symbols, (int)max_value, out num_unique_symbols);
            long table_bits =
                ApproximateRAnsFrequencyTableBits((int)max_value, num_unique_symbols);
            out_num_unique_symbols = num_unique_symbols;
            return table_bits + data_bits;
        }

        public static void EncodeSymbols(Span<int> symbols, int numValues, int numComponents, DracoEncodeOptions options,
            EncoderBuffer targetBuffer)
        {
            if (symbols.Length == 0)
                return ;
            if (numComponents == 0)
                numComponents = 1;
            int maxValue;
            int[] bitLengths = ComputeBitLengths(symbols, numComponents, out maxValue);

            // Approximate number of bits needed for storing the symbols using the tagged
            // scheme.
            long tagged_scheme_total_bits =
                ApproximateTaggedSchemeBits(bitLengths.AsSpan(), numComponents);

            // Approximate number of bits needed for storing the symbols using the raw
            // scheme.
            int num_unique_symbols = 0;
            long raw_scheme_total_bits = ApproximateRawSchemeBits(symbols, numValues, (uint)maxValue, out num_unique_symbols);

            // The maximum bit length of a single entry value that we can encode using
            // the raw scheme.
            int max_value_bit_length =
                DracoUtils.MostSignificantBit((uint)Math.Max(1, maxValue)) + 1;

            int method;
            if (options != null && options.SymbolEncodingMethod.HasValue)
            {
                method = options.SymbolEncodingMethod.Value;
            }
            else
            {
                if (tagged_scheme_total_bits < raw_scheme_total_bits ||
                    max_value_bit_length > kMaxRawEncodingBitLength)
                {
                    method = Tagged;
                }
                else
                {
                    method = Raw;
                }
            }

            // Use the tagged scheme.
            targetBuffer.Encode((byte) method);
            if (method == Tagged)
            {
                EncodeTaggedSymbols(
                    symbols, /*numValues, */numComponents, bitLengths, targetBuffer);
            }
            else if (method == Raw)
            {
                EncodeRawSymbols(symbols, numValues, (uint)maxValue,
                    num_unique_symbols, options,
                    targetBuffer);
            }
            else
                // Unknown method selected.
                throw DracoUtils.Failed();
        }

        static bool EncodeTaggedSymbols(Span<int> symbols, int numComponents, int[] bitLengths,
            EncoderBuffer targetBuffer)
        {
            // Create entries for entropy coding. Each entry corresponds to a different
            // number of bits that are necessary to encode a given value. Every value
            // has at most 32 bits. Therefore, we need 32 different entries (for
            // bitLengts [1-32]). For each entry we compute the frequency of a given
            // bit-length in our data set.
            ulong[] frequencies = new ulong[kMaxTagSymbolBitLength];

            // Compute the frequencies from input data.
            // Maximum integer value for the values across all components.
            for (int i = 0; i < bitLengths.Length; ++i)
            {
                // Update the frequency of the associated entry id.
                ++frequencies[bitLengths[i]];
            }

            // Create one extra buffer to store raw value.
            EncoderBuffer valueBuffer = new EncoderBuffer();
            // Number of expected bits we need to store the values (can be optimized if
            // needed).
            int valueBits = kMaxTagSymbolBitLength * symbols.Length;

            // Create encoder for encoding the bit tags.
            var tagEncoder = new RAnsSymbolEncoder(5, frequencies, targetBuffer);

            // Start encoding bit tags.
            tagEncoder.StartEncoding(targetBuffer);

            // Also start encoding the values.
            valueBuffer.StartBitEncoding(valueBits, false);

            // Encoder needs the values to be encoded in the reverse order.
            for (int i = symbols.Length - numComponents; i >= 0; i -= numComponents)
            {
                int bitLength = bitLengths[i / numComponents];
                tagEncoder.EncodeSymbol(bitLength);

                // Values are always encoded in the normal order
                int j = symbols.Length - numComponents - i;
                int valueBitLength = bitLengths[j / numComponents];
                for (int c = 0; c < numComponents; ++c)
                {
                    valueBuffer.EncodeLeastSignificantBits32(valueBitLength, (uint)symbols[j + c]);
                }
            }
            tagEncoder.EndEncoding(targetBuffer);
            valueBuffer.EndBitEncoding();

            // Append the values to the end of the target buffer.
            targetBuffer.Encode(valueBuffer.Data, valueBuffer.Bytes);
            return true;
        }

        static bool EncodeRawSymbols(Span<int> symbols, int num_values,
            uint max_entry_value, int num_unique_symbols,
            DracoEncodeOptions options, EncoderBuffer target_buffer)
        {
            int symbol_bits = 0;
            if (num_unique_symbols > 0)
            {
                symbol_bits = DracoUtils.MostSignificantBit((uint)num_unique_symbols);
            }

            int unique_symbols_bit_length = symbol_bits + 1;
            // Currently, we don't support encoding of more than 2^18 unique symbols.
            if (unique_symbols_bit_length > kMaxRawEncodingBitLength)
                return false;
            int compression_level = options.GetCompressionLevel();

            // Adjust the bit_length based on compression level. Lower compression levels
            // will use fewer bits while higher compression levels use more bits. Note
            // that this is going to work for all valid bit_lengths because the actual
            // number of bits allocated for rANS encoding is hard coded as:
            // std::max(12, 3 * bit_length / 2) , therefore there will be always a
            // sufficient number of bits available for all symbols.
            // See ComputeRAnsPrecisionFromUniqueSymbolsBitLength() for the formula.
            // This hardcoded equation cannot be changed without changing the bitstream.
            if (compression_level < 4)
            {
                unique_symbols_bit_length -= 2;
            }
            else if (compression_level < 6)
            {
                unique_symbols_bit_length -= 1;
            }
            else if (compression_level > 9)
            {
                unique_symbols_bit_length += 2;
            }
            else if (compression_level > 7)
            {
                unique_symbols_bit_length += 1;
            }

            // Clamp the bit_length to a valid range.
            unique_symbols_bit_length = Math.Min(Math.Max(1, unique_symbols_bit_length),
                kMaxRawEncodingBitLength);
            target_buffer.Encode((byte) unique_symbols_bit_length);
            // Use appropriate symbol encoder based on the maximum symbol bit length.

            return EncodeRawSymbolsInternal(unique_symbols_bit_length, symbols, num_values, max_entry_value,
                target_buffer);
        }

        static bool EncodeRawSymbolsInternal(int unique_symbols_bit_length, Span<int> symbols, int num_values, uint max_entry_value, EncoderBuffer target_buffer)
        {
            // Count the frequency of each entry value.
            ulong[] frequencies = new ulong[max_entry_value + 1];
            for (int i = 0; i < num_values; ++i)
            {
                ++frequencies[symbols[i]];
            }

            RAnsSymbolEncoder encoder = new RAnsSymbolEncoder(unique_symbols_bit_length, frequencies, target_buffer);
            encoder.StartEncoding(target_buffer);
            // Encode all values.
            const bool needsReverseEncoding = true;
#if true
            if (needsReverseEncoding)
            {
                for (int i = num_values - 1; i >= 0; --i)
                {
                    encoder.EncodeSymbol(symbols[i]);
                }
            }
#else
            else
            {
                for (int i = 0; i < num_values; ++i)
                {
                    encoder.EncodeSymbol(symbols[i]);
                }
            }
#endif

            encoder.EndEncoding(target_buffer);
            return true;
        }

        static bool EncodeRawSymbols(Span<int> symbols, int maxValue, EncoderBuffer targetBuffer)
        {
            int maxEntryValue = maxValue;
            // If the maxValue is not provided, find it.
            int maxValueBits = 0;
            if (maxEntryValue > 0)
            {
                maxValueBits = DracoUtils.MostSignificantBit((uint)maxValue);
            }
            int maxValueBitLength = maxValueBits + 1;
            // Currently, we don't support encoding of values larger than 2^20.
            if (maxValueBitLength > kMaxRawEncodingBitLength)
                return false;
            targetBuffer.Encode((byte) (maxValueBitLength));
            // Use appropriate symbol encoder based on the maximum symbol bit length.

            // Count the frequency of each entry value.
            ulong[] frequencies = new ulong[maxEntryValue + 1];
            for (int i = 0; i < symbols.Length; ++i)
            {
                ++frequencies[symbols[i]];
            }

            RAnsSymbolEncoder encoder = new RAnsSymbolEncoder(maxValueBitLength, frequencies, targetBuffer);

            encoder.StartEncoding(targetBuffer);
            // Encode all values.
            for (int i = symbols.Length - 1; i >= 0; --i)
            {
                encoder.EncodeSymbol(symbols[i]);
            }
            encoder.EndEncoding(targetBuffer);
            return true;
        }

        public static bool EncodeVarint(ulong val, EncoderBuffer buffer)
        {

            // Coding of unsigned values.
            // 0-6 bit - data
            // 7 bit - next byte?
            byte out_ = 0;
            out_ |= (byte) (val & ((1 << 7) - 1));
            if (val >= (1 << 7))
            {
                out_ |= (byte) (1 << 7);
                if (!buffer.Encode(out_))
                    return false;
                if (!EncodeVarint(val >> 7, buffer))
                    return false;
                return true;
            }

            if (!buffer.Encode(out_))
                return false;
            return true;
        }

        public static bool EncodeVarint(int val, EncoderBuffer buffer)
        {
            return EncodeVarint((uint) val, buffer);
        }
        public static bool EncodeVarint(uint val, EncoderBuffer buffer)
        {

            // Coding of unsigned values.
            // 0-6 bit - data
            // 7 bit - next byte?
            byte out_ = 0;
            out_ |= (byte) (val & ((1 << 7) - 1));
            if (val >= (1 << 7))
            {
                out_ |= (byte) (1 << 7);
                if (!buffer.Encode(out_))
                    return false;
                if (!EncodeVarint(val >> 7, buffer))
                    return false;
                return true;
            }

            if (!buffer.Encode(out_))
                return false;
            return true;
        }
    }
}
