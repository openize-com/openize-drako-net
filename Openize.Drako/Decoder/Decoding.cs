using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Decoder
{
    static class Decoding
    {

        internal static void DecodeSymbols(int numValues, int numComponents,
            DecoderBuffer srcBuffer, Span<int> outValues)
        {
            if (numValues < 0)
                throw DracoUtils.Failed();
            if (numValues == 0)
                return;
            // Decode which scheme to use.
            byte scheme = srcBuffer.DecodeU8();
            if (scheme == 0)
            {
                DecodeTaggedSymbols(numValues, numComponents, srcBuffer, outValues);
            }
            else if (scheme == 1)
            {
                DecodeRawSymbols(numValues, srcBuffer, outValues);
            }
            else
                throw DracoUtils.Failed();
        }

        static bool DecodeTaggedSymbols(int numValues, int numComponents,
            DecoderBuffer srcBuffer, Span<int> outValues)
        {
            // Decode the encoded data.
            RAnsSymbolDecoder tagDecoder = new RAnsSymbolDecoder(5);
            tagDecoder.Create(srcBuffer);

            tagDecoder.StartDecoding(srcBuffer);

            if (numValues > 0 && tagDecoder.NumSymbols == 0)
                throw DracoUtils.Failed(); // Wrong number of symbols.

            // srcBuffer now points behind the encoded tag data (to the place where the
            // values are encoded).
            long tmp = srcBuffer.StartBitDecoding(false);
            int valueId = 0;
            for (int i = 0; i < numValues; i += numComponents)
            {
                // Decode the tag.
                int bitLength = tagDecoder.DecodeSymbol();
                // Decode the actual value.
                for (int j = 0; j < numComponents; ++j)
                {
                    uint val;
                    if (!srcBuffer.DecodeLeastSignificantBits32(bitLength, out val))
                        throw DracoUtils.Failed();
                    outValues[valueId++] = (int) val;
                }
            }
            tagDecoder.EndDecoding();
            srcBuffer.EndBitDecoding();
            return true;
        }

        static void DecodeRawSymbols(int numValues, DecoderBuffer srcBuffer,
            Span<int> outValues)
        {
            byte maxBitLength = srcBuffer.DecodeU8();

            RAnsSymbolDecoder decoder = new RAnsSymbolDecoder(maxBitLength);
            decoder.Create(srcBuffer);

            if (numValues > 0 && decoder.NumSymbols == 0)
                throw DracoUtils.Failed(); // Wrong number of symbols.

            decoder.StartDecoding(srcBuffer);
            for (int i = 0; i < numValues; ++i)
            {
                // Decode a symbol into the value.
                int value = decoder.DecodeSymbol();
                outValues[i] = value;
            }
            decoder.EndDecoding();
        }

        public static void ConvertSymbolsToSignedInts(Span<int> symbols, Span<int> result)
        {
            for (int i = 0; i < symbols.Length; ++i)
            {
                int val = symbols[i];
                bool isNegative = (val & 1) != 0;//lowest bit is sign bit
                val >>= 1;
                int ret = val;
                if (isNegative)
                    ret = -ret - 1;
                result[i] = ret;
            }
        }

        /// <summary>
        /// Decodes a specified integer as varint. Note that the IntTypeT must be the
        /// same as the one used in the corresponding EncodeVarint() call.
        /// </summary>
        /// <param name="out_val"></param>
        /// <param name="buffer"></param>
        /// <returns></returns>
        public static uint DecodeVarintU32(this DecoderBuffer buffer)
        {
            // Coding of unsigned values.
            // 0-6 bit - data
            // 7 bit - next byte?
            byte in_ = buffer.DecodeU8();
            uint out_val;

            if ((in_ & (1 << 7)) != 0)
            {
                // Next byte is available, decode it first.
                out_val = DecodeVarintU32(buffer);
                // Append decoded info from this byte.
                out_val <<= 7;
                out_val |= (uint)(in_ & ((1 << 7) - 1));
            }
            else
            {
                // Last byte reached
                out_val = in_;
            }
            return out_val;
        }
        public static ushort DecodeVarintU16(this DecoderBuffer buffer)
        {
            // Coding of unsigned values.
            // 0-6 bit - data
            // 7 bit - next byte?
            byte in_ = buffer.DecodeU8();
            ushort out_val;

            if ((in_ & (1 << 7)) != 0)
            {
                // Next byte is available, decode it first.
                out_val = DecodeVarintU16(buffer);
                // Append decoded info from this byte.
                out_val <<= 7;
                out_val |= (ushort)(in_ & ((1 << 7) - 1));
            }
            else
            {
                // Last byte reached
                out_val = in_;
            }
            return out_val;
        }
        public static ulong DecodeVarintU64(this DecoderBuffer buffer)
        {
            // Coding of unsigned values.
            // 0-6 bit - data
            // 7 bit - next byte?
            byte in_ = buffer.DecodeU8();
            ulong out_val;

            if ((in_ & (1 << 7)) != 0)
            {
                // Next byte is available, decode it first.

                out_val = DecodeVarintU64(buffer);
                // Append decoded info from this byte.
                out_val <<= 7;
                out_val |= (uint)(in_ & ((1 << 7) - 1));
            }
            else
            {
                // Last byte reached
                out_val = in_;
            }
            return out_val;
        }
    }
}
