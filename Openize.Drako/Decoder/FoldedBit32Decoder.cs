using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Decoder
{
    class FoldedBit32Decoder : IBitDecoder
    {
        RAnsBitDecoder[] folded_number_decoders_ = new RAnsBitDecoder[32];
        RAnsBitDecoder bit_decoder_;

        public FoldedBit32Decoder()
        {
            for (int i = 0; i < folded_number_decoders_.Length; i++)
                folded_number_decoders_[i] = new RAnsBitDecoder();
            bit_decoder_ = new RAnsBitDecoder();
        }


        // Sets |source_buffer| as the buffer to decode bits from.
        public void StartDecoding(DecoderBuffer source_buffer)
        {
            for (int i = 0; i < 32; i++)
            {
                folded_number_decoders_[i].StartDecoding(source_buffer);
            }

            bit_decoder_.StartDecoding(source_buffer);
        }

        // Decode one bit. Returns true if the bit is a 1, otherwise false.
        public bool DecodeNextBit()
        {
            return bit_decoder_.DecodeNextBit();
        }

        // Decode the next |nbits| and return the sequence in |value|. |nbits| must be
        // > 0 and <= 32.
        public uint DecodeLeastSignificantBits32(int nbits)
        {
            uint result = 0;
            for (int i = 0; i < nbits; ++i)
            {
                bool bit = folded_number_decoders_[i].DecodeNextBit();
                result = (result << 1) + (bit ? 1U : 0U);
            }

            return result;
        }

        public void EndDecoding()
        {
            for (int i = 0; i < 32; i++)
            {
                folded_number_decoders_[i].EndDecoding();
            }

            bit_decoder_.EndDecoding();
        }

        public void Clear()
        {
            for (int i = 0; i < 32; i++)
            {
                folded_number_decoders_[i].Clear();
            }

            bit_decoder_.Clear();
        }
    }
}
