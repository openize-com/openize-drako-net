using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako.Encoder
{
    class FoldedBit32Encoder : IBitEncoder
    {

        public FoldedBit32Encoder()
        {
            folded_number_encoders_ = new RAnsBitEncoder[32];
            for (int i = 0; i < folded_number_encoders_.Length; i++)
                folded_number_encoders_[i] = new RAnsBitEncoder();
            bit_encoder_ = new RAnsBitEncoder();
        }

        // Must be called before any Encode* function is called.
        public void StartEncoding()
        {
            for (int i = 0; i < 32; i++)
            {
                folded_number_encoders_[i].StartEncoding();
            }

            bit_encoder_.StartEncoding();
        }

        // Encode one bit. If |bit| is true encode a 1, otherwise encode a 0.
        public void EncodeBit(bool bit)
        {
            bit_encoder_.EncodeBit(bit);
        }


        // Encode |nbits| of |value|, starting from the least significant bit.
        // |nbits| must be > 0 and <= 32.
        public void EncodeLeastSignificantBits32(int nbits, uint value)
        {
            int selector = 1 << (nbits - 1);
            for (int i = 0; i < nbits; i++)
            {
                bool bit = (value & selector) != 0;
                folded_number_encoders_[i].EncodeBit(bit);
                selector = selector >> 1;
            }
        }

        // Ends the bit encoding and stores the result into the target_buffer.
        public void EndEncoding(EncoderBuffer target_buffer)
        {
            for (int i = 0; i < 32; i++)
            {
                folded_number_encoders_[i].EndEncoding(target_buffer);
            }

            bit_encoder_.EndEncoding(target_buffer);
        }

        public void Clear()
        {
            for (int i = 0; i < 32; i++)
            {
                folded_number_encoders_[i].Clear();
            }

            bit_encoder_.Clear();
        }

        private RAnsBitEncoder[] folded_number_encoders_;
        private RAnsBitEncoder bit_encoder_;
    }
}
