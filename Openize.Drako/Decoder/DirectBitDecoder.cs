using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Decoder
{
    class DirectBitDecoder : IBitDecoder
    {

        uint[] bits_;
        int pos_;
        int num_used_bits_;

        public void Clear()
        {
            bits_ = null;
            num_used_bits_ = 0;
            pos_ = 0;
        }

        // Sets |source_buffer| as the buffer to decode bits from.
        public void StartDecoding(DecoderBuffer source_buffer)
        {

            Clear();
            int size_in_bytes = source_buffer.DecodeI32();

            // Check that size_in_bytes is > 0 and a multiple of 4 as the encoder always
            // encodes 32 bit elements.
            if (size_in_bytes == 0 || (size_in_bytes & 0x3) != 0)
                throw DracoUtils.Failed();
            if (size_in_bytes > source_buffer.RemainingSize)
                throw DracoUtils.Failed();
            int num_32bit_elements = size_in_bytes / 4;
            bits_ = new uint[num_32bit_elements];
            if (!source_buffer.Decode(bits_))
                throw DracoUtils.Failed();
            pos_ = 0;
            num_used_bits_ = 0;
        }

        // Decode one bit. Returns true if the bit is a 1, otherwise false.
        public bool DecodeNextBit()
        {
            int selector = 1 << (31 - num_used_bits_);
            if (pos_ == bits_.Length)
            {
                return false;
            }

            bool bit = (bits_[pos_] & selector) != 0;
            ++num_used_bits_;
            if (num_used_bits_ == 32)
            {
                ++pos_;
                num_used_bits_ = 0;
            }

            return bit;
        }

        // Decode the next |nbits| and return the sequence in |value|. |nbits| must be
        // > 0 and <= 32.
        public uint DecodeLeastSignificantBits32(int nbits)
        {
            //DRACO_DCHECK_EQ(true, nbits <= 32);
            //DRACO_DCHECK_EQ(true, nbits > 0);
            int remaining = 32 - num_used_bits_;
            uint value = 0;
            if (nbits <= remaining)
            {
                if (pos_ == bits_.Length)
                {
                    return 0;
                }

                value = (uint) ((bits_[pos_] << num_used_bits_) >> (32 - nbits));
                num_used_bits_ += nbits;
                if (num_used_bits_ == 32)
                {
                    ++pos_;
                    num_used_bits_ = 0;
                }
            }
            else
            {
                if (pos_ + 1 == bits_.Length)
                {
                    return 0;
                }

                uint value_l = ((bits_[pos_]) << num_used_bits_);
                num_used_bits_ = nbits - remaining;
                ++pos_;
                uint value_r = (bits_[pos_]) >> (32 - num_used_bits_);
                value = (value_l >> (32 - num_used_bits_ - remaining)) | value_r;
            }

            return value;
        }

        public void EndDecoding()
        {
        }
    }
}
