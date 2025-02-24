using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Encoder
{
    class DirectBitEncoder : IBitEncoder
    {

        // Must be called before any Encode* function is called.
        public void StartEncoding()
        {
            Clear();
        }

        // Encode one bit. If |bit| is true encode a 1, otherwise encode a 0.
        public void EncodeBit(bool bit)
        {
            if (bit)
            {
                local_bits_ |= 1U << (31 - num_local_bits_);
            }

            num_local_bits_++;
            if (num_local_bits_ == 32)
            {
                bits_.Add((int)local_bits_);
                num_local_bits_ = 0;
                local_bits_ = 0;
            }
        }

        // Encode |nbits| of |value|, starting from the least significant bit.
        // |nbits| must be > 0 and <= 32.
        public void EncodeLeastSignificantBits32(int nbits, uint value)
        {
            //DRACO_DCHECK_EQ(true, nbits <= 32);
            //DRACO_DCHECK_EQ(true, nbits > 0);

            int remaining = 32 - num_local_bits_;

            // Make sure there are no leading bits that should not be encoded and
            // start from here.
            value = value << (32 - nbits);
            if (nbits <= remaining)
            {
                value = value >> num_local_bits_;
                local_bits_ = local_bits_ | value;
                num_local_bits_ += nbits;
                if (num_local_bits_ == 32)
                {
                    bits_.Add((int)local_bits_);
                    local_bits_ = 0;
                    num_local_bits_ = 0;
                }
            }
            else
            {
                value = value >> (32 - nbits);
                num_local_bits_ = nbits - remaining;
                uint value_l = value >> num_local_bits_;
                local_bits_ = local_bits_ | value_l;
                bits_.Add((int)local_bits_);
                local_bits_ = value << (32 - num_local_bits_);
            }
        }

        // Ends the bit encoding and stores the result into the target_buffer.
        public void EndEncoding(EncoderBuffer target_buffer)
        {

            bits_.Add((int)local_bits_);
            var size_in_byte = bits_.Count * 4;
            target_buffer.Encode(size_in_byte);
            target_buffer.Encode(bits_);
            Clear();
        }

        public void Clear()
        {

            bits_.Clear();
            local_bits_ = 0;
            num_local_bits_ = 0;
        }

        private IntList bits_ = new IntList();
        uint local_bits_;
        int num_local_bits_;
    }
}
