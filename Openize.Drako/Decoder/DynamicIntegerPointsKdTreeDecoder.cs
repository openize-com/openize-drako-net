using Openize.Drako.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Openize.Drako.Decoder
{
    class DynamicIntegerPointsKdTreeDecoder
    {


        class Status
        {
            public Status(int num_remaining_points_, int last_axis_,
                int stack_pos_)
            {
                this.num_remaining_points = num_remaining_points_;
                this.last_axis = last_axis_;
                this.stack_pos = stack_pos_;
            }

            public int num_remaining_points;
            public int last_axis;
            public int stack_pos; // used to get base and levels
        }

        int bit_length_;
        int num_points_;
        int num_decoded_points_;
        int dimension_;
        IBitDecoder numbers_decoder_;
        DirectBitDecoder remaining_bits_decoder_;
        DirectBitDecoder axis_decoder_;
        DirectBitDecoder half_decoder_;
        int[] p_;
        int[] axes_;
        int[][] base_stack_;
        int[][] levels_stack_;
        private int compression_level_t;

        public DynamicIntegerPointsKdTreeDecoder(int compression_level, int dimension)
        {
            this.compression_level_t = compression_level;
            dimension_ = dimension;

            p_ = new int[dimension];
            axes_ = new int[dimension];
            // Init the stack with the maximum depth of the tree.
            // +1 for a second leaf.
            base_stack_ = new int[32 * dimension + 1][];
            levels_stack_ = new int[32 * dimension + 1][];
            for (int i = 0; i < base_stack_.Length; i++)
            {
                base_stack_[i] = new int[dimension];
                levels_stack_[i] = new int[dimension];
            }

            switch (compression_level)
            {
                case 0:
                case 1:
                    numbers_decoder_ = new DirectBitDecoder();
                    break;
                case 2:
                case 3:
                    numbers_decoder_ = new RAnsBitDecoder();
                    break;
                case 4:
                case 5:
                case 6:
                    numbers_decoder_ = new FoldedBit32Decoder();
                    break;
                default:
                    throw new InvalidOperationException("Invalid compression level.");
            }
            remaining_bits_decoder_ = new DirectBitDecoder();
            axis_decoder_ = new DirectBitDecoder();
            half_decoder_ = new DirectBitDecoder();
        }

        public void DecodePoints(DecoderBuffer buffer, PointAttributeVectorOutputIterator oit)
        {
            bit_length_ = buffer.DecodeI32();
            if (bit_length_ > 32)
                throw DracoUtils.Failed();
            num_points_ = buffer.DecodeI32(); 
            if (num_points_ == 0)
                return;
            num_decoded_points_ = 0;

            numbers_decoder_.StartDecoding(buffer);
            remaining_bits_decoder_.StartDecoding(buffer);
            axis_decoder_.StartDecoding(buffer);
            half_decoder_.StartDecoding(buffer);

            DecodeInternal(num_points_, oit);

            numbers_decoder_.EndDecoding();
            remaining_bits_decoder_.EndDecoding();
            axis_decoder_.EndDecoding();
            half_decoder_.EndDecoding();

        }

        int GetAxis(int num_remaining_points, int[] levels, int last_axis)
        {
            var select_axis = compression_level_t == 6;
            if (!select_axis)
                return DracoUtils.IncrementMod(last_axis, dimension_);

            int best_axis = 0;
            if (num_remaining_points < 64)
            {
                for (int axis = 1; axis < dimension_; ++axis)
                {
                    if (levels[best_axis] > levels[axis])
                    {
                        best_axis = axis;
                    }
                }
            }
            else
            {
                best_axis = (int) axis_decoder_.DecodeLeastSignificantBits32(4);
            }

            return best_axis;
        }

        void DecodeInternal(int num_points, PointAttributeVectorOutputIterator oit)
        {
            base_stack_[0] = new int[dimension_];
            levels_stack_[0] = new int [dimension_];
            Status init_status = new Status(num_points, 0, 0);
            Stack<Status> status_stack = new Stack<Status>();
            status_stack.Push(init_status);

            // TODO(hemmer): use preallocated vector instead of stack.
            while (status_stack.Count > 0)
            {
                var status = status_stack.Peek();
                status_stack.Pop();

                int num_remaining_points = status.num_remaining_points;
                int last_axis = status.last_axis;
                int stack_pos = status.stack_pos;
                var old_base = base_stack_[stack_pos];
                var levels = levels_stack_[stack_pos];

                if (num_remaining_points > num_points)
                    throw DracoUtils.Failed();

                int axis = GetAxis(num_remaining_points, levels, last_axis);
                if (axis >= dimension_)
                    throw DracoUtils.Failed();

                int level = levels[axis];

                // All axes have been fully subdivided, just output points.
                if ((bit_length_ - level) == 0)
                {
                    for (int i = 0; i < num_remaining_points; i++)
                    {
                        oit.Set(old_base);
                        oit.Next();
                        ++num_decoded_points_;
                    }

                    continue;
                }

                //DRACO_DCHECK_EQ(true, num_remaining_points != 0);

                // Fast decoding of remaining bits if number of points is 1 or 2.
                int num_remaining_bits;
                if (num_remaining_points <= 2)
                {
                    // TODO(hemmer): axes_ not necessary, remove would change bitstream!
                    axes_[0] = axis;
                    for (int i = 1; i < dimension_; i++)
                    {
                        axes_[i] = DracoUtils.IncrementMod(axes_[i - 1], dimension_);
                    }

                    for (int i = 0; i < num_remaining_points; ++i)
                    {
                        for (int j = 0; j < dimension_; j++)
                        {
                            p_[axes_[j]] = 0;
                            num_remaining_bits = bit_length_ - levels[axes_[j]];
                            if (num_remaining_bits != 0)
                                p_[axes_[j]] =
                                    (int) remaining_bits_decoder_.DecodeLeastSignificantBits32(num_remaining_bits);
                            p_[axes_[j]] = old_base[axes_[j]] | p_[axes_[j]];
                        }

                        oit.Set(p_);
                        oit.Next();
                        ++num_decoded_points_;
                    }

                    continue;
                }

                if (num_decoded_points_ > num_points_)
                    throw DracoUtils.Failed();

                num_remaining_bits = bit_length_ - level;
                int modifier = 1 << (num_remaining_bits - 1);
                Array.Copy(old_base, base_stack_[stack_pos + 1], dimension_); // copy
                base_stack_[stack_pos + 1][axis] += modifier; // new base

                int incoming_bits = DracoUtils.MostSignificantBit((uint)num_remaining_points);

                int number = DecodeNumber(incoming_bits);

                int first_half = num_remaining_points / 2 - number;
                int second_half = num_remaining_points - first_half;

                if (first_half != second_half)
                {
                    if (!half_decoder_.DecodeNextBit())
                    {
                        var t = first_half;
                        first_half = second_half;
                        second_half = t;
                    }
                }


                levels_stack_[stack_pos][axis] += 1;
                Array.Copy(levels_stack_[stack_pos], levels_stack_[stack_pos + 1],dimension_); // copy
                if (first_half != 0)
                    status_stack.Push(new Status(first_half, axis, stack_pos));
                if (second_half != 0)
                    status_stack.Push(new Status(second_half, axis, stack_pos + 1));
            }

        }

        private int DecodeNumber(int nbits)
        {
            return (int) numbers_decoder_.DecodeLeastSignificantBits32(nbits);
        }
    }

}
