using Openize.Drako.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako.Encoder
{
    class DynamicIntegerPointsKdTreeEncoder
    {

        struct EncodingStatus
        {
            public EncodingStatus(int begin, int end, int last_axis_, int stack_pos_)
            {
                this.begin = begin;
                this.end = end;
                this.last_axis = last_axis_;
                this.stack_pos = stack_pos_;
                num_remaining_points = end - begin;
            }

            public int begin;
            public int end;
            public int last_axis;
            public int num_remaining_points;
            public int stack_pos; // used to get base and levels
        };

        private int compression_level_t;
          uint bit_length_;
          int num_points_;
          int dimension_;
          IBitEncoder numbers_encoder_;
          private IBitEncoder remaining_bits_encoder_ = new DirectBitEncoder();
          private IBitEncoder axis_encoder_ = new DirectBitEncoder();
          private IBitEncoder half_encoder_ = new DirectBitEncoder();
          uint[] deviations_;
          uint[] num_remaining_bits_;
          int[] axes_;
          uint[][] base_stack_;
          uint[][] levels_stack_;
        public DynamicIntegerPointsKdTreeEncoder(int compression_level, int dimension)
        {
            this.compression_level_t = compression_level;
            this.dimension_ = dimension;
            this.deviations_ = new uint[dimension];
            this.num_remaining_bits_ = new uint[dimension];
            this.axes_ = new int[dimension];
            this.base_stack_ = new uint[32 * dimension + 1][];
            for (int i = 0; i < base_stack_.Length; i++)
                base_stack_[i] = new uint[dimension];
            this.levels_stack_ = new uint[32 * dimension + 1][];
            for (int i = 0; i < levels_stack_.Length; i++)
                levels_stack_[i] = new uint[dimension];

            numbers_encoder_ = CreateNumbersEncoder();
        }

        private IBitEncoder CreateNumbersEncoder()
        {
            switch (compression_level_t)
            {
                case 0:
                case 1:
                    return new DirectBitEncoder();
                case 2:
                case 3:
                    return new RAnsBitEncoder();
                case 4:
                case 5:
                case 6:
                    return new FoldedBit32Encoder();
                default:
                    throw new InvalidOperationException("Invalid compression level");
            }
        }

        public void EncodePoints(int[][] array, int bit_length, EncoderBuffer buffer)
        {
            bit_length_ = (uint)bit_length;
            num_points_ = array.Length;

            buffer.Encode(bit_length_);
            buffer.Encode(num_points_);
            if (num_points_ == 0)
                return ;

            numbers_encoder_.StartEncoding();
            remaining_bits_encoder_.StartEncoding();
            axis_encoder_.StartEncoding();
            half_encoder_.StartEncoding();

            EncodeInternal(array);

            numbers_encoder_.EndEncoding(buffer);
            remaining_bits_encoder_.EndEncoding(buffer);
            axis_encoder_.EndEncoding(buffer);
            half_encoder_.EndEncoding(buffer);

        }

        public void EncodeNumber(int nbits, uint value)
        {
            numbers_encoder_.EncodeLeastSignificantBits32(nbits, value);
        }

        public int GetAndEncodeAxis(int[][] array, int begin, int end, uint[] old_base, uint[] levels, int last_axis)
        {
            bool select_axis = compression_level_t == 6;
            if (!select_axis)
                return DracoUtils.IncrementMod(last_axis, dimension_);

            // For many points this function selects the axis that should be used
            // for the split by keeping as many points as possible bundled.
            // In the best case we do not split the point cloud at all.
            // For lower number of points, we simply choose the axis that is refined the
            // least so far.

            //DRACO_DCHECK_EQ(true, end - begin != 0);

            int best_axis = 0;
            if (array.Length < 64)
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
                uint size = (uint)array.Length;
                for (int i = 0; i < dimension_; i++)
                {
                    deviations_[i] = 0;
                    num_remaining_bits_[i] = bit_length_ - levels[i];
                    if (num_remaining_bits_[i] > 0)
                    {
                        var split = (int)(old_base[i] + (1 << (int)(num_remaining_bits_[i] - 1)));
                        uint deviation = 0;
                        for (int it = 0; it < size; it++)
                        {
                            //deviation += (array[it][i] < split) ? 1U : 0U;
                            if (array[it][i] < split)
                                deviation++;
                        }

                        deviations_[i] = Math.Max(size - deviation, deviation);
                    }
                }

                uint max_value = 0;
                best_axis = 0;
                for (int i = 0; i < dimension_; i++)
                {
                    // If axis can be subdivided.
                    if (num_remaining_bits_[i] != 0)
                    {
                        // Check if this is the better axis.
                        if (max_value < deviations_[i])
                        {
                            max_value = deviations_[i];
                            best_axis = i;
                        }
                    }
                }

                axis_encoder_.EncodeLeastSignificantBits32(4, (uint)best_axis);
            }
            return best_axis;
        }

        void EncodeInternal(int[][] array)
        {

            base_stack_[0] = new uint[dimension_];
            levels_stack_[0] = new uint[dimension_];
            var init_status = new EncodingStatus(0, array.Length, 0, 0);
            Stack<EncodingStatus> status_stack = new Stack<EncodingStatus>();
            status_stack.Push(init_status);

            // TODO(hemmer): use preallocated vector instead of stack.
            while (status_stack.Count > 0)
            {
                var status = status_stack.Peek();
                status_stack.Pop();

                var begin = status.begin;
                var end = status.end;
                int last_axis = status.last_axis;
                int stack_pos = status.stack_pos;
                var old_base = base_stack_[stack_pos];
                var levels = levels_stack_[stack_pos];

                int axis = GetAndEncodeAxis(array, begin, end, old_base, levels, last_axis);
                uint level = levels[axis];
                int num_remaining_points = end - begin;

                // If this happens all axis are subdivided to the end.
                if ((bit_length_ - level) == 0)
                    continue;

                // Fast encoding of remaining bits if number of points is 1 or 2.
                // Doing this also for 2 gives a slight additional speed up.
                if (num_remaining_points <= 2)
                {
                    // TODO(hemmer): axes_ not necessary, remove would change bitstream!
                    axes_[0] = axis;
                    for (uint i = 1; i < dimension_; i++)
                    {
                        axes_[i] = DracoUtils.IncrementMod(axes_[i - 1], dimension_);
                    }

                    for (uint i = 0; i < num_remaining_points; ++i)
                    {
                        var p = array[begin + i];
                        for (uint j = 0; j < dimension_; j++)
                        {
                            var num_remaining_bits = (int)(bit_length_ - levels[axes_[j]]);
                            if (num_remaining_bits != 0)
                            {
                                remaining_bits_encoder_.EncodeLeastSignificantBits32(
                                    num_remaining_bits, (uint)p[axes_[j]]);
                            }
                        }
                    }

                    continue;
                }

                uint modifier = 1U << (int)((bit_length_ - level) - 1);
                Array.Copy(old_base, base_stack_[stack_pos + 1], old_base.Length); // copy
                base_stack_[stack_pos + 1][axis] += modifier;
                uint[] new_base = base_stack_[stack_pos + 1];

                int split = Partition(array, begin, end, axis, new_base[axis]);

                //DRACO_DCHECK_EQ(true, (end - begin) > 0);

                // Encode number of points in first and second half.
                int required_bits = DracoUtils.MostSignificantBit((uint)num_remaining_points);

                int first_half = split - begin;
                int second_half = end - split;
                bool left = first_half < second_half;

                if (first_half != second_half)
                    half_encoder_.EncodeBit(left);

                if (left)
                {
                    EncodeNumber(required_bits, (uint)(num_remaining_points / 2 - first_half));
                }
                else
                {
                    EncodeNumber(required_bits, (uint)(num_remaining_points / 2 - second_half));
                }

                levels_stack_[stack_pos][axis] += 1;
                Array.Copy(levels_stack_[stack_pos], levels_stack_[stack_pos + 1], dimension_) ; // copy
                if (split != begin)
                    status_stack.Push(new EncodingStatus(begin, split, axis, stack_pos));
                if (split != end)
                    status_stack.Push(new EncodingStatus(split, end, axis, stack_pos + 1));
            }
        }

        private int Partition(int[][] points, int first, int last, int axis, uint v)
        {

            for (;;)
            {
                for (;;)
                {
                    if (first == last)
                        return first;
                    if(!(points[first][axis] < v))
                        break;
                    ++first;
                }

                do
                {
                    --last;
                    if (first == last)
                        return first;
                } while (!(points[last][axis] < v));

                var t = points[first];
                points[first] = points[last];
                points[last] = t;
                ++first;
            }
        }
    }
}
