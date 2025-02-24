using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako
{
    class OctahedronToolBox
    {
        int quantization_bits_ = -1;
        private int max_quantized_value_ = -1;
        private int max_value_ = -1;
        private int center_value_ = -1;


        public void SetQuantizationBits(int q)
        {
            if (q < 2 || q > 30)
                throw new ArgumentException("Invalid quantization parameters");
            quantization_bits_ = q;
            max_quantized_value_ = (1 << quantization_bits_) - 1;
            max_value_ = max_quantized_value_ - 1;
            center_value_ = max_value_ / 2;
        }

        public bool IsInitialized => quantization_bits_ != -1;

  // |s| and |t| are expected to be signed values.
        public bool IsInDiamond(int s, int t)
        {
            // Expect center already at origin.
            return Math.Abs(s) + Math.Abs(t) <= center_value_;
        }
        public void InvertDiamond(ref IntVector v)
        {
            // Expect center already at origin.
            int sign_s = 0;
            int sign_t = 0;
            if (v.x >= 0 && v.y >= 0)
            {
                sign_s = 1;
                sign_t = 1;
            }
            else if (v.x <= 0 && v.y <= 0)
            {
                sign_s = -1;
                sign_t = -1;
            }
            else
            {
                sign_s = (v.x > 0) ? 1 : -1;
                sign_t = (v.y > 0) ? 1 : -1;
            }

            int corner_point_s = sign_s * center_value_;
            int corner_point_t = sign_t * center_value_;
            v.x = 2 * v.x - corner_point_s;
            v.y = 2 * v.y - corner_point_t;
            if (sign_s * sign_t >= 0)
            {
                int temp = v.x;
                v.x = -v.y;
                v.y = -temp;
            }
            else
            {
                int temp = v.x;
                v.x = v.y;
                v.y = temp;
            }

            v.x = (v.x + corner_point_s) / 2;
            v.y = (v.y + corner_point_t) / 2;
        }

        public void InvertDiamond(ref int s, ref int t)
        {
            // Expect center already at origin.
            int sign_s = 0;
            int sign_t = 0;
            if (s >= 0 && t >= 0)
            {
                sign_s = 1;
                sign_t = 1;
            }
            else if (s <= 0 && t <= 0)
            {
                sign_s = -1;
                sign_t = -1;
            }
            else
            {
                sign_s = (s > 0) ? 1 : -1;
                sign_t = (t > 0) ? 1 : -1;
            }

            int corner_point_s = sign_s * center_value_;
            int corner_point_t = sign_t * center_value_;
            s = 2 * s - corner_point_s;
            t = 2 * t - corner_point_t;
            if (sign_s * sign_t >= 0)
            {
                int temp = s;
                s = -t;
                t = -temp;
            }
            else
            {
                int temp = s;
                s = t;
                t = temp;
            }

            s = (s + corner_point_s) / 2;
            t = (t + corner_point_t) / 2;
        }

        public void InvertDirection(ref int s, ref int t)
        {
            s *= -1;
            t *= -1;
            InvertDiamond(ref s, ref t);
        }


        // For correction values.
        public int ModMax(int x)
        {
            if (x > center_value_)
                return x - max_quantized_value_;
            if (x < -center_value_)
                return x + max_quantized_value_;
            return x;
        }

        // For correction values.
        public int MakePositive(int x)
        {
            if (x < 0)
                return x + max_quantized_value_;
            return x;
        }

        public int QuantizationBits => quantization_bits_;
        public int MaxQuantizedValue => max_quantized_value_;
        public int MaxValue => max_value_;
        public int CenterValue => center_value_;


  // Normalize |vec| such that its abs sum is equal to the center value;
        public void CanonicalizeIntegerVector(int[] vec)
        {
            long abs_sum = Math.Abs(vec[0]) +
                           Math.Abs(vec[1]) +
                           Math.Abs(vec[2]);

            if (abs_sum == 0)
            {
                vec[0] = center_value_; // vec[1] == v[2] == 0
            }
            else
            {
                vec[0] = (int) (((long)vec[0] * center_value_) / abs_sum);
                vec[1] = (int) (((long)vec[1] * center_value_) / abs_sum);
                if (vec[2] >= 0)
                {
                    vec[2] = center_value_ - Math.Abs(vec[0]) - Math.Abs(vec[1]);
                }
                else
                {
                    vec[2] = -(center_value_ - Math.Abs(vec[0]) - Math.Abs(vec[1]));
                }
            }
        }

  // Converts an integer vector to octahedral coordinates.
  // Precondition: |int_vec| abs sum must equal center value.
        public void IntegerVectorToQuantizedOctahedralCoords(int[] int_vec, out int out_s,  out int out_t) {
            int s, t;
            if (int_vec[0] >= 0)
            {
                // Right hemisphere.
                s = (int_vec[1] + center_value_);
                t = (int_vec[2] + center_value_);
            }
            else
            {
                // Left hemisphere.
                if (int_vec[1] < 0)
                {
                    s = Math.Abs(int_vec[2]);
                }
                else
                {
                    s = (max_value_ - Math.Abs(int_vec[2]));
                }

                if (int_vec[2] < 0)
                {
                    t = Math.Abs(int_vec[1]);
                }
                else
                {
                    t = (max_value_ - Math.Abs(int_vec[1]));
                }
            }

            CanonicalizeOctahedralCoords(s, t, out out_s, out out_t);
        }

  // Convert all edge points in the top left and bottom right quadrants to
  // their corresponding position in the bottom left and top right quadrants.
  // Convert all corner edge points to the top right corner.
        private void CanonicalizeOctahedralCoords(int s, int t, out int out_s, out int out_t)
        {
            if ((s == 0 && t == 0) || (s == 0 && t == max_value_) ||
                (s == max_value_ && t == 0))
            {
                s = max_value_;
                t = max_value_;
            }
            else if (s == 0 && t > center_value_)
            {
                t = center_value_ - (t - center_value_);
            }
            else if (s == max_value_ && t < center_value_)
            {
                t = center_value_ + (center_value_ - t);
            }
            else if (t == max_value_ && s < center_value_)
            {
                s = center_value_ + (center_value_ - s);
            }
            else if (t == 0 && s > center_value_)
            {
                s = center_value_ - (s - center_value_);
            }

            out_s = s;
            out_t = t;
        }

        public void FloatVectorToQuantizedOctahedralCoords(Span<float> vector, out int out_s, out int out_t)
        {
            double abs_sum = Math.Abs(vector[0]) +
                             Math.Abs(vector[1]) +
                             Math.Abs(vector[2]);

            // Adjust values such that abs sum equals 1.
            double[] scaled_vector = new double[3];
            if (abs_sum > 1e-6)
            {
                // Scale needed to project the vector to the surface of an octahedron.
                double scale = 1.0 / abs_sum;
                scaled_vector[0] = vector[0] * scale;
                scaled_vector[1] = vector[1] * scale;
                scaled_vector[2] = vector[2] * scale;
            }
            else
            {
                scaled_vector[0] = 1.0;
                scaled_vector[1] = 0;
                scaled_vector[2] = 0;
            }

            // Scale vector such that the sum equals the center value.
            int[] int_vec = new int[3];
            int_vec[0] = (int) Math.Floor(scaled_vector[0] * center_value_ + 0.5);
            int_vec[1] = (int) Math.Floor(scaled_vector[1] * center_value_ + 0.5);
            // Make sure the sum is exactly the center value.
            int_vec[2] = center_value_ - Math.Abs(int_vec[0]) - Math.Abs(int_vec[1]);
            if (int_vec[2] < 0)
            {
                // If the sum of first two coordinates is too large, we need to decrease
                // the length of one of the coordinates.
                if (int_vec[1] > 0)
                {
                    int_vec[1] += int_vec[2];
                }
                else
                {
                    int_vec[1] -= int_vec[2];
                }

                int_vec[2] = 0;
            }

            // Take care of the sign.
            if (scaled_vector[2] < 0)
                int_vec[2] *= -1;

            IntegerVectorToQuantizedOctahedralCoords(int_vec, out out_s, out out_t);
        }
    }
}
