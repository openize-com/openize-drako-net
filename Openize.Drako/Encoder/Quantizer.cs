using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Encoder
{

    /// <summary>
    /// Class for quantizing single precision floating point values. The values must
    /// be centered around zero and be within interval (-range, +range), where the
    /// range is specified in the Init() method.
    /// </summary>
    struct Quantizer
    {
        private float inverse_delta_;

        public Quantizer(float range, int maxQuantizedValue)
        {
            inverse_delta_ = (float)(maxQuantizedValue) / range;
        }

        public int QuantizeFloat(float val)
        {
            val *= inverse_delta_;
            return (int)(Math.Floor(val + 0.5f));
        }

    }
}
