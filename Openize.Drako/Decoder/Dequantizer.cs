using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako.Decoder
{

    /// <summary>
    /// Class for dequantizing values that were previously quantized using the
    /// Quantizer class.
    /// </summary>
    class Dequantizer
    {
        private float delta;


        /// <summary>
        /// Initializes the dequantizer. Both parameters must correspond to the values
        /// provided to the initializer of the Quantizer class.
        /// </summary>
        public Dequantizer(float range, int maxQuantizedValue)
        {
            if (maxQuantizedValue > 0)
            {
                delta = range / (float) maxQuantizedValue;
            }
        }

        public float DequantizeFloat(int v)
        {
            return v * delta;
        }
    }
}
