using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako
{
    class RAnsBitCodec
    {
        protected const int ansP8Precision = 256;
        protected const int ansP8Shift = 8;
        protected const int ansP10Precision = 1024;

        protected const int lBase = (ansP10Precision*4); // lBase % precision must be 0
        protected const int ioBase = 256;

        public struct RansSym
        {
            public uint prob;
            public uint cumProb; // not-inclusive
            public override string ToString()
            {
                return string.Format("prob = {0}, cumProb={1}", prob, cumProb);
            }
        };

        /// <summary>
        /// Computes the desired precision of the rANS method for the specified maximal
        /// symbol bit length of the input data.
        /// </summary>
        /// <param name="maxBitLength"></param>
        /// <returns></returns>
        public static int ComputeRAnsUnclampedPrecision(int maxBitLength)
        {
            return (3 * maxBitLength) / 2;
        }

        /// <summary>
        /// Computes the desired precision clamped to guarantee a valid funcionality of
        /// our rANS library (which is between 12 to 20 bits).
        /// </summary>
        /// <param name="maxBitLength"></param>
        /// <returns></returns>
        public static int ComputeRAnsPrecisionFromMaxSymbolBitLength(int maxBitLength)
        {
            return ComputeRAnsUnclampedPrecision(maxBitLength) < 12
                ? 12
                : ComputeRAnsUnclampedPrecision(maxBitLength) > 20
                    ? 20
                    : ComputeRAnsUnclampedPrecision(maxBitLength);
        }
    }
}
