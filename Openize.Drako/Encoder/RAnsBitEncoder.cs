using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Encoder
{

    /// <summary>
    /// Class for encoding a sequence of bits using rANS. The probability table used
    /// to encode the bits is based off the total counts of bits.
    /// TODO(fgalligan): Investigate using an adaptive table for more compression.
    /// </summary>
    class RAnsBitEncoder : RAnsBitCodec, IBitEncoder
    {
        private ulong[] bitCounts = new ulong[2];
        private List<uint> bits = new List<uint>();
        private uint localBits;
        private uint numLocalBits;


        struct AnsCoder
        {
            internal byte[] buf;
            internal int bufOffset;
            internal uint state;
        }

        /// <summary>
        /// Must be called before any Encode* function is called.
        /// </summary>
        public void StartEncoding()
        {
            Clear();
        }

        /// <summary>
        /// Encode one bit. If |bit| is true encode a 1, otherwise encode a 0.
        /// </summary>
        public void EncodeBit(bool bit)
        {

            if (bit)
            {
                bitCounts[1]++;
                localBits |= (uint)(1 << (int)numLocalBits);
            }
            else
            {
                bitCounts[0]++;
            }
            numLocalBits++;

            if (numLocalBits == 32)
            {
                bits.Add(localBits);
                numLocalBits = 0;
                localBits = 0;
            }
        }

        /// <summary>
        /// Encode |nibts| of |value|, starting from the least significant bit.
        /// |nbits| must be &gt; 0 and &lt;= 32.
        /// </summary>
        public void EncodeLeastSignificantBits32(int nbits, uint value)
        {
            uint reversed = DracoUtils.ReverseBits32(value) >> (32 - nbits);
            int ones = DracoUtils.CountOnes32(reversed);
            bitCounts[0] += (ulong)(nbits - ones);
            bitCounts[1] += (ulong)ones;

            int remaining = (int)(32 - numLocalBits);

            if (nbits <= remaining)
            {
                DracoUtils.CopyBits32(ref localBits, (int)numLocalBits, reversed, 0, nbits);
                numLocalBits += (uint)(nbits);
                if (numLocalBits == 32)
                {
                    bits.Add(localBits);
                    localBits = 0;
                    numLocalBits = 0;
                }
            }
            else
            {
                DracoUtils.CopyBits32(ref localBits, (int)numLocalBits, reversed, 0, remaining);
                bits.Add(localBits);
                localBits = 0;
                DracoUtils.CopyBits32(ref localBits, 0, reversed, remaining, nbits - remaining);
                numLocalBits = (uint)(nbits - remaining);
            }
        }

        /// <summary>
        /// Ends the bit encoding and stores the result into the targetBuffer.
        /// </summary>
        public void EndEncoding(EncoderBuffer targetBuffer)
        {

            ulong total = bitCounts[1] + bitCounts[0];
            if (total == 0)
                total++;

            // The probability interval [0,1] is mapped to values of [0, 256]. However,
            // the coding scheme can not deal with probabilities of 0 or 1, which is why
            // we must clamp the values to interval [1, 255]. Specifically 128
            // corresponds to 0.5 exactly. And the value can be given as uint8T.
            uint zeroProbRaw = (uint) (
                ((bitCounts[0] / (double) total) * 256.0) + 0.5);

            byte zeroProb = 255;
            if (zeroProbRaw < 255)
                zeroProb = (byte) (zeroProbRaw);

            zeroProb += (byte) ((zeroProb == 0) ? 1 : 0);

            // Space for 32 bit integer and some extra space.
            byte[] buffer = new byte[(bits.Count + 8) * 8];
            AnsCoder ansCoder = new AnsCoder();
            ansWriteInit(ref ansCoder, buffer);

            for (int i = (int)numLocalBits - 1; i >= 0; --i)
            {
                uint bit = (localBits >> i) & 1;
                rabsDescWrite(ref ansCoder, bit, zeroProb);
            }
            for(int j = bits.Count - 1; j >= 0; j--)
            {
                uint nbits = this.bits[j];
                for (int i = 31; i >= 0; --i)
                {
                    uint bit = (nbits >> i) & 1;
                    rabsDescWrite(ref ansCoder, bit, zeroProb);
                }
            }

            int sizeInBytes = ansWriteEnd(ref ansCoder);
            targetBuffer.Encode(zeroProb);
            Encoding.EncodeVarint(sizeInBytes, targetBuffer);
            targetBuffer.Encode(buffer, sizeInBytes);

            Clear();
        }

        /// <summary>
        /// rABS with descending spread
        /// p or p0 takes the place of lS from the paper
        /// ansP8Precision is m
        /// </summary>
        static void rabsDescWrite(ref AnsCoder ans, uint val, byte p0)
        {
            var p = ansP8Precision - p0;
            var lS = val == 1 ? p : p0;
            uint quot, rem;
            if (ans.state >= lBase / ansP8Precision * ioBase * lS)
            {
                ans.buf[ans.bufOffset++] = (byte)(ans.state % ioBase);
                ans.state /= ioBase;
            }

            quot = (uint)(ans.state / lS);
            rem = (uint)(ans.state - quot * lS);

            ans.state = (uint)(quot * ansP8Precision + rem + (val == 1 ? 0 : p));
        }

        static void ansWriteInit(ref AnsCoder ans, byte[] buf)
        {
            ans.buf = buf;
            ans.bufOffset = 0;
            ans.state = lBase;
        }

        static int ansWriteEnd(ref AnsCoder ans)
        {
            uint state;
            //assert(ans->state >= lBase);
            //assert(ans->state < lBase * ioBase);
            state = ans.state - lBase;
            if (state < (1 << 6))
            {
                ans.buf[ans.bufOffset] = (byte)((0x00 << 6) + state);
                return ans.bufOffset + 1;
            }
            else if (state < (1 << 14))
            {
                Unsafe.PutLE16(ans.buf, ans.bufOffset, (ushort)((0x01 << 14) + state));
                return ans.bufOffset + 2;
            }
            else if (state < (1 << 22))
            {
                Unsafe.PutLE24(ans.buf, ans.bufOffset, (uint)((0x02 << 22) + state));
                return ans.bufOffset + 3;
            }
            else
            {
                throw new Exception("Invalid RAns state");
            }
        }

        public void Clear()
        {

            bitCounts[0] = 0;
            bitCounts[1] = 0;
            bits.Clear();
            localBits = 0;
            numLocalBits = 0;
        }
    }
}
