using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Encoder
{

    /// <summary>
    /// Class for performing rANS encoding using a desired number of precision bits.
    /// The max number of precision bits is currently 19. The actual number of
    /// symbols in the input alphabet should be (much) smaller than that, otherwise
    /// the compression rate may suffer.
    /// </summary>
    class RAnsEncoder : RAnsBitCodec
    {
        private int precisionBits;
        private int ransPrecision;
        private uint lRansBase;
        private BytePointer buf;
        private int bufOffset;
        private uint state;

        public RAnsEncoder(int precisionBits)
        {
            this.precisionBits = precisionBits;
            ransPrecision = 1 << precisionBits;
            lRansBase = (uint)ransPrecision * 4;
        }


        public void Reset(BytePointer ptr)
        {
            this.buf = ptr;
            this.bufOffset = 0;
            this.state = lRansBase ;
        }

        // Needs to be called after all symbols are encoded.
        public int writeEnd()
        {
            uint state;
            //assert(ans.state >= lRansBase);
            //assert(ans.state < lRansBase * ioBase);
            state = this.state - lRansBase;
            if (state < (1 << 6))
            {
                this.buf[this.bufOffset] = (byte)((0x00 << 6) + state);
                return this.bufOffset + 1;
            }
            else if (state < (1 << 14))
            {
                Unsafe.PutLE16(this.buf.BaseData, buf.Offset + bufOffset, (ushort)((0x01 << 14) + state));
                return this.bufOffset + 2;
            }
            else if (state < (1 << 22))
            {
                Unsafe.PutLE24(this.buf.BaseData, buf.Offset + bufOffset, (0x02 << 22) + state);
                return this.bufOffset + 3;
            }
            else if (state < (1 << 30))
            {
                Unsafe.PutLE32(this.buf.BaseData, buf.Offset + bufOffset, (uint) ((0x03 << 30) + state));
                return bufOffset + 4;
            }
            else
            {
                throw new Exception("Invalid rANS state.");
            }
        }

        /// <summary>
        /// rANS with normalization
        /// sym->prob takes the place of lS from the paper
        /// ransPrecision is m
        /// </summary>
        public void Write(RansSym sym)
        {
            uint p = sym.prob;
            while (this.state >= lRansBase / ransPrecision * ioBase * p)
            {
                buf[bufOffset++] = (byte) (this.state % ioBase);
                state /= ioBase;
            }
            // TODO(ostava): The division and multiplication should be optimized.
            state = (uint)((state / p) * ransPrecision + state % p + sym.cumProb);
        }
    }
}
