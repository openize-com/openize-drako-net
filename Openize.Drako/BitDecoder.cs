using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako
{
    sealed class BitDecoder
    {
        private BytePointer data;
        private int dataEnd;
        private int bitOffset;

        public void CopyFrom(BitDecoder bitDecoder)
        {
            this.data = bitDecoder.data;
            this.dataEnd = bitDecoder.dataEnd;
            this.bitOffset = bitDecoder.bitOffset;
        }
        public void Load(BytePointer data, int count)
        {
            this.data = data;
            dataEnd = count;
            bitOffset = 0;
        }

        public int BitsDecoded
        {
            get { return bitOffset; }
        }

        public void Consume(int k)
        {
            bitOffset += k;
        }


        public int GetBit()
        {
            int off = bitOffset;
            int byteOffset = (off >> 3);
            int bitShift = (int) (off & 0x7);
            if (byteOffset < dataEnd)
            {
                int bit = (data[byteOffset] >> bitShift) & 1;
                bitOffset++;
                return bit;
            }
            return 0;
        }

        public int PeekBit(int offset)
        {
            int off = bitOffset + offset;
            int byteOffset = (off >> 3);
            int bitShift = (int) (off & 0x7);
            if (byteOffset < dataEnd)
            {
                int bit = (data[byteOffset] >> bitShift) & 1;
                return bit;
            }
            return 0;
        }

        public uint GetBits(int nbits)
        {
            uint ret = 0;
            for (int bit = 0; bit < nbits; ++bit)
                ret |= (uint) (GetBit() << bit);
            return ret;
        }

    }
}
