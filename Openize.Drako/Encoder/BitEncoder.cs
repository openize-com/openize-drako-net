using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Encoder
{
    class BitEncoder
    {
        private BytePointer buffer;
        private int bitOffset;

        public BitEncoder(BytePointer buffer)
        {
            this.buffer = buffer;
        }


        /// <summary>
        /// Write |nbits| of |data| into the bit buffer.
        /// </summary>
        public void PutBits(uint data, int nbits)
        {
            for (int bit = 0; bit < nbits; ++bit)
                PutBit((byte) ((data >> bit) & 1));
        }

        /// <summary>
        /// Return number of bits encoded so far.
        /// </summary>
        public int Bits
        {
            get { return bitOffset; }
        }
        private void PutBit(byte value)
        {
            int byteSize = 8;
            int off = bitOffset;
            int byteOffset = off / byteSize;
            int bitShift = off % byteSize;

            // TODO(fgalligan): Check performance if we add a branch and only do one
            // memory write if bitShift is 7. Also try using a temporary variable to
            // hold the bits before writing to the buffer.
            buffer[byteOffset] = (byte)((buffer[byteOffset] & ~(1 << bitShift)) | (value << bitShift));

            byte t = buffer[byteOffset];
            t = (byte) (t & ~(1 << bitShift));
            t = (byte) (t | value << bitShift);
            buffer[byteOffset] = t;
            bitOffset++;
        }
    }
}
