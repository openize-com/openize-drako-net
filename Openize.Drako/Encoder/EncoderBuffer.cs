using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Encoder
{

    /// <summary>
    /// Class representing a buffer that can be used for either for byte-aligned
    /// encoding of arbitrary data structures or for encoding of varialble-length
    /// bit data.
    /// </summary>
    class EncoderBuffer
    {
        private BitEncoder bitEncoder;
        private DataBuffer buffer = new DataBuffer();

        /// <summary>
        /// The number of bytes reserved for bit encoder.
        /// Values > 0 indicate we are in the bit encoding mode.
        /// </summary>
        private long bitEncoderReservedBytes;

        /// <summary>
        /// Flag used indicating that we need to store the length of the currently
        /// processed bit sequence.
        /// </summary>
        private bool encodeBitSequenceSize;

        public void Encode(ushort val)
        {
            int offset = buffer.Length;
            DebugBreak(2);
            buffer.Length += 2;
            Unsafe.PutLE16(buffer.GetBuffer(), offset, val);
        }
        public bool Encode(byte val)
        {
            int offset = buffer.Length;
            DebugBreak(1);
            buffer.Length += 1;
            buffer[offset] = val;
            return true;
        }
        public void Encode(uint val)
        {
            int offset = buffer.Length;
            DebugBreak(4);
            buffer.Length += 4;
            Unsafe.PutLE32(buffer.GetBuffer(), offset, val);
        }

        public void Encode(float val)
        {
            int offset = buffer.Length;
            DebugBreak(4);
            buffer.Length += 4;
            Unsafe.PutLE32(buffer.GetBuffer(), offset, Unsafe.FloatToUInt32(val));
        }
        public void Encode(float[] val)
        {
            int offset = buffer.Length;
            DebugBreak(4 * val.Length);
            buffer.Length += 4 * val.Length;
            for (int i = 0; i < val.Length; i++)
            {
                Unsafe.PutLE32(buffer.GetBuffer(), offset, Unsafe.FloatToUInt32(val[i]));
                offset += 4;
            }
        }
        public void Encode(Span<int> val)
        {
            int offset = buffer.Length;
            DebugBreak(4 * val.Length);
            buffer.Length += 4 * val.Length;
            for (int i = 0; i < val.Length; i++)
            {
                Unsafe.PutLE32(buffer.GetBuffer(), offset, (uint)val[i]);
                offset += 4;
            }
        }

        public void Encode(IntList val)
        {
            Encode(val.data, val.Count);
        }

        public void Encode(int[] val)
        {
            Encode(val, val.Length);
        }
        public void Encode(int[] val, int len)
        {
            int offset = buffer.Length;
            DebugBreak(4 * len);
            buffer.Length += 4 * len;
            for (int i = 0; i < len; i++)
            {
                Unsafe.PutLE32(buffer.GetBuffer(), offset, (uint)val[i]);
                offset += 4;
            }
        }

        public void Encode(int val)
        {
            Encode((uint)val);
        }

        public void Encode(byte[] buffer, int length)
        {
            int offset = this.buffer.Length;
            DebugBreak(length);
            this.buffer.Length += length;
            Array.Copy(buffer, 0, this.buffer.GetBuffer(), offset, length);
        }
        public void Encode(Span<byte> buffer, int start, int length)
        {
            int offset = this.buffer.Length;
            DebugBreak(length);
            this.buffer.Length += length;
            buffer.Slice(start, length).CopyTo(this.buffer.GetBuffer().AsSpan().Slice(offset));
        }

        [Conditional("DEBUG")]
        private void DebugBreak(int len)
        {
            /*
            int debugOffset = 29;
            int offset = this.buffer.Length;
            if (debugOffset >= offset && debugOffset < offset + len)
                Debugger.Break();
            */
        }

        public void Clear()
        {

            buffer.Clear();
            bitEncoderReservedBytes = 0;
        }

        public void Resize(int nbytes)
        {
            DebugBreak(nbytes - buffer.Length);
            buffer.Length = nbytes;
        }

        /// <summary>
        /// Start encoding a bit sequence. A maximum size of the sequence needs to
        /// be known upfront.
        /// If encodeSize is true, the size of encoded bit sequence is stored before
        /// the sequence. Decoder can then use this size to skip over the bit sequence
        /// if needed.
        /// Returns false on error.
        /// </summary>
        public bool StartBitEncoding(int requiredBits, bool encodeSize)
        {

            if (BitEncoderActive)
                return false; // Bit encoding mode already active.
            if (requiredBits <= 0)
                return false; // Invalid size.
            encodeBitSequenceSize = encodeSize;
            int requiredBytes = (requiredBits + 7) / 8;
            bitEncoderReservedBytes = requiredBytes;
            int bufferStartSize = buffer.Length;
            if (encodeSize)
            {
                // Reserve memory for storing the encoded bit sequence size. It will be
                // filled once the bit encoding ends.
                bufferStartSize += 8; //sizeof(long)
            }
            // Resize buffer to fit the maximum size of encoded bit data.
            DebugBreak(requiredBytes);
            buffer.Length = bufferStartSize + requiredBytes;
            // Get the buffer data pointer for the bit encoder.
            BytePointer data = new BytePointer(Data, (int)bufferStartSize);
            bitEncoder = new BitEncoder(data);
            return true;
        }

        /// <summary>
        /// End the encoding of the bit sequence and return to the default byte-aligned
        /// encoding.
        /// </summary>
        public void EndBitEncoding()
        {

            if (!BitEncoderActive)
                return;
            // Get the number of encoded bits and bytes (rounded up).
            long encodedBits = bitEncoder.Bits;
            long encodedBytes = (encodedBits + 7) / 8;
            // Encode size if needed.
            if (encodeBitSequenceSize)
            {
                int out_mem = (int)(this.Bytes - (bitEncoderReservedBytes + 8));
                EncoderBuffer var_size_buffer = new EncoderBuffer();
                Encoding.EncodeVarint((ulong)encodedBytes, var_size_buffer);

                int size_len = var_size_buffer.Bytes;
                int dst = out_mem + size_len;
                int src = out_mem + 8 /*sizeof(uint64_t)*/;
                Array.Copy(Data, src, Data, dst, encodedBytes);
                // Store the size of the encoded data.
                Array.Copy(var_size_buffer.buffer.GetBuffer(), 0, Data, out_mem, size_len);

                // We need to account for the difference between the preallocated and actual
                // storage needed for storing the encoded length. This will be used later to
                // compute the correct size of |buffer_|.
                bitEncoderReservedBytes += 8 /*sizeof uint64_t*/ - size_len;
            }

            // Resize the underlying buffer to match the number of encoded bits.
            buffer.Length = (int) (buffer.Length - bitEncoderReservedBytes + encodedBytes);
            bitEncoderReservedBytes = 0;
        }

        /// <summary>
        /// Encode up to 32 bits into the buffer. Can be called only in between
        /// StartBitEncoding and EndBitEncoding. Otherwise returns false.
        /// </summary>
        public bool EncodeLeastSignificantBits32(int nbits, uint value)
        {
            if (!BitEncoderActive)
                return false;
            bitEncoder.PutBits(value, nbits);
            return true;
        }

        public bool BitEncoderActive { get { return bitEncoderReservedBytes > 0; }}

        internal void Encode(Span<int> ints, int bytesOffset, int bytes)
        {
            int offset = this.buffer.Length;
            this.buffer.Length += bytes;
            //Buffer.BlockCopy(ints.data, ints.ByteOffset * 4, this.buffer.GetBuffer(), offset, bytes);
            var dst = MemoryMarshal.Cast<byte, int>(buffer.GetBuffer().AsSpan(offset, bytes));
            ints.Slice(bytesOffset / 4, bytes / 4).CopyTo(dst);

        }
        internal void Encode(Span<int> ints, int bytes)
        {
            Encode(ints, 0, bytes);
        }

        public BitEncoder BitEncoder
        {
            get { return bitEncoder; }
        }

        public int Bytes { get { return buffer.Length; } }
        public byte[] Data { get { return buffer.GetBuffer(); } }
    }
}
