using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Utils
{
    /// <summary>
    /// This simulates a byte pointer used in Draco implementation, also makes it easier to be ported to Java using CsPorter made by Lex Chou. 
    /// I've benchmarked this, it's okay to be used, I'll replace this by Span later.
    /// .NET Version: 8.0
    /// Build Type: Release
    /// Test memory block: 1GB
    /// Random read/write counts: 10M
    /// Test PC Setup: AMD 3700X, 96G Mem
    ///
    /// All benchmark tests are pre-warmed in JIT and memory access.
    /// The benchmark result:
    ///
    /// Array Write         00:00:00.2786953
    /// Unsafe Write        00:00:00.2074192
    /// Span Write          00:00:00.1917591
    /// BytePointer Write   00:00:00.2552571
    /// Memory Slice Write  00:00:00.2477228 
    /// Memory Span Write   00:00:00.2477730
    ///
    /// Array Read          00:00:00.3126299 
    /// Unsafe Read         00:00:00.3272010 
    /// Span Read           00:00:00.3320909
    /// BytePointer Read    00:00:00.3381226
    /// Memory Span Read    00:00:00.3725757
    /// Memory Slice Read   00:00:00.6809910
    /// </summary>
    struct BytePointer
    {
        private byte[] data;
        private int offset;

        public BytePointer(byte[] data)
        {
            this.data = data;
            this.offset = 0;
        }
        public BytePointer(byte[] data, int offset)
        {
            this.data = data;
            this.offset = offset;
        }

        public int Offset
        {
            get { return offset; }
        }

        public byte[] BaseData
        {
            get { return data; }
        }

        public byte this[int offset]
        {
            get { return data[this.offset + offset];}
            set { data[this.offset + offset] = value;}
        }

        public byte ToByte()
        {
            return data[this.offset];
        }
        public ushort ToUInt16LE()
        {
            return Unsafe.GetLE16(data, this.offset);
        }
        public ushort ToUInt16LE(int offset)
        {
            return Unsafe.GetLE16(data, this.offset + offset);
        }
        public uint ToUInt24LE(int offset)
        {
            return Unsafe.GetLE24(data, this.offset + offset);
        }
        public uint ToUInt32LE(int offset)
        {
            return Unsafe.GetLE32(data, this.offset + offset);
        }
        public ulong ToUInt64LE(int offset)
        {
            return Unsafe.GetLE64(data, this.offset + offset);
        }
        public float ToSingle(int offset)
        {
            return Unsafe.GetFloat(data, this.offset + offset);
        }

        public bool IsOverflow(int offset)
        {
            int p = offset + this.offset;
            return (p >= data.Length) || p < 0;
        }

        public static BytePointer operator +(BytePointer ptr, int offset)
        {
            return new BytePointer(ptr.data, ptr.offset + offset);
        }

        public void Copy(int srcOffset, byte[] dst, int dstOffset, int len)
        {
            Array.Copy(data, this.offset + srcOffset, dst, dstOffset, len);
        }

        public override string ToString()
        {
            return string.Format("byte[{0}]+{1}", data == null ? 0 : data.Length, offset);
        }
    }
}
