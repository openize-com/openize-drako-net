using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Decoder
{
    sealed class DecoderBuffer
    {
        private readonly BitDecoder bitDecoder = new BitDecoder();
        private bool bitMode = false;
        private int pos;
        private int length;
        private BytePointer data;
        private byte[] tmp = new byte[8];
        public int BitstreamVersion { get; set; }

        public DecoderBuffer(byte[] data)
        {
            Initialize(new BytePointer(data), data.Length);
        }

        private DecoderBuffer(BytePointer data, int length)
        {
            Initialize(data, length);
        }

        public DecoderBuffer()
        {

        }

        public bool BitDecoderActive => bitMode;

        public void CopyFrom(DecoderBuffer src)
        {
            bitDecoder.CopyFrom(src.bitDecoder);
            bitMode = src.bitMode;
            pos = src.pos;
            length = src.length;
            data = src.data;
            BitstreamVersion = src.BitstreamVersion;
        }
        public DecoderBuffer Clone()
        {
            var ret = new DecoderBuffer();
            ret.CopyFrom(this);
            return ret;
        }

        private void Initialize(BytePointer data, int length)
        {
            pos = 0;
            bitMode = false;
            
            this.data = data;
            this.length = length;
        }

        public int RemainingSize
        {
            get { return length -pos; }
        }

        public DecoderBuffer SubBuffer(int offset)
        {
            int length = this.length - this.pos - offset;
            var ret = new DecoderBuffer(data + (pos + offset), length);
            ret.BitstreamVersion = BitstreamVersion;
            return ret;
        }


        public int DecodedSize
        {
            get { return pos; }
        }

        public int BufferSize
        {
            get { return length; }
        }

        public BytePointer Pointer
        {
            get { return data;}
        }
        /// <summary>
        /// Discards #bytes from the input buffer.
        /// </summary>
        public void Advance(int bytes)
        {
            pos += bytes;
        }

        public long StartBitDecoding(bool decodeSize)
        {
            long outSize = 0;
            if (decodeSize)
            {
                if (BitstreamVersion < 22)
                {
                    outSize = DecodeI64();
                }
                else
                {
                    ulong n = Decoding.DecodeVarintU64(this);
                    outSize = (long)n;
                }
            }
            bitMode = true;
            bitDecoder.Load(data + pos, length - pos);
            return outSize;
        }

        public void EndBitDecoding()
        {
            bitMode = false;
            int bitsDecoded = bitDecoder.BitsDecoded;
            //some bits in last byte are discard
            int bytesDecoded = (bitsDecoded + 7)/8;
            pos += bytesDecoded;
        }

        public bool Decode(byte[] buf, int len)
        {
            return Decode(buf, 0, len);
        }

        public bool Decode(uint[] values)
        {
            int bytes = values.Length * 4;

            if (!RemainingIsEnough(bytes))
                return false;
            int n = data.Offset + pos;
            var buf = data.BaseData;
            for (int i = 0; i < values.Length; i++)
            {
                values[i] = Unsafe.GetLE32(buf, n);
                n += 4;
                pos += 4;
            }
            return true;
        }
        public bool Decode(float[] floats)
        {
            int bytes = floats.Length * 4;
            if(bytes >= tmp.Length)
                tmp = new byte[bytes];
            if (!Decode(tmp, bytes))
                return false;
            Unsafe.ToFloatArray(tmp, floats);
            return true;
        }
        public bool Decode(byte[] buf)
        {
            return Decode(buf, 0, buf.Length);
        }

        public bool Decode(byte[] buf, int start, int len)
        {
            if (!RemainingIsEnough(len))
                return false;
            this.data.Copy(pos, buf, start, len);
            pos += len;
            return true;
        }

        public float DecodeF32()
        {
            if (!Decode(tmp, 4))
                throw DracoUtils.Failed();
            return BitConverter.ToSingle(tmp, 0);
        }
        public sbyte DecodeI8()
        {
            return (sbyte)DecodeU8();
        }
        public byte DecodeU8()
        {
            if (!RemainingIsEnough(1))
                throw DracoUtils.Failed();
            var val = this.data[pos];
            pos++;
            return val;
        }
        public ushort DecodeU16()
        {
            if (!RemainingIsEnough(2))
                throw DracoUtils.Failed();
            var val = data.ToUInt16LE(pos);
            pos += 2;
            return val;
        }
        public long DecodeI64()
        {
            if (!RemainingIsEnough(8))
                throw DracoUtils.Failed();
            var val = (long)data.ToUInt64LE(pos);
            pos += 8;
            return val;
        }

        public uint DecodeU32()
        {
            return (uint)DecodeI32();
        }
        public int DecodeI32()
        {
            if (!RemainingIsEnough(4))
                throw DracoUtils.Failed();
            var val = (int)data.ToUInt32LE(pos);
            pos += 4;
            return val;
        }

        private bool RemainingIsEnough(int size)
        {
            return length >= pos + size;
        }
        public bool Peek(byte[] result, int size)
        {
            if (data.IsOverflow(pos + size))
                return false;//overflow
            this.data.Copy(pos, result, 0, size);
            return true;
        }

        /// <summary>
        /// Decodes up to 32 bits into outVal. Can be called only in between
        /// StartBitDecoding and EndBitDeoding. Otherwise returns false.
        /// </summary>
        /// <param name="nbits"></param>
        /// <param name="value"></param>
        /// <returns></returns>
        public bool DecodeLeastSignificantBits32(int nbits, out uint value)
        {
            if (!bitMode)
            {
                value = 0;
                return false;
            }
            value = bitDecoder.GetBits(nbits);
            return true;
        }

        public override string ToString()
        {
            int preview = Math.Min(16, RemainingSize);
            StringBuilder sb = new StringBuilder();
            sb.AppendFormat("0x{0:x4} ", pos + data.Offset);
            for (int i = 0, p = pos; i < preview; i++, p++)
            {
                sb.AppendFormat("{0:x2} ", data[p]);
            }
            int rest = RemainingSize - preview;
            if (rest > 0)
            {
                sb.AppendFormat("...({0} bytes rest)", rest);
            }
            sb.AppendFormat(", pos = {0}", pos);
            return sb.ToString();
        }
    }
}
