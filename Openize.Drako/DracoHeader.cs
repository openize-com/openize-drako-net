using FileFormat.Drako.Decoder;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako
{
    class DracoHeader
    {
        public const int MetadataFlagMask = 0x8000;
        private static readonly byte[] HEADER = { (byte)'D', (byte)'R', (byte)'A', (byte)'C', (byte)'O' };
        private byte major;
        private byte minor;
        public int version;
        public ushort flags;

        public EncodedGeometryType encoderType;
        public DracoEncodingMethod method;

        public static DracoHeader Parse(DecoderBuffer buffer)
        {
            byte[] header = new byte[5];
            if (!buffer.Decode(header))
                return null;
            if (DracoUtils.Compare(header, HEADER, HEADER.Length) != 0)
                return null;
            DracoHeader ret = new DracoHeader();
            ret.encoderType = EncodedGeometryType.Invalid;
            ret.method = DracoEncodingMethod.Sequential;
            ret.major = buffer.DecodeU8();
            ret.minor = buffer.DecodeU8();
            ret.version = ret.major * 10 + ret.minor;
            byte t = buffer.DecodeU8();
            ret.encoderType = (EncodedGeometryType)t;
            t = buffer.DecodeU8();
            ret.method = (DracoEncodingMethod)t;
            ret.flags = buffer.DecodeU16();
            return ret;
        }
    }
}
