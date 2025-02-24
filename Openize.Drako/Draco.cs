using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using FileFormat.Drako;
using FileFormat.Drako.Decoder;
using FileFormat.Drako.Encoder;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako
{
    /// <summary>
    /// Google Draco
    /// </summary>
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    class Draco
    {

        /// <summary>
        /// Decode a <see cref="DracoPointCloud"/> or <see cref="DracoMesh"/> from bytes
        /// </summary>
        /// <param name="data">Raw draco bytes.</param>
        /// <returns>a <see cref="DracoPointCloud"/> or <see cref="DracoMesh"/> instance</returns>
        public static DracoPointCloud Decode(byte[] data)
        {
            if (data == null)
                throw new ArgumentNullException("data");
            DecoderBuffer buffer = new DecoderBuffer(data);
            return DracoMesh.Decode(buffer);
        }

        /// <summary>
        /// Encode the point cloud or mesh and get the encoded bytes in draco format.
        /// </summary>
        /// <param name="m">The <see cref="DracoPointCloud"/> or <see cref="DracoMesh"/> instance </param>
        /// <returns>Bytes in draco format</returns>
        public static byte[] Encode(DracoPointCloud m)
        {
            return Encode(m, new DracoEncodeOptions());
        }
        /// <summary>
        /// Encode the point cloud or mesh and get the encoded bytes in draco format.
        /// </summary>
        /// <param name="m">The <see cref="DracoPointCloud"/> or <see cref="DracoMesh"/> instance </param>
        /// <param name="options">Encode options</param>
        /// <returns>Bytes in draco format</returns>
        public static byte[] Encode(DracoPointCloud m, DracoEncodeOptions options)
        {
            if (m == null)
                throw new ArgumentNullException("m");
            if (options == null)
                throw new ArgumentNullException("options");
            var buf = EncodeImpl(m, options);
            if (buf.Data.Length == buf.Bytes)
                return buf.Data;
            else
            {
                var ret = new byte[buf.Bytes];
                Array.Copy(buf.Data, ret, buf.Bytes);
                return ret;
            }
        }
        public static void Encode(DracoPointCloud m, DracoEncodeOptions options, Stream stream)
        {
            if (m == null)
                throw new ArgumentNullException("m");
            if (options == null)
                throw new ArgumentNullException("options");
            if (stream == null)
                throw new ArgumentNullException("stream");

            var buf = EncodeImpl(m, options);
            stream.Write(buf.Data, 0, buf.Bytes);
        }
        internal static EncoderBuffer EncodeImpl(DracoPointCloud m, DracoEncodeOptions options)
        {
            EncoderBuffer ret = new EncoderBuffer();

            var encoder = CreateEncoder(m, options);
            //Encode header

            // Encode the header according to our v1 specification.
            // Five bytes for Draco format.
            ret.Encode(new byte[] {(byte) 'D', (byte) 'R', (byte) 'A', (byte) 'C', (byte) 'O'}, 5);
            // Version (major, minor).
            byte majorVersion;
            byte minorVersion;
            if(m is DracoMesh)
            {
                majorVersion = 2;
                minorVersion = 2;
            }
            else
            {
                //point cloud
                majorVersion = 2;
                minorVersion = 3;
            }
            ret.Encode(majorVersion);
            ret.Encode(minorVersion);
            // Type of the encoder (point cloud, mesh, ...).
            ret.Encode((byte) encoder.GeometryType);
            // Unique identifier for the selected encoding method (edgebreaker, etc...).
            ret.Encode((byte) encoder.GetEncodingMethod());
            // Reserved for flags.
            ret.Encode((ushort) 0);

            //encode body
            encoder.Encode(options, ret);

            return ret;
        }

        private static PointCloudEncoder CreateEncoder(DracoPointCloud pc, DracoEncodeOptions options)
        {
            if (pc is DracoMesh && ((DracoMesh)pc).NumFaces > 0)
            {
                MeshEncoder encoder;
                if (options.CompressionLevel == DracoCompressionLevel.NoCompression)
                    encoder = new MeshSequentialEncoder();
                else
                    encoder = new MeshEdgeBreakerEncoder();
                encoder.Mesh = (DracoMesh) pc;
                return encoder;
            }
            else
            {
                //check if kd-tree is possible
                // Speed < 10, use POINT_CLOUD_KD_TREE_ENCODING if possible.
                if (!IsKdTreePossible(pc, options))
                    throw new InvalidOperationException("KD Tree encoder is not supported on this point cloud.");
                var ret = new PointCloudKdTreeEncoder();
                ret.PointCloud = pc;
                return ret;
            }
        }

        private static bool IsKdTreePossible(DracoPointCloud pc, DracoEncodeOptions options)
        {

            // Kd-Tree encoder can be currently used only when the following conditions
            // are satisfied for all attributes:
            //     -data type is float32 and quantization is enabled, OR
            //     -data type is uint32, uint16, uint8 or int32, int16, int8
            for (int i = 0; i < pc.NumAttributes; ++i)
            {
                PointAttribute att = pc.Attribute(i);
                if (att.DataType != DataType.FLOAT32 &&
                    att.DataType != DataType.UINT32 && att.DataType != DataType.UINT16 &&
                    att.DataType != DataType.UINT8 && att.DataType != DataType.INT32 &&
                    att.DataType != DataType.INT16 && att.DataType != DataType.INT8)
                    return false;
                if (att.DataType == DataType.FLOAT32 &&
                    options.GetQuantizationBits(att) <= 0)
                    return false; // Quantization not enabled.
            }

            return true;
        }
    }
}
