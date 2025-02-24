using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Encoder
{

    /// <summary>
    /// Class that encodes mesh data using a simple binary representation of mesh's
    /// connectivity and geometry.
    /// </summary>
    class MeshSequentialEncoder : MeshEncoder
    {
        public override DracoEncodingMethod GetEncodingMethod()
        {
            return DracoEncodingMethod.Sequential;
        }

        protected override void EncodeConnectivity()
        {
            // Serialize indices.
            int numFaces = Mesh.NumFaces;
            Encoding.EncodeVarint(numFaces, Buffer);
            Encoding.EncodeVarint(Mesh.NumPoints, Buffer);

            // We encode all attributes in the original (possibly duplicated) format.
            // TODO(ostava): This may not be optimal if we have only one attribute or if
            // all attributes share the same index mapping.
#if COMPRESS_CONNECTIVITY
            const bool compressConnectivity = false;
            if (compressConnectivity)
            {
                // 0 = Encode compressed indices.
                Buffer.Encode((byte) (0));
                if (!CompressAndEncodeIndices())
                    return false;
            }
            else
#endif
            {
                // 1 = Encode indices directly.
                Buffer.Encode((byte) 1);
                // Store vertex indices using a smallest datatype that fits their range.
                // TODO(ostava): This can be potentially improved by using a tighter
                // fit that is not bound by a bit-length of any particular data type.
                Span<int> face = stackalloc int[3];
                if (Mesh.NumPoints < 256)
                {
                    // Serialize indices as uint8T.
                    for (int i = 0; i < numFaces; ++i)
                    {
                        Mesh.ReadFace(i, face);
                        Buffer.Encode((byte) (face[0]));
                        Buffer.Encode((byte) (face[1]));
                        Buffer.Encode((byte) (face[2]));
                    }
                }
                else if (Mesh.NumPoints < (1 << 16))
                {
                    // Serialize indices as uint16T.
                    for (int i = 0; i < numFaces; ++i)
                    {
                        Mesh.ReadFace(i, face);
                        Buffer.Encode((ushort) (face[0]));
                        Buffer.Encode((ushort) (face[1]));
                        Buffer.Encode((ushort) (face[2]));
                    }
                }
                else if (Mesh.NumPoints < (1 << 21))
                {
                    // Serialize indices as varint.
                    for (int i = 0; i < numFaces; ++i)
                    {
                        Mesh.ReadFace(i, face);
                        Encoding.EncodeVarint(face[0], buffer);
                        Encoding.EncodeVarint(face[1], buffer);
                        Encoding.EncodeVarint(face[2], buffer);
                    }
                }
                else
                {
                    // Serialize faces as uint (default).
                    for (int i = 0; i < numFaces; ++i)
                    {
                        Mesh.ReadFace(i, face);
                        Buffer.Encode(face);
                    }
                }
            }
        }

        protected override void GenerateAttributesEncoder(int attId)
        {
            // Create only one attribute encoder that is going to encode all points in a
            // linear sequence.
            if (attId == 0)
            {
                // Create a new attribute encoder only for the first attribute.
                AddAttributesEncoder(new SequentialAttributeEncodersController(
                    new LinearSequencer(PointCloud.NumPoints), attId));
            }
            else
            {
                // Reuse the existing attribute encoder for other attributes.
                AttributesEncoder(0).AddAttributeId(attId);
            }
        }

        private bool CompressAndEncodeIndices()
        {
            // Collect all indices to a buffer and encode them.
            // Each new indice is a difference from the previous value.
            int numFaces = Mesh.NumFaces;
            Span<int> indicesBuffer = new int[3 * numFaces];
            int lastIndexValue = 0;
            int p = 0;
            Span<int> face = stackalloc int[3];
            for (int i = 0; i < numFaces; ++i)
            {
                Mesh.ReadFace(i, face);
                for (int j = 0; j < 3; ++j)
                {
                    int indexValue = face[j];
                    int indexDiff = indexValue - lastIndexValue;
                    // Encode signed value to an unsigned one (put the sign to lsb pos).
                    int encodedVal = (Math.Abs(indexDiff) << 1) | (indexDiff < 0 ? 1 : 0);
                    indicesBuffer[p++] = encodedVal;
                    lastIndexValue = indexValue;
                }
            }
            Encoding.EncodeSymbols(indicesBuffer, indicesBuffer.Length, 1, null, Buffer);
            return true;
        }
    }
}
