using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Decoder
{
    class MeshSequentialDecoder : MeshDecoder
    {
        protected override void DecodeConnectivity()
        {
            uint numFaces;
            uint numPoints;
            if (BitstreamVersion < 22)
            {
                numFaces = Buffer.DecodeU32();
                numPoints = Buffer.DecodeU32();
            }
            else
            {
                numFaces = Decoding.DecodeVarintU32(buffer);
                numPoints = buffer.DecodeVarintU32();
            }

            byte connectivityMethod = Buffer.DecodeU8();
            if (connectivityMethod == 0)
            {
                DecodeAndDecompressIndices((int)numFaces);
            }
            else
            {
                if (numPoints < 256)
                {
                    // Decode indices as uint8T.
                    int[] face = new int[3];
                    for (int i = 0; i < numFaces; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                        {
                            byte val = buffer.DecodeU8();
                            face[j] = val;
                        }
                        Mesh.AddFace(face);
                    }
                }
                else if (numPoints < (1 << 16))
                {
                    // Decode indices as uint16T.
                    int[] face = new int[3];
                    for (int i = 0; i < numFaces; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                        {
                            ushort val = Buffer.DecodeU16();
                            face[j] = val;
                        }
                        Mesh.AddFace(face);
                    }
                }
                else if (Mesh.NumPoints < (1 << 21) && BitstreamVersion >= 22)
                {
                    // Decode indices as uint32_t.
                    int[] face = new int[3];
                    for (int i = 0; i < numFaces; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                        {
                            uint val = buffer.DecodeVarintU32();
                            face[j] = (int)val;
                        }
                        Mesh.AddFace(face);
                    }
                }
                else
                {
                    // Decode faces as uint (default).
                    int[] face = new int[3];
                    for (int i = 0; i < numFaces; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                        {
                            int val = Buffer.DecodeI32();
                            face[j] = val;
                        }
                        Mesh.AddFace(face);
                    }
                }
            }
            PointCloud.NumPoints = (int)numPoints;
        }

        protected override void CreateAttributesDecoder(int attrDecoderId)
        {

            // Always create the basic attribute decoder.
            SetAttributesDecoder(
                attrDecoderId,
                new SequentialAttributeDecodersController(
                    new LinearSequencer(PointCloud.NumPoints)));
        }

        /// <summary>
        /// Decodes face indices that were compressed with an entropy code.
        /// Returns false on error.
        /// </summary>
        void DecodeAndDecompressIndices(int numFaces)
        {

            // Get decoded indices differences that were encoded with an entropy code.
            Span<int> indicesBuffer = stackalloc int[numFaces * 3];
            Decoding.DecodeSymbols(numFaces * 3, 1, Buffer, indicesBuffer);
            // Reconstruct the indices from the differences.
            // See MeshSequentialEncoder::CompressAndEncodeIndices() for more details.
            int lastIndexValue = 0;
            int vertexIndex = 0;
            int[] face = new int[3];
            for (int i = 0; i < numFaces; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    int encodedVal = indicesBuffer[vertexIndex++];
                    int indexDiff = (encodedVal >> 1);
                    if ((encodedVal & 1) != 0)
                        indexDiff = -indexDiff;
                    int indexValue = indexDiff + lastIndexValue;
                    face[j] = indexValue;
                    lastIndexValue = indexValue;
                }
                Mesh.AddFace(face);
            }
        }

    }
}
