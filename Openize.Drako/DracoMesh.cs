using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;
using Openize.Drako.Decoder;

namespace Openize.Drako
{

    /// <summary>
    /// List of different variants of mesh attributes.
    /// </summary>
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    enum MeshAttributeElementType
    {
        /// <summary>
        /// All corners attached to a vertex share the same attribute value. A typical
        /// example are the vertex positions and often vertex colors.
        /// </summary>
        Vertex = 0,
        /// <summary>
        /// The most general attribute where every corner of the mesh can have a
        /// different attribute value. Often used for texture coordinates or normals.
        /// </summary>
        Corner,
        /// <summary>
        /// All corners of a single face share the same value.
        /// </summary>
        Face
    }

#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    class DracoMesh : DracoPointCloud
    {


        class AttributeData
        {
            public AttributeData()
            {
                elementType = MeshAttributeElementType.Corner;
            }

            internal MeshAttributeElementType elementType;
        }

        /// <summary>
        /// Mesh specific per-attribute data.
        /// </summary>
        private List<AttributeData> attributeData = new List<AttributeData>();
        private IntList faces = new IntList();

        public IntList Indices => faces;
        public void SetCorner(int corner, int value)
        {
            faces[corner] = value;
            //var tmp = faces[corner / 3];
            //tmp[corner % 3] = value;
        }

        public int ReadCorner(int corner)
        {
            return faces[corner];
            //var tmp = faces[corner / 3];
            //return tmp[corner % 3];
        }

        public void ReadFace(int faceId, int[] face)
        {
            var ptr = faceId * 3;
            face[0] = faces[ptr++];
            face[1] = faces[ptr++];
            face[2] = faces[ptr++];

        }
        public void ReadFace(int faceId, Span<int> face)
        {
            var ptr = faceId * 3;
            face[0] = faces[ptr++];
            face[1] = faces[ptr++];
            face[2] = faces[ptr++];

        }

        public void AddFace(int[] face)
        {
            faces.Add(face[0]);
            faces.Add(face[1]);
            faces.Add(face[2]);
        }

        public override int AddAttribute(PointAttribute pa)
        {
            var ad = new AttributeData();
            attributeData.Add(ad);
            return base.AddAttribute(pa);
        }

        public void SetFace(int faceId, int[] face)
        {
            if(faceId >= NumFaces)
                faces.Resize((faceId + 1) * 3);
            var p = faceId * 3;
            var data = faces.data;
            if(p + 2 < data.Length)
            {
                data[p + 0] = face[0];
                data[p + 1] = face[1];
                data[p + 2] = face[2];
            }

            //if(faceId >= NumFaces)
            //    A3DUtils.Resize(NumFaces * 3, faceId + 1);
            //faces[faceId] = face;
        }

        public int NumFaces
        {
            get { return faces.Count / 3; }
            set {
                faces.Resize(value * 3);
                    //A3DUtils.Resize(faces, value);
            }
        }

        internal override void ApplyPointIdDeduplication(int[] idMap, IntList uniquePointIds)
        {
            base.ApplyPointIdDeduplication(idMap, uniquePointIds);

            int p = 0;
            for (int f = 0; f < NumFaces; ++f)
            {
                for (int c = 0; c < 3; ++c)
                {
                    faces[p] = idMap[faces[p]];
                    p++;
                }
            }
        }

        private static MeshDecoder CreateMeshDecoder(DracoEncodingMethod method)
        {
            if (method == DracoEncodingMethod.Sequential)
                return new MeshSequentialDecoder();
            if(method == DracoEncodingMethod.EdgeBreaker)
                return new MeshEdgeBreakerDecoder();
            return null;
        }
        private static PointCloudDecoder CreatePointCloudDecoder(DracoEncodingMethod method)
        {
            if (method == DracoEncodingMethod.Sequential)
                return new PointCloudSequentialDecoder();
            if(method == DracoEncodingMethod.KdTree)
                return new PointCloudKdTreeDecoder();
            return null;
        }
        internal static DracoPointCloud Decode(DecoderBuffer buffer, bool decodeData = true)
        {
            DracoHeader header = DracoHeader.Parse(buffer);
            if (header == null)
                return null;
            if (header.encoderType == EncodedGeometryType.TriangularMesh)
            {
                return DecodeMesh(buffer, header, decodeData);
            }
            else if (header.encoderType == EncodedGeometryType.PointCloud)
            {
                return DecodePointCloud(buffer, header, decodeData);
            }
            return null;
        }

        private static DracoPointCloud DecodePointCloud(DecoderBuffer buffer, DracoHeader header, bool decodeData)
        {
            buffer.BitstreamVersion = header.version;
            PointCloudDecoder decoder = CreatePointCloudDecoder(header.method);
            if (decoder == null)
                return null;
            try
            {
                DracoPointCloud ret = new DracoPointCloud();
                decoder.Decode(header, buffer, ret, decodeData);
                return ret;
            }
            catch(Exception)
            {
                return null;
            }
        }
        private static DracoMesh DecodeMesh(DecoderBuffer buffer, DracoHeader header, bool decodeData)
        {
            buffer.BitstreamVersion = header.version;
            MeshDecoder decoder = CreateMeshDecoder(header.method);
            if (decoder == null)
                return null;
            try
            {
                DracoMesh ret = new DracoMesh();
                decoder.Decode(header, buffer, ret, decodeData);
                return ret;
            }
            catch(Exception)
            {
                return null;
            }
        }

        public MeshAttributeElementType GetAttributeElementType(int attId)
        {
            return attributeData[attId].elementType;
        }
    }
}
