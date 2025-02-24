using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako
{
    class MeshTraversalSequencer<TCornerTable> : PointsSequencer where TCornerTable : ICornerTable
    {
        private ICornerTableTraverser<TCornerTable> traverser;
        private DracoMesh mesh;
        private MeshAttributeIndicesEncodingData encodingData;
        private IntList cornerOrder;

        public MeshTraversalSequencer(DracoMesh mesh,
            MeshAttributeIndicesEncodingData encodingData)
        {
            this.mesh = mesh;
            this.encodingData = encodingData;
        }

        public void SetTraverser(ICornerTableTraverser<TCornerTable> t)
        {
            traverser = t;
        }

        /// <summary>
        /// Function that can be used to set an order in which the mesh corners should
        /// be processed. This is an optional flag used usually only by the encoder
        /// to match the same corner order that is going to be used by the decoder.
        /// Note that |cornerOrder| should contain only one corner per face (it can
        /// have all corners but only the first encountered corner for each face is
        /// going to be used to start a traversal). If the corner order is not set, the
        /// corners are processed sequentially based on their ids.
        /// </summary>
        public void SetCornerOrder(IntList cornerOrder)
        {
            this.cornerOrder = cornerOrder;
        }

        public override void UpdatePointToAttributeIndexMapping(PointAttribute attribute)
        {
            TCornerTable cornerTable = traverser.CornerTable;
            attribute.SetExplicitMapping(mesh.NumPoints);
            int numFaces = mesh.NumFaces;
            int numPoints = mesh.NumPoints;
            Span<int> face = stackalloc int[3];
            for (int f = 0; f < numFaces; ++f)
            {
                mesh.ReadFace(f, face);
                for (int p = 0; p < 3; ++p)
                {
                    int pointId = face[p];
                    int vertId = cornerTable.Vertex(3*f + p);
                    int attEntryId = encodingData.vertexToEncodedAttributeValueIndexMap[vertId];
                    if (attEntryId >= numPoints)
                    {
                        // There cannot be more attribute values than the number of points.
                        throw DracoUtils.Failed();
                    }
                    attribute.SetPointMapEntry(pointId, attEntryId);
                }
            }
        }


        protected override void GenerateSequenceInternal()
        {
            traverser.OnTraversalStart();
            if (cornerOrder != null)
            {
                for (int i = 0; i < cornerOrder.Count; ++i)
                {
                    ProcessCorner(cornerOrder[i]);
                }
            }
            else
            {
                int num_faces = traverser.CornerTable.NumFaces;
                for (int i = 0; i < num_faces; ++i)
                {
                    ProcessCorner(3 * i);
                }
            }
            traverser.OnTraversalEnd();
        }

        private void ProcessCorner(int cornerId)
        {
            traverser.TraverseFromCorner(cornerId);
        }
    }
}
