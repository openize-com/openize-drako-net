using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako
{
    class MeshAttributeIndicesEncodingObserver<TCornerTable> where TCornerTable : ICornerTable
    {
        private TCornerTable attConnectivity;
        private MeshAttributeIndicesEncodingData encodingData;
        private DracoMesh mesh;
        private PointsSequencer sequencer;

        public MeshAttributeIndicesEncodingObserver(TCornerTable cornerTable, DracoMesh mesh, PointsSequencer sequencer,
            MeshAttributeIndicesEncodingData encodingData)
        {
            this.encodingData = encodingData;
            this.mesh = mesh;
            this.attConnectivity = cornerTable;
            this.sequencer = sequencer;
        }

        public TCornerTable CornerTable
        {
            get { return attConnectivity; }
        }

        public void OnNewFaceVisited(int face)
        {
        }

        public void OnNewVertexVisited(int vertex, int corner)
        {
            //int pointId = mesh.Face(corner/3)[corner%3];
            int pointId = mesh.ReadCorner(corner);
            // Append the visited attribute to the encoding order.
            sequencer.AddPointId(pointId);

            // Keep track of visited corners.
            encodingData.encodedAttributeValueIndexToCornerMap.Add(
                corner);

            encodingData.vertexToEncodedAttributeValueIndexMap[vertex] =
                encodingData.numValues;

            encodingData.numValues++;
        }

    }
}
