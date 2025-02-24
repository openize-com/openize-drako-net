using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako
{
    /// <summary>
    /// Class providing the basic traversal funcionality needed by traversers (such
    /// as the EdgeBreakerTraverser, see edgebreakerTraverser.h). It is used to
    /// return the corner table that is used for the traversal, plus it provides a
    /// basic book-keeping of visited faces and vertices during the traversal.
    /// </summary>
    class CornerTableTraversalProcessor<TCornerTable> where TCornerTable : ICornerTable
    {
        private TCornerTable cornerTable;
        private bool[] isFaceVisited;
        private bool[] isVertexVisited;

        public CornerTableTraversalProcessor(TCornerTable cornerTable)
        {
            //Contract.Assert(cornerTable != null);
            this.cornerTable = cornerTable;
            isFaceVisited = new bool[cornerTable.NumFaces];
            ResetVertexData();
        }

        public bool IsFaceVisited(int faceId)
        {
            if (faceId < 0)
                return true; // Invalid faces are always considered as visited.
            return isFaceVisited[faceId];
        }

        public void MarkFaceVisited(int faceId)
        {
            isFaceVisited[faceId] = true;
        }

        public bool IsVertexVisited(int vertId)
        {
            return isVertexVisited[vertId];
        }

        public void MarkVertexVisited(int vertId)
        {
            isVertexVisited[vertId] = true;
        }

        protected virtual void ResetVertexData()
        {
            InitVertexData(cornerTable.NumVertices);
        }

        protected void InitVertexData(int numVerts)
        {
            isVertexVisited = new bool[numVerts];
        }

        public TCornerTable CornerTable
        {
            get { return cornerTable; }
        }
    }
}
