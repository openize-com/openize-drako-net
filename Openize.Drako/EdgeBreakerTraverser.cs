using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Encoder;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako
{
    class EdgeBreakerTraverser<TCornerTable> : TraverserBase<TCornerTable> where TCornerTable:ICornerTable
    {
        //private CornerTableTraversalProcessor<TCornerTable> processor;
        private MeshAttributeIndicesEncodingObserver<TCornerTable> traversalObserver;
        private IntList cornerTraversalStack = new IntList();


        public EdgeBreakerTraverser(TCornerTable corner_table,
            MeshAttributeIndicesEncodingObserver<TCornerTable> observer)
        {

            this.traversalObserver = observer;
            Init(corner_table, observer);
        }

        public EdgeBreakerTraverser(MeshAttributeIndicesEncodingObserver<TCornerTable> traversalObserver)
        {
            this.corner_table_ = traversalObserver.CornerTable;
            this.traversalObserver = traversalObserver;
            Init(corner_table_, traversalObserver);
        }

        public override void TraverseFromCorner(int cornerId)
        {
            cornerTraversalStack.Clear();
            cornerTraversalStack.Add(cornerId);
            // For the first face, check the remaining corners as they may not be
            // processed yet.
            int nextVert = corner_table_.Vertex(corner_table_.Next(cornerId));
            int prevVert = corner_table_.Vertex(corner_table_.Previous(cornerId));
            if (!IsVertexVisited(nextVert))
            {
                MarkVertexVisited(nextVert);
                traversalObserver.OnNewVertexVisited(nextVert,corner_table_.Next(cornerId));
            }
            if (!IsVertexVisited(prevVert))
            {
                MarkVertexVisited(prevVert);
                traversalObserver.OnNewVertexVisited(prevVert, corner_table_.Previous(cornerId));
            }

            // Start the actual traversal.
            while (cornerTraversalStack.Count > 0)
            {
                // Currently processed corner.
                cornerId = cornerTraversalStack[cornerTraversalStack.Count - 1];
                int faceId  = cornerId/3;
                // Make sure the face hasn't been visited yet.
                if (cornerId < 0 || IsFaceVisited(faceId))
                {
                    // This face has been already traversed.
                    cornerTraversalStack.RemoveAt(cornerTraversalStack.Count - 1);
                    continue;
                }
                while (true)
                {
                    faceId = cornerId/3;
                    MarkFaceVisited(faceId);
                    traversalObserver.OnNewFaceVisited(faceId);
                    int vertId = corner_table_.Vertex(cornerId);
                    bool onBoundary = corner_table_.IsOnBoundary(vertId);
                    if (!IsVertexVisited(vertId))
                    {
                        MarkVertexVisited(vertId);
                        traversalObserver.OnNewVertexVisited(vertId, cornerId);
                        if (!onBoundary)
                        {
                            cornerId = corner_table_.GetRightCorner(cornerId);
                            continue;
                        }
                    }
                    // The current vertex has been already visited or it was on a boundary.
                    // We need to determine whether we can visit any of it's neighboring
                    // faces.
                    int rightCornerId =
                        corner_table_.GetRightCorner(cornerId);
                    int leftCornerId =
                        corner_table_.GetLeftCorner(cornerId);
                    int rightFaceId = rightCornerId < 0 ? -1 : rightCornerId/3;
                    int leftFaceId = leftCornerId < 0 ? -1 : leftCornerId/3;
                    if (IsFaceVisited(rightFaceId))
                    {
                        // Right face has been already visited.
                        if (IsFaceVisited(leftFaceId))
                        {
                            // Both neighboring faces are visited. End reached.
                            cornerTraversalStack.RemoveAt(cornerTraversalStack.Count - 1);
                            break; // Break from the while (true) loop.
                        }
                        else
                        {
                            // Go to the left face.
                            cornerId = leftCornerId;
                        }
                    }
                    else
                    {
                        // Right face was not visited.
                        if (IsFaceVisited(leftFaceId))
                        {
                            // Left face visited, go to the right one.
                            cornerId = rightCornerId;
                        }
                        else
                        {
                            // Both neighboring faces are unvisited, we need to visit both of
                            // them.

                            // Split the traversal.
                            // First make the top of the current corner stack point to the left
                            // face (this one will be processed second).
                            cornerTraversalStack[cornerTraversalStack.Count - 1] = leftCornerId;
                            // Add a new corner to the top of the stack (right face needs to
                            // be traversed first).
                            cornerTraversalStack.Add(rightCornerId);
                            // Break from the while (true) loop.
                            break;
                        }
                    }
                }
            }
        }

        public override void OnTraversalStart()
        {
        }

        public override void OnTraversalEnd()
        {
        }
    }
}
