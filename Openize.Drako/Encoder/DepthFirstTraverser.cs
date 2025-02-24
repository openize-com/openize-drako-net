using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Encoder
{
    internal class DepthFirstTraverser : TraverserBase<CornerTable>
    {
        private IntList corner_traversal_stack_ = new IntList();

        public DepthFirstTraverser(CornerTable cornerTable, MeshAttributeIndicesEncodingObserver<CornerTable> traversal_observer)
        {
            Init(cornerTable, traversal_observer);
        }

        public override void TraverseFromCorner(int corner_id)
        {

            if (this.IsCornerVisited(corner_id))
                return;  // Already traversed.

            corner_traversal_stack_.Clear();
            corner_traversal_stack_.Add(corner_id);
            // For the first face, check the remaining corners as they may not be
            // processed yet.
            int next_vert =
                this.corner_table_.Vertex(this.corner_table_.Next(corner_id));
            int prev_vert =
                this.corner_table_.Vertex(this.corner_table_.Previous(corner_id));
            if (next_vert == -1 || prev_vert == -1)
                throw DracoUtils.Failed();
            if (!this.IsVertexVisited(next_vert))
            {
                this.MarkVertexVisited(next_vert);
                this.traversal_observer_.OnNewVertexVisited(
                    next_vert, this.corner_table_.Next(corner_id));
            }
            if (!this.IsVertexVisited(prev_vert))
            {
                this.MarkVertexVisited(prev_vert);
                this.traversal_observer_.OnNewVertexVisited(
                    prev_vert, this.corner_table_.Previous(corner_id));
            }

            // Start the actual traversal.
            while (corner_traversal_stack_.Count > 0)
            {
                // Currently processed corner.
                corner_id = corner_traversal_stack_[corner_traversal_stack_.Count - 1];
                int face_id = corner_id / 3;
                // Make sure the face hasn't been visited yet.
                if (corner_id == -1 || this.IsFaceVisited(face_id))
                {
                    // This face has been already traversed.
                    corner_traversal_stack_.RemoveAt(corner_traversal_stack_.Count - 1);
                    continue;
                }
                while (true)
                {
                    this.MarkFaceVisited(face_id);
                    this.traversal_observer_.OnNewFaceVisited(face_id);
                    int vert_id = this.corner_table_.Vertex(corner_id);
                    if (vert_id == -1)
                        throw DracoUtils.Failed();
                    if (!this.IsVertexVisited(vert_id))
                    {
                        bool on_boundary = this.corner_table_.IsOnBoundary(vert_id);
                        this.MarkVertexVisited(vert_id);
                        this.traversal_observer_.OnNewVertexVisited(vert_id, corner_id);
                        if (!on_boundary)
                        {
                            corner_id = this.corner_table_.GetRightCorner(corner_id);
                            face_id = corner_id / 3;
                            continue;
                        }
                    }
                    // The current vertex has been already visited or it was on a boundary.
                    // We need to determine whether we can visit any of it's neighboring
                    // faces.
                    int right_corner_id =
                        this.corner_table_.GetRightCorner(corner_id);
                    int left_corner_id =
                        this.corner_table_.GetLeftCorner(corner_id);
                    int right_face_id = (right_corner_id == -1 ? -1 : right_corner_id / 3);
                    int left_face_id =
                        (left_corner_id == -1
                             ? -1
                             : left_corner_id / 3);
                    if (this.IsFaceVisited(right_face_id))
                    {
                        // Right face has been already visited.
                        if (this.IsFaceVisited(left_face_id))
                        {
                            // Both neighboring faces are visited. End reached.
                            corner_traversal_stack_.RemoveAt(corner_traversal_stack_.Count - 1);
                            break;  // Break from the while (true) loop.
                        }
                        else
                        {
                            // Go to the left face.
                            corner_id = left_corner_id;
                            face_id = left_face_id;
                        }
                    }
                    else
                    {
                        // Right face was not visited.
                        if (this.IsFaceVisited(left_face_id))
                        {
                            // Left face visited, go to the right one.
                            corner_id = right_corner_id;
                            face_id = right_face_id;
                        }
                        else
                        {
                            // Both neighboring faces are unvisited, we need to visit both of
                            // them.

                            // Split the traversal.
                            // First make the top of the current corner stack point to the left
                            // face (this one will be processed second).
                            corner_traversal_stack_[corner_traversal_stack_.Count - 1] = left_corner_id;
                            // Add a new corner to the top of the stack (right face needs to
                            // be traversed first).
                            corner_traversal_stack_.Add(right_corner_id);
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