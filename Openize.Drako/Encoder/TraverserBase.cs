namespace FileFormat.Drako.Encoder
{
    abstract class TraverserBase<TCornerTable> : ICornerTableTraverser<TCornerTable> where TCornerTable:ICornerTable
    {
        protected TCornerTable corner_table_;
        protected bool[] is_face_visited_;
        protected bool[] is_vertex_visited_;
        protected MeshAttributeIndicesEncodingObserver<TCornerTable> traversal_observer_;

        public void Init(TCornerTable corner_table, MeshAttributeIndicesEncodingObserver<TCornerTable> traversal_observer)
        {
            corner_table_ = corner_table;
            is_face_visited_ = new bool[corner_table_.NumFaces];
            is_vertex_visited_ = new bool[corner_table_.NumVertices];
            this.traversal_observer_ = traversal_observer;
        }
        public TCornerTable CornerTable { get => corner_table_; }

        protected bool IsFaceVisited(int face_id)
        {
            if (face_id == -1)
                return true;  // Invalid faces are always considered as visited.
            return is_face_visited_[face_id];
        }

        // Returns true if the face containing the given corner was visited.
        protected bool IsCornerVisited(int corner_id)
        {
            if (corner_id == -1)
                return true;  // Invalid faces are always considered as visited.
            return is_face_visited_[corner_id / 3];
        }

        protected void MarkFaceVisited(int face_id)
        {
            is_face_visited_[face_id] = true;
        }
        protected bool IsVertexVisited(int vert_id)
        {
            return is_vertex_visited_[vert_id];
        }
        protected void MarkVertexVisited(int vert_id)
        {
            is_vertex_visited_[vert_id] = true;
        }

        public abstract void TraverseFromCorner(int cornerId);

        public abstract void OnTraversalStart();

        public abstract void OnTraversalEnd();
    }
}