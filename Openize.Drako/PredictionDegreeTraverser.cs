using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Encoder;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako
{
    class PredictionDegreeTraverser : TraverserBase<CornerTable>
    {

        // For efficiency reasons, the priority traversal is implemented using buckets
        // where each buckets represent a stack of available corners for a given
        // priority. Corners with the highest priority are always processed first.
        const int kMaxPriority = 3;
        private Stack<int>[] traversal_stacks_ = new Stack<int>[kMaxPriority];

        // Keep the track of the best available priority to improve the performance
        // of PopNextCornerToTraverse() method.
        int best_priority_;

        // Prediction degree available for each vertex.
        private int[] prediction_degree_;
        CornerTable cornerTable;
        private MeshAttributeIndicesEncodingObserver<CornerTable> observer;

        public PredictionDegreeTraverser(MeshAttributeIndicesEncodingObserver<CornerTable> observer)
        {
            this.cornerTable = observer.CornerTable;
            this.observer = observer;
            Init(cornerTable, observer);
            for (int i = 0; i < kMaxPriority; i++)
                traversal_stacks_[i] = new Stack<int>();
        }


        // Called before any traversing starts.
        public override void OnTraversalStart()
        {
            prediction_degree_ = new int[cornerTable.NumVertices];
        }

        // Called when all the traversing is done.
        public override void OnTraversalEnd()
        {
        }

        public override void TraverseFromCorner(int corner_id)
        {
            if (prediction_degree_.Length == 0)
                return;

            // Traversal starts from the |corner_id|. It's going to follow either the
            // right or the left neighboring faces to |corner_id| based on their
            // prediction degree.
            traversal_stacks_[0].Push(corner_id);
            best_priority_ = 0;
            // For the first face, check the remaining corners as they may not be
            // processed yet.
            int next_vert =
                cornerTable.Vertex(cornerTable.Next(corner_id));
            int prev_vert =
                cornerTable.Vertex(cornerTable.Previous(corner_id));
            if (!IsVertexVisited(next_vert))
            {
                MarkVertexVisited(next_vert);
                observer.OnNewVertexVisited(next_vert,
                    cornerTable.Next(corner_id));
            }

            if (!IsVertexVisited(prev_vert))
            {
                MarkVertexVisited(prev_vert);
                observer.OnNewVertexVisited(
                    prev_vert, cornerTable.Previous(corner_id));
            }

            int tip_vertex = cornerTable.Vertex(corner_id);
            if (!IsVertexVisited(tip_vertex))
            {
                MarkVertexVisited(tip_vertex);
                observer.OnNewVertexVisited(tip_vertex, corner_id);
            }

            // Start the actual traversal.
            while ((corner_id = PopNextCornerToTraverse()) != CornerTable.kInvalidCornerIndex)
            {
                int face_id = corner_id / 3;
                // Make sure the face hasn't been visited yet.
                if (IsFaceVisited(face_id))
                {
                    // This face has been already traversed.
                    continue;
                }

                while (true)
                {
                    face_id = corner_id / 3;
                    MarkFaceVisited(face_id);
                    observer.OnNewFaceVisited(face_id);

                    // If the newly reached vertex hasn't been visited, mark it and notify
                    // the observer.
                    int vert_id = cornerTable.Vertex(corner_id);
                    if (!IsVertexVisited(vert_id))
                    {
                        MarkVertexVisited(vert_id);
                        observer.OnNewVertexVisited(vert_id, corner_id);
                    }

                    // Check whether we can traverse to the right and left neighboring
                    // faces.
                    int right_corner_id =
                        cornerTable.GetRightCorner(corner_id);
                    int left_corner_id =
                        cornerTable.GetLeftCorner(corner_id);
                    int right_face_id =
                        right_corner_id == CornerTable.kInvalidCornerIndex
                            ? CornerTable.kInvalidFaceIndex
                            : right_corner_id / 3;
                    int left_face_id = left_corner_id == CornerTable.kInvalidCornerIndex
                        ? CornerTable.kInvalidFaceIndex
                        : left_corner_id / 3;
                    bool is_right_face_visited =
                        IsFaceVisited(right_face_id);
                    bool is_left_face_visited =
                        IsFaceVisited(left_face_id);

                    if (!is_left_face_visited)
                    {
                        // We can go to the left face.
                        int priority = ComputePriority(left_corner_id);
                        if (is_right_face_visited && priority <= best_priority_)
                        {
                            // Right face has been already visited and the priority is equal or
                            // better than the best priority. We are sure that the left face
                            // would be traversed next so there is no need to put it onto the
                            // stack.
                            corner_id = left_corner_id;
                            continue;
                        }
                        else
                        {
                            AddCornerToTraversalStack(left_corner_id, priority);
                        }
                    }

                    if (!is_right_face_visited)
                    {
                        // Go to the right face.
                        int priority = ComputePriority(right_corner_id);
                        if (priority <= best_priority_)
                        {
                            // We are sure that the right face would be traversed next so there
                            // is no need to put it onto the stack.
                            corner_id = right_corner_id;
                            continue;
                        }
                        else
                        {
                            AddCornerToTraversalStack(right_corner_id, priority);
                        }
                    }

                    // Couldn't proceed directly to the next corner
                    break;
                }
            }
        }


        // Retrieves the next available corner (edge) to traverse. Edges are processed
        // based on their priorities.
        // Returns kInvalidCornerIndex when there is no edge available.
        private int PopNextCornerToTraverse()
        {
            for (int i = best_priority_; i < kMaxPriority; ++i)
            {
                if (traversal_stacks_[i].Count > 0)
                {
                    int ret = traversal_stacks_[i].Peek();
                    traversal_stacks_[i].Pop();
                    best_priority_ = i;
                    return ret;
                }
            }

            return CornerTable.kInvalidCornerIndex;
        }

        private void AddCornerToTraversalStack(int ci, int priority)
        {
            traversal_stacks_[priority].Push(ci);
            // Make sure that the best available priority is up to date.
            if (priority < best_priority_)
                best_priority_ = priority;
        }

        // Returns the priority of traversing edge leading to |corner_id|.
        private int ComputePriority(int corner_id)
        {
            int v_tip = cornerTable.Vertex(corner_id);
            // Priority 0 when traversing to already visited vertices.
            int priority = 0;
            if (!IsVertexVisited(v_tip))
            {
                int degree = ++prediction_degree_[v_tip];
                // Priority 1 when prediction degree > 1, otherwise 2.
                priority = (degree > 1 ? 1 : 2);
            }

            // Clamp the priority to the maximum number of buckets.
            if (priority >= kMaxPriority)
                priority = kMaxPriority - 1;
            return priority;
        }

    }
}
