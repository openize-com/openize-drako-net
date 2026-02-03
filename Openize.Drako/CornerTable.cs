using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;
using Openize.Drako.Compression;

namespace Openize.Drako
{
    /// <summary>
    /// CornerTable is used to represent connectivity of triangular meshes.
    /// For every corner of all faces, the corner table stores the index of the
    /// opposite corner in the neighboring face (if it exists) as illustrated in the
    /// figure below (see corner |c| and it's opposite corner |o|).
    ///
    ///     *
    ///    /c\
    ///   /   \
    ///  /n   p\
    /// *-------*
    ///  \     /
    ///   \   /
    ///    \o/
    ///     *
    ///
    /// All corners are defined by unique CornerIndex and each triplet of corners
    /// that define a single face id always ordered consecutively as:
    ///     { 3 * FaceIndex, 3 * FaceIndex + 1, 3 * FaceIndex +2 }.
    /// This representation of corners allows CornerTable to easily retrieve Next and
    /// Previous corners on any face (see corners |n| and |p| in the figure above).
    /// Using the Next, Previous, and Opposite corners then enables traversal of any
    /// 2-manifold surface.
    /// If the CornerTable is constructed from a non-manifold surface, the input
    /// non-manifold edges and vertices are automatically split.
    /// </summary>
    class CornerTable : ICornerTable
    {
        public const int kInvalidFaceIndex = -1;
        public const int kInvalidCornerIndex = -1;
        public const int kInvalidVertexIndex = -1;
        private int[] oppositeCorners;
        private IntList vertexCorners = new IntList();
        private int[] cornerToVertexMap;

        private int numOriginalVertices;
        private int numDegeneratedFaces;
        private int numIsolatedVertices;
        private IntList nonManifoldVertexParents = new IntList();
        private ValenceCache valenceCache;

        struct VertexEdgePair
        {
            public int sinkVert;
            public int edgeCorner;

            public VertexEdgePair(int sinkVert, int edgeCorner)
            {
                this.sinkVert = sinkVert;
                this.edgeCorner = edgeCorner;
            }
        }

        public CornerTable()
        {
            this.valenceCache = new ValenceCache(this);
        }

        public void Initialize(int[,] faces)
        {

            valenceCache.ClearValenceCache();
            valenceCache.ClearValenceCacheInaccurate();
            int numFaces = faces.GetLength(0);
            cornerToVertexMap = new int[numFaces * 3];
            for (int fi = 0; fi < numFaces; ++fi)
            {
                for (int i = 0; i < 3; ++i)
                {
                    int corner = FirstCorner(fi);
                    cornerToVertexMap[corner + i] = faces[fi, i];
                }
            }

            int numVertices = -1;
            ComputeOppositeCorners(out numVertices);
            ComputeVertexCorners(numVertices);
        }

        public int[] AllCorners(int face)
        {
            int ci = face * 3;
            return new int[] {ci, ci + 1, ci + 2};
        }


        public override int Opposite(int corner)
        {
            if (corner < 0)
                return corner;
            return oppositeCorners[corner];
        }

        public override int Vertex(int corner)
        {
            if (corner < 0)
                return kInvalidVertexIndex;
            return cornerToVertexMap[corner];
        }

        public int Face(int corner)
        {
            if (corner < 0)
                return kInvalidFaceIndex;
            return corner / 3;
        }

        public int FirstCorner(int face)
        {
            if (face < 0)
                return kInvalidCornerIndex;
            return face * 3;
        }

        /// <summary>
        ///     *-------*
        ///    / \     / \
        ///   /   \   /   \
        ///  /   sl\c/sr   \
        /// *-------v-------*
        /// Returns the corner on the adjacent face on the right that maps to
        /// the same vertex as the given corner (sr in the above diagram).
        /// </summary>
        public override int SwingRight(int corner)
        {
            return Previous(Opposite(Previous(corner)));
        }


        /// <summary>
        /// Returns the number of new vertices that were created as a result of
        /// spliting of non-manifold vertices of the input geometry.
        /// </summary>
        public int NumNewVertices
        {
            get { return NumVertices - numOriginalVertices; }
        }

        public int NumOriginalVertices
        {
            get { return numOriginalVertices; }
        }

        /// <summary>
        /// Returns the number of faces with duplicated vertex indices.
        /// </summary>
        public int NumDegeneratedFaces
        {
            get { return numDegeneratedFaces; }
        }

        /// <summary>
        /// Returns the number of isolated vertices (vertices that have
        /// vertexCorners mapping set to kInvalidCornerIndex.
        /// </summary>
        public int NumIsolatedVertices
        {
            get { return numIsolatedVertices; }
        }

        public override int NumVertices
        {
            get { return vertexCorners.Count; }
        }

        public override int NumCorners
        {
            get { return cornerToVertexMap.Length; }
        }

        public override int NumFaces
        {
            get { return cornerToVertexMap.Length / 3; }
        }


        /// <summary>
        /// Returns the corner on the left face that maps to the same vertex as the
        /// given corner (sl in the above diagram).
        /// </summary>
        public override int SwingLeft(int corner)
        {
            return Next(Opposite(Next(corner)));
        }

        private void ComputeOppositeCorners(out int numVertices)
        {
            oppositeCorners = new int[NumCorners];
            for (int i = 0; i < oppositeCorners.Length; i++)
                oppositeCorners[i] = kInvalidCornerIndex;

            // Out implementation for finding opposite corners is based on keeping track
            // of outgoing half-edges for each vertex of the mesh. Half-edges (defined by
            // their opposite corners) are processed one by one and whenever a new
            // half-edge (corner) is processed, we check whether the sink vertex of
            // this half-edge contains its sibling half-edge. If yes, we connect them and
            // remove the sibling half-edge from the sink vertex, otherwise we add the new
            // half-edge to its source vertex.

            // First compute the number of outgoing half-edges (corners) attached to each
            // vertex.
            IntList numCornersOnVertices = new IntList();
            numCornersOnVertices.Capacity = NumCorners;
            for (int c = 0; c < NumCorners; ++c)
            {
                int v1 = Vertex(c);
                if (v1 >= numCornersOnVertices.Count)
                    numCornersOnVertices.Resize(v1 + 1, 0);
                // For each corner there is always exactly one outgoing half-edge attached
                // to its vertex.
                numCornersOnVertices[v1]++;
            }

            // Create a storage for half-edges on each vertex. We store all half-edges in
            // one array, where each entry is identified by the half-edge's sink vertex id
            // and the associated half-edge corner id (corner opposite to the half-edge).
            // Each vertex will be assigned storage for up to
            // |numCornersOnVertices[vertId]| half-edges. Unused half-edges are marked
            // with |sinkVert| == -1.
            VertexEdgePair[] vertexEdges = new VertexEdgePair[NumCorners];
            for(int i = 0; i < vertexEdges.Length; i++)
                vertexEdges[i] = new VertexEdgePair(-1, -1);

            // For each vertex compute the offset (location where the first half-edge
            // entry of a given vertex is going to be stored). This way each vertex is
            // guaranteed to have a non-overlapping storage with respect to the other
            // vertices.
            int[] vertexOffset = new int[numCornersOnVertices.Count];
            int offset = 0;
            for (int i = 0; i < numCornersOnVertices.Count; ++i)
            {
                vertexOffset[i] = offset;
                offset += numCornersOnVertices[i];
            }

            // Now go over the all half-edges (using their opposite corners) and either
            // insert them to the |vertexEdge| array or connect them with existing
            // half-edges.
            for (int c = 0; c < NumCorners; ++c)
            {
                int sourceV = Vertex(Next(c));
                int sinkV = Vertex(Previous(c));

                int faceIndex = Face(c);
                if (c == FirstCorner(faceIndex))
                {
                    // Check whether the face is degenerated, if so ignore it.
                    int v0 = Vertex(c);
                    if (v0 == sourceV || v0 == sinkV || sourceV == sinkV)
                    {
                        ++numDegeneratedFaces;
                        c += 2; // Ignore the next two corners of the same face.
                        continue;
                    }
                }

                int oppositeC = -1;
                // The maximum number of half-edges attached to the sink vertex.
                int numCornersOnVert = numCornersOnVertices[sinkV];
                // Where to look for the first half-edge on the sink vertex.
                offset = vertexOffset[sinkV];
                for (int i = 0; i < numCornersOnVert; ++i, ++offset)
                {
                    int otherV = vertexEdges[offset].sinkVert;
                    if (otherV < 0)
                        break; // No matching half-edge found on the sink vertex.
                    if (otherV == sourceV)
                    {
                        // A matching half-edge was found on the sink vertex. Mark the
                        // half-edge's opposite corner.
                        oppositeC = vertexEdges[offset].edgeCorner;
                        // Remove the half-edge from the sink vertex. We remap all subsequent
                        // half-edges one slot down.
                        // TODO(ostava): This can be optimized a little bit, by remaping only
                        // the half-edge on the last valid slot into the deleted half-edge's
                        // slot.
                        for (int j = i + 1; j < numCornersOnVert; ++j, ++offset)
                        {
                            vertexEdges[offset] = vertexEdges[offset + 1];
                            if (vertexEdges[offset].sinkVert < 0)
                                break; // Unused half-edge reached.
                        }
                        // Mark the last entry as unused.
                        vertexEdges[offset].sinkVert = -1;
                        break;
                    }
                }
                if (oppositeC < 0)
                {
                    // No opposite corner found. Insert the new edge
                    int numCornersOnSourceVert = numCornersOnVertices[sourceV];
                    offset = vertexOffset[sourceV];
                    for (int i = 0; i < numCornersOnSourceVert; ++i, ++offset)
                    {
                        // Find the first unused half-edge slot on the source vertex.
                        if (vertexEdges[offset].sinkVert < 0)
                        {
                            vertexEdges[offset].sinkVert = sinkV;
                            vertexEdges[offset].edgeCorner = c;
                            break;
                        }
                    }
                }
                else
                {
                    // Opposite corner found.
                    oppositeCorners[c] = oppositeC;
                    oppositeCorners[oppositeC] = c;
                }
            }
            numVertices = numCornersOnVertices.Count;
        }

        void ComputeVertexCorners(int numVertices)
        {
            numOriginalVertices = numVertices;
            vertexCorners.Resize(numVertices, kInvalidCornerIndex);
            // Arrays for marking visited vertices and corners that allow us to detect
            // non-manifold vertices.
            bool[] visitedVertices = new bool[numVertices];
            int numVisitedVertices = numVertices;

            bool[] visitedCorners = new bool[NumCorners];

            for (int f = 0; f < NumFaces; ++f)
            {
                int firstFaceCorner = FirstCorner(f);
                // Check whether the face is degenerated. If so ignore it.
                if (IsDegenerated(f))
                    continue;

                for (int k = 0; k < 3; ++k)
                {
                    int c = firstFaceCorner + k;
                    if (visitedCorners[c])
                        continue;
                    int v = cornerToVertexMap[c];
                    // Note that one vertex maps to many corners, but we just keep track
                    // of the vertex which has a boundary on the left if the vertex lies on
                    // the boundary. This means that all the related corners can be accessed
                    // by iterating over the SwingRight() operator.
                    // In case of a vertex inside the mesh, the choice is arbitrary.
                    bool isNonManifoldVertex = false;
                    if (visitedVertices[v])
                    {
                        // A visited vertex of an unvisited corner found. Must be a non-manifold
                        // vertex.
                        // Create a new vertex for it.
                        vertexCorners.Add(kInvalidCornerIndex);
                        nonManifoldVertexParents.Add(v);
                        if(numVisitedVertices >= visitedVertices.Length)
                        {
                            //resize 
                            Array.Resize(ref visitedVertices, visitedVertices.Length * 2);
                        }
                        visitedVertices[numVisitedVertices++] = false;
                        v = numVertices++;
                        isNonManifoldVertex = true;
                    }
                    // Mark the vertex as visited.
                    visitedVertices[v] = true;

                    // First swing all the way to the left and mark all corners on the way.
                    int actC = c;
                    while (actC != kInvalidCornerIndex)
                    {
                        visitedCorners[actC] = true;
                        // Vertex will eventually point to the left most corner.
                        vertexCorners[v] = actC;
                        if (isNonManifoldVertex)
                        {
                            // Update vertex index in the corresponding face.
                            cornerToVertexMap[actC] = v;
                        }
                        actC = SwingLeft(actC);
                        if (actC == c)
                            break; // Full circle reached.
                    }
                    if (actC == kInvalidCornerIndex)
                    {
                        // If we have reached an open boundary we need to swing right from the
                        // initial corner to mark all corners in the opposite direction.
                        actC = SwingRight(c);
                        while (actC != kInvalidCornerIndex)
                        {
                            visitedCorners[actC] = true;
                            if (isNonManifoldVertex)
                            {
                                // Update vertex index in the corresponding face.
                                int actF = Face(actC);
                                cornerToVertexMap[actC] = v;
                            }
                            actC = SwingRight(actC);
                        }
                    }
                }
            }

            // Count the number of isolated (unprocessed) vertices.
            numIsolatedVertices = 0;
            for (int i = 0; i < numVisitedVertices; i++)
            {
                if (!visitedVertices[i])
                    ++numIsolatedVertices;
            }
        }

        public bool IsDegenerated(int face)
        {
            if (face == kInvalidFaceIndex)
                return true;
            int firstFaceCorner = FirstCorner(face);
            int v0 = Vertex(firstFaceCorner);
            int v1 = Vertex(Next(firstFaceCorner));
            int v2 = Vertex(Previous(firstFaceCorner));
            if (v0 == v1 || v0 == v2 || v1 == v2)
                return true;
            return false;
        }

        /// <summary>
        /// Returns the left-most corner of a single vertex 1-ring. If a vertex is not
        /// on a boundary (in which case it has a full 1-ring), this function returns
        /// any of the corners mapped to the given vertex.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public override int LeftMostCorner(int v)
        {
            return vertexCorners[v];
        }

        /// <summary>
        /// Returns true if the specified vertex is on a boundary.
        /// </summary>
        /// <param name="vert"></param>
        /// <returns></returns>
        public override bool IsOnBoundary(int vert)
        {
            int corner = LeftMostCorner(vert);
            if (SwingLeft(corner) < 0)
                return true;
            return false;
        }

        /// <summary>
        /// Get opposite corners on the left and right faces respecitively (see image
        /// below, where L and R are the left and right corners of a corner X.
        ///
        /// *-------*-------*
        ///  \L    /X\    R/
        ///   \   /   \   /
        ///    \ /     \ /
        ///     *-------*
        /// </summary>
        /// <param name="cornerId"></param>
        /// <returns></returns>
        public override int GetLeftCorner(int cornerId)
        {
            if (cornerId < 0)
                return kInvalidCornerIndex;
            return Opposite(Previous(cornerId));
        }

        public override int GetRightCorner(int cornerId)
        {
            if (cornerId < 0)
                return kInvalidCornerIndex;
            return Opposite(Next(cornerId));
        }

        /// <summary>
        /// Methods that modify an existing corner table.
        /// Sets the opposite corner mapping between two corners. Caller must ensure
        /// that the indices are valid.
        /// </summary>
        /// <param name="cornerId"></param>
        /// <param name="oppCornerId"></param>
        public void SetOppositeCorner(int cornerId, int oppCornerId)
        {
            oppositeCorners[cornerId] = oppCornerId;
        }

        /// <summary>
        /// Updates mapping betweeh a corner and a vertex.
        /// </summary>
        /// <param name="cornerId"></param>
        /// <param name="vertId"></param>
        public void MapCornerToVertex(int cornerId, int vertId)
        {
            if (vertId >= 0)
            {
                if (vertexCorners.Count <= vertId)
                    vertexCorners.Resize(vertId + 1);
                cornerToVertexMap[cornerId] = vertId;
            }
        }

        /// <summary>
        /// Sets a new left most corner for a given vertex.
        /// </summary>
        /// <param name="vert"></param>
        /// <param name="corner"></param>
        public void SetLeftMostCorner(int vert, int corner)
        {
            if (vert != kInvalidVertexIndex)
                vertexCorners[vert] = corner;
        }

        /// <summary>
        /// Makes a vertex isolated (not attached to any corner).
        /// </summary>
        /// <param name="vert"></param>
        public void MakeVertexIsolated(int vert)
        {
            vertexCorners[vert] = kInvalidCornerIndex;
        }

        /// <summary>
        /// Returns true if a vertex is not attached to any face.
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public bool IsVertexIsolated(int v)
        {
            return LeftMostCorner(v) < 0;
        }

        /// <summary>
        /// Updates the vertex to corner map on a specified vertex. This should be
        /// called in cases where the mapping may be invalid (e.g. when the corner
        /// table was constructed manually).
        /// </summary>
        public void UpdateVertexToCornerMap(int vert)
        {
            int firstC = vertexCorners[vert];
            if (firstC < 0)
                return; // Isolated vertex.
            int actC = SwingLeft(firstC);
            int c = firstC;
            while (actC >= 0 && actC != firstC)
            {
                c = actC;
                actC = SwingLeft(actC);
            }
            if (actC != firstC)
            {
                vertexCorners[vert] = c;
            }
        }

        // Resets the corner table to the given number of invalid faces.
        public void Reset(int numFaces, int numVertices)
        {
            if (numFaces < 0 || numVertices < 0)
                throw new ArgumentException();
            if (numFaces > int.MaxValue / 3)
                throw new ArgumentException();
            cornerToVertexMap = new int[numFaces * 3];
            oppositeCorners = new int[numFaces * 3];
            for (int i = 0; i < cornerToVertexMap.Length; i++)
            {
                cornerToVertexMap[i] = -1;
                oppositeCorners[i] = -1;
            }
            vertexCorners.Capacity = numVertices;
            valenceCache.ClearValenceCache();
            valenceCache.ClearValenceCacheInaccurate();
        }

        public int ConfidentVertex(int corner)
        {
            //DRACO_DCHECK_GE(corner.value(), 0);
            //DRACO_DCHECK_LT(corner.value(), num_corners());
            return cornerToVertexMap[corner];
        }
        public int Valence(int v)
        {
            if (v == kInvalidVertexIndex)
                return -1;
            // iterating over vertices in a 1-ring around the specified vertex.
            int valence = 0;
            /*
            VertexRingIterator vi(this , v );
            for (; !vi.End(); vi.Next())
            {
                ++valence;
            }
            */
            int startCorner = LeftMostCorner(v);
            int corner = startCorner;
            bool leftTraversal = true;

            while (corner >= 0)
            {
                valence++;

                if (leftTraversal)
                {
                    corner = SwingLeft(corner);
                    if (corner < 0)
                    {
                        // Open boundary reached.
                        corner = startCorner;
                        leftTraversal = false;
                    }
                    else if (corner == startCorner)
                    {
                        // End reached.
                        corner = kInvalidCornerIndex;
                    }
                }
                else
                {
                    // Go to the right until we reach a boundary there (no explicit check
                    // is needed in this case).
                    corner = SwingRight(corner);
                }
            }
            return valence;
        }

        public int AddNewVertex()
        {
            vertexCorners.Add(kInvalidCornerIndex);
            return vertexCorners.Count - 1;
        }
    }
}
