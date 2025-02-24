using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;
using FileFormat.Drako.Compression;

namespace FileFormat.Drako.Encoder
{

    /// <summary>
    /// Class implementing the edgebreaker encoding as described in "3D Compression
    /// Made Simple: Edgebreaker on a Corner-Table" by Rossignac at al.'01.
    /// http://www.cc.gatech.edu/~jarek/papers/CornerTableSMI.pdf
    /// </summary>
    class MeshEdgeBreakerEncoderImpl : IMeshEdgeBreakerEncoder
    {
        /// <summary>
        /// The main encoder that own's this class.
        /// </summary>
        private MeshEdgeBreakerEncoder encoder;
        /// <summary>
        /// Mesh that's being encoded.
        /// </summary>
        private DracoMesh mesh;
        /// <summary>
        /// Corner table stores the mesh face connectivity data.
        /// </summary>
        private CornerTable cornerTable;
        /// <summary>
        /// Stack used for storing corners that need to be traversed when encoding
        /// the connectivity. New corner is added for each initial face and a split
        /// symbol, and one corner is removed when the end symbol is reached.
        /// Stored as member variable to prevent frequent memory reallocations when
        /// handling meshes with lots of disjoint components. Originally, we used
        /// recursive functions to handle this behavior, but that can cause stack
        /// memory overflow when compressing huge meshes.
        /// </summary>
        private IntList cornerTraversalStack = new IntList();
        /// <summary>
        /// Array for marking visited faces.
        /// </summary>
        private List<bool> visitedFaces = new List<bool>();

        /// <summary>
        /// Attribute data for position encoding.
        /// </summary>
        private MeshAttributeIndicesEncodingData posEncodingData = new MeshAttributeIndicesEncodingData();

        /// <summary>
        /// Array storing corners in the order they were visited during the
        /// connectivity encoding (always storing the tip corner of each newly visited
        /// face).
        /// </summary>
        private IntList processedConnectivityCorners = new IntList();

        /// <summary>
        /// Array for storing visited vertex ids of all input vertices.
        /// </summary>
        private bool[] visitedVertexIds;

        /// <summary>
        /// For each traversal, this array stores the number of visited vertices.
        /// </summary>
        private IntList vertexTraversalLength = new IntList();
        /// <summary>
        /// Array for storing all topology split events encountered during the mesh
        /// traversal.
        /// </summary>
        private List<TopologySplitEventData> topologySplitEventData = new List<TopologySplitEventData>();
        /// <summary>
        /// Map between faceId and symbolId. Contains entries only for faces that
        /// were encoded with TOPOLOGYS symbol.
        /// </summary>
        private Dictionary<int, int> faceToSplitSymbolMap = new Dictionary<int, int>();

        /// <summary>
        /// Array for marking holes that has been reached during the traversal.
        /// </summary>
        private List<bool> visitedHoles = new List<bool>();

        /// <summary>
        /// Array for mapping vertices to hole ids. If a vertex is not on a hole, the
        /// stored value is -1.
        /// </summary>
        private int[] vertexHoleId;
        /// <summary>
        /// Array of hole events encountered during the traversal. There will be always
        /// exactly one hole event for each hole in the input mesh.
        /// </summary>
        private List<HoleEventData> holeEventData = new List<HoleEventData>();

        /// <summary>
        /// Id of the last encoded symbol.
        /// </summary>
        private int lastEncodedSymbolId;

        /// <summary>
        /// The number of encoded split symbols.
        /// </summary>
        private int numSplitSymbols;

        /// <summary>
        /// Struct holding data used for encoding each non-position attribute.
        /// TODO(ostava): This should be probably renamed to something better.
        /// </summary>
        class AttributeData
        {
            public int attributeIndex = -1;
            public MeshAttributeCornerTable connectivityData;
            /// <summary>
            /// Flag that can mark the connectivityData invalid. In such case the base
            /// corner table of the mesh should be used instead.
            /// </summary>
            public bool isConnectivityUsed = true;
            /// <summary>
            /// Data about attribute encoding order.
            /// </summary>
            public MeshAttributeIndicesEncodingData encodingData = new MeshAttributeIndicesEncodingData();

            public MeshTraversalMethod traversalMethod;
        };

        private AttributeData[] attributeData = null;
        private MeshTraversalMethod posTraversalMethod;

        /// <summary>
        /// Array storing mapping between attribute encoder id and attribute data id.
        /// </summary>
        private IntList attributeEncoderToDataIdMap = new IntList();

        private ITraversalEncoder traversalEncoder;

        /// <summary>
        /// Initializes data needed for encoding non-position attributes.
        /// Returns false on error.
        /// </summary>
        public void InitAttributeData()
        {

            int numAttributes = mesh.NumAttributes;
            // Ignore the position attribute. It's decoded separately.
            attributeData = new AttributeData[numAttributes - 1];
            for (int i = 0; i < attributeData.Length; i++)
            {
                attributeData[i] = new AttributeData();
            }
            if (numAttributes == 1)
                return;
            int dataIndex = 0;
            for (int i = 0; i < numAttributes; ++i)
            {
                int attIndex = i;
                if (mesh.Attribute(attIndex).AttributeType ==
                    AttributeType.Position)
                    continue;
                PointAttribute att = mesh.Attribute(attIndex);
                attributeData[dataIndex].attributeIndex = attIndex;
                attributeData[dataIndex]
                    .encodingData.encodedAttributeValueIndexToCornerMap.Clear();

                attributeData[dataIndex]
                    .encodingData.encodedAttributeValueIndexToCornerMap.Capacity = cornerTable.NumCorners;
                DracoUtils.Fill(attributeData[dataIndex].encodingData.vertexToEncodedAttributeValueIndexMap =
                    new int[cornerTable.NumCorners], -1);

                attributeData[dataIndex].encodingData.numValues = 0;
                attributeData[dataIndex].connectivityData = new MeshAttributeCornerTable(mesh, cornerTable, att);
                ++dataIndex;
            }
        }

        private void Assign(IntList list, int size, int val)
        {
            //Make list to the specified size and all content to specified value
            //remove unnecessary values
            while (list.Count > size)
                list.RemoveAt(list.Count - 1);
            //replace old values
            for (int i = 0; i < list.Count; i++)
                list[i] = val;
            //extend to specified size
            while (list.Count < size)
                list.Add(val);
        }
        private void Assign<T>(List<T> list, int size, T val)
        {
            //Make list to the specified size and all content to specified value
            //remove unnecessary values
            while (list.Count > size)
                list.RemoveAt(list.Count - 1);
            //replace old values
            for (int i = 0; i < list.Count; i++)
                list[i] = val;
            //extend to specified size
            while (list.Count < size)
                list.Add(val);
        }

        /// <summary>
        /// Finds the configuration of the initial face that starts the traversal.
        /// Configurations are determined by location of holes around the init face
        /// and they are described in meshEdgebreakerShared.h.
        /// Returns true if the face configuration is interior and false if it is
        /// exterior.
        /// </summary>
        public bool FindInitFaceConfiguration(int faceId, out int outCorner)
        {

            int cornerIndex = 3 * faceId;
            for (int i = 0; i < 3; ++i)
            {
                if (cornerTable.Opposite(cornerIndex) == -1)
                {
                    // If there is a boundary edge, the configuration is exterior and return
                    // the int opposite to the first boundary edge.
                    outCorner = cornerIndex;
                    return false;
                }
                if (vertexHoleId[cornerTable.Vertex(cornerIndex)] != -1)
                {
                    // Boundary vertex found. Find the first boundary edge attached to the
                    // point and return the corner opposite to it.
                    int rightCorner = cornerIndex;
                    while (rightCorner >= 0)
                    {
                        cornerIndex = rightCorner;
                        rightCorner = cornerTable.SwingRight(rightCorner);
                    }
                    // |cornerIndex| now lies on a boundary edge and its previous corner is
                    // guaranteed to be the opposite corner of the boundary edge.
                    outCorner = cornerTable.Previous(cornerIndex);
                    return false;
                }
                cornerIndex = cornerTable.Next(cornerIndex);
            }
            // Else we have an interior configuration. Return the first corner id.
            outCorner = cornerIndex;
            return true;
        }

        /// <summary>
        /// Encodes the connectivity between vertices.
        /// </summary>
        public void EncodeConnectivityFromCorner(int cornerId)
        {

            cornerTraversalStack.Clear();
            cornerTraversalStack.Add(cornerId);
            int numFaces = mesh.NumFaces;
            while (cornerTraversalStack.Count > 0)
            {
                // Currently processed corner.
                cornerId = cornerTraversalStack[cornerTraversalStack.Count - 1];
                // Make sure the face hasn't been visited yet.
                if (cornerId < 0 ||
                    visitedFaces[cornerTable.Face(cornerId)])
                {
                    // This face has been already traversed.
                    cornerTraversalStack.RemoveAt(cornerTraversalStack.Count - 1);
                    continue;
                }
                int numVisitedFaces = 0;
                while (numVisitedFaces < numFaces)
                {
                    // Mark the current face as visited.
                    ++numVisitedFaces;
                    ++lastEncodedSymbolId;

                    int faceId = cornerTable.Face(cornerId);
                    visitedFaces[faceId] = true;
                    processedConnectivityCorners.Add(cornerId);
                    traversalEncoder.NewCornerReached(cornerId);
                    int vertId = cornerTable.Vertex(cornerId);
                    bool onBoundary = (vertexHoleId[vertId] != -1);
                    if (!IsVertexVisited(vertId))
                    {
                        // A new unvisited vertex has been reached. We need to store its
                        // position difference using next,prev, and opposite vertices.
                        visitedVertexIds[vertId] = true;
                        if (!onBoundary)
                        {
                            // If the vertex is on boundary it must correspond to an unvisited
                            // hole and it will be encoded with TOPOLOGYS symbol later).
                            traversalEncoder.EncodeSymbol(EdgeBreakerTopologyBitPattern.C);
                            // Move to the right triangle.
                            cornerId = GetRightCorner(cornerId);
                            continue;
                        }
                    }
                    // The current vertex has been already visited or it was on a boundary.
                    // We need to determine whether we can visit any of it's neighboring
                    // faces.
                    int rightCornerId = GetRightCorner(cornerId);
                    int leftCornerId = GetLeftCorner(cornerId);
                    int rightFaceId = cornerTable.Face(rightCornerId);
                    int leftFaceId = cornerTable.Face(leftCornerId);
                    if (IsRightFaceVisited(cornerId))
                    {
                        // Right face has been already visited.
                        // Check whether there is a topology split event.
                        if (rightFaceId != -1)
                            CheckAndStoreTopologySplitEvent(lastEncodedSymbolId,
                                faceId, EdgeFaceName.RightFaceEdge,
                                rightFaceId);
                        if (IsLeftFaceVisited(cornerId))
                        {
                            // Both neighboring faces are visited. End reached.
                            // Check whether there is a topology split event on the left face.
                            if (leftFaceId != -1)
                                CheckAndStoreTopologySplitEvent(lastEncodedSymbolId,
                                    faceId, EdgeFaceName.LeftFaceEdge,
                                    leftFaceId);
                            traversalEncoder.EncodeSymbol(EdgeBreakerTopologyBitPattern.E);
                            cornerTraversalStack.RemoveAt(cornerTraversalStack.Count - 1);
                            break; // Break from the while (numVisitedFaces < numFaces) loop.
                        }
                        else
                        {
                            traversalEncoder.EncodeSymbol(EdgeBreakerTopologyBitPattern.R);
                            // Go to the left face.
                            cornerId = leftCornerId;
                        }
                    }
                    else
                    {
                        // Right face was not visited.
                        if (IsLeftFaceVisited(cornerId))
                        {
                            // Check whether there is a topology split event on the left face.
                            if (leftFaceId != -1)
                                CheckAndStoreTopologySplitEvent(lastEncodedSymbolId,
                                    faceId, EdgeFaceName.LeftFaceEdge,
                                    leftFaceId);
                            traversalEncoder.EncodeSymbol(EdgeBreakerTopologyBitPattern.L);
                            // Left face visited, go to the right one.
                            cornerId = rightCornerId;
                        }
                        else
                        {
                            traversalEncoder.EncodeSymbol(EdgeBreakerTopologyBitPattern.S);
                            ++numSplitSymbols;
                            // Both neighboring faces are unvisited, we need to visit both of
                            // them.
                            if (onBoundary)
                            {
                                // The tip vertex is on a hole boundary. If the hole hasn't been
                                // visited yet we need to encode it.
                                int holeId = vertexHoleId[vertId];
                                if (!visitedHoles[holeId])
                                {
                                    EncodeHole(cornerId, false);
                                    holeEventData.Add(
                                        new HoleEventData(lastEncodedSymbolId));
                                }
                            }
                            faceToSplitSymbolMap[faceId] = lastEncodedSymbolId;
                            // Split the traversal.
                            // First make the top of the current corner stack point to the left
                            // face (this one will be processed second).
                            cornerTraversalStack[cornerTraversalStack.Count - 1] = leftCornerId;
                            // Add a new corner to the top of the stack (right face needs to
                            // be traversed first).
                            cornerTraversalStack.Add(rightCornerId);
                            // Break from the while (numVisitedFaces < numFaces) loop.
                            break;
                        }
                    }
                }
            }
            // All corners have been processed.
        }

        /// <summary>
        /// Encodes all vertices of a hole starting at startCornerId.
        /// The vertex associated with the first corner is encoded only if
        /// |encodeFirstVertex| is true.
        /// Returns the number of encoded hole vertices.
        /// </summary>
        public int EncodeHole(int startCornerId, bool encodeFirstVertex)
        {

            // We know that the start corner lies on a hole but we first need to find the
            // boundary edge going from that vertex. It is the first edge in CW
            // direction.
            int cornerId = startCornerId;
            cornerId = cornerTable.Previous(cornerId);
            while (cornerTable.Opposite(cornerId) != -1)
            {
                cornerId = cornerTable.Opposite(cornerId);
                cornerId = cornerTable.Next(cornerId);
            }
            int startVertexId = cornerTable.Vertex(startCornerId);

            int numEncodedHoleVerts = 0;
            if (encodeFirstVertex)
            {
                visitedVertexIds[startVertexId] = true;
                ++numEncodedHoleVerts;
            }

            // cornerId is now opposite to the boundary edge.
            // Mark the hole as visited.
            visitedHoles[vertexHoleId[startVertexId]] = true;
            // Get the start vertex of the edge and use it as a reference.
            int startVertId =
                cornerTable.Vertex(cornerTable.Next(cornerId));
            // Get the end vertex of the edge.
            int actVertexId =
                cornerTable.Vertex(cornerTable.Previous(cornerId));
            while (actVertexId != startVertexId)
            {
                // Encode the end vertex of the boundary edge.

                startVertId = actVertexId;

                // Mark the vertex as visited.
                visitedVertexIds[actVertexId] = true;
                ++numEncodedHoleVerts;
                cornerId = cornerTable.Next(cornerId);
                // Look for the next attached open boundary edge.
                while (cornerTable.Opposite(cornerId) != -1)
                {
                    cornerId = cornerTable.Opposite(cornerId);
                    cornerId = cornerTable.Next(cornerId);
                }
                actVertexId = cornerTable.Vertex(cornerTable.Previous(cornerId));
            }
            return numEncodedHoleVerts;
        }

        public int GetRightCorner(int cornerId)
        {

            int nextCornerId = cornerTable.Next(cornerId);
            return cornerTable.Opposite(nextCornerId);
        }

        public int GetLeftCorner(int cornerId)
        {

            int prevCornerId = cornerTable.Previous(cornerId);
            return cornerTable.Opposite(prevCornerId);
        }

        public bool IsRightFaceVisited(int cornerId)
        {

            int nextCornerId = cornerTable.Next(cornerId);
            int oppCornerId = cornerTable.Opposite(nextCornerId);
            if (oppCornerId != -1)
                return visitedFaces[cornerTable.Face(oppCornerId)];
            // Else we are on a boundary.
            return true;
        }

        public bool IsLeftFaceVisited(int cornerId)
        {

            int prevCornerId = cornerTable.Previous(cornerId);
            int oppCornerId = cornerTable.Opposite(prevCornerId);
            if (oppCornerId != -1)
                return visitedFaces[cornerTable.Face(oppCornerId)];
            // Else we are on a boundary.
            return true;
        }

        public bool IsVertexVisited(int vertId)
        {
            return visitedVertexIds[vertId];
        }

        /// <summary>
        /// Finds and stores data about all holes in the input mesh.
        /// </summary>
        public bool FindHoles()
        {

            // TODO(ostava): Add more error checking for invalid geometry data.
            int numCorners = cornerTable.NumCorners;
            // Go over all corners and detect non-visited open boundaries
            for (int i = 0; i < numCorners; ++i)
            {
                if (cornerTable.IsDegenerated(cornerTable.Face(i)))
                    continue; // Don't process corners assigned to degenerated faces.
                if (cornerTable.Opposite(i) == -1)
                {
                    // No opposite corner means no opposite face, so the opposite edge
                    // of the corner is an open boundary.
                    // Check whether we have already traversed the boundary.
                    int boundaryVertId = cornerTable.Vertex(cornerTable.Next(i));
                    if (vertexHoleId[boundaryVertId] != -1)
                    {
                        // The start vertex of the boundary edge is already assigned to an
                        // open boundary. No need to traverse it again.
                        continue;
                    }
                    // Else we found a new open boundary and we are going to traverse along it
                    // and mark all visited vertices.
                    int boundaryId = visitedHoles.Count;
                    visitedHoles.Add(false);

                    int cornerId = i;
                    while (vertexHoleId[boundaryVertId] == -1)
                    {
                        // Mark the first vertex on the open boundary.
                        vertexHoleId[boundaryVertId] = boundaryId;
                        cornerId = cornerTable.Next(cornerId);
                        // Look for the next attached open boundary edge.
                        while (cornerTable.Opposite(cornerId) != -1)
                        {
                            cornerId = cornerTable.Opposite(cornerId);
                            cornerId = cornerTable.Next(cornerId);
                        }
                        // Id of the next vertex in the vertex on the hole.
                        boundaryVertId =
                            cornerTable.Vertex(cornerTable.Next(cornerId));
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// For faces encoded with symbol TOPOLOGYS (split), this method returns
        /// the encoded symbol id or -1 if the face wasn't encoded by a split symbol.
        /// </summary>
        public int GetSplitSymbolIdOnFace(int faceId)
        {
            int ret;
            if (faceToSplitSymbolMap.TryGetValue(faceId, out ret))
                return ret;
            return -1;
        }

        /// <summary>
        /// Checks whether there is a topology split event on a neighboring face and
        /// stores the event data if necessary. For more info about topology split
        /// events, see description of TopologySplitEventData in
        /// meshEdgebreakerShared.h.
        /// </summary>
        public void CheckAndStoreTopologySplitEvent(int srcSymbolId, int srcFaceId,
            EdgeFaceName srcEdge,
            int neighborFaceId)
        {

            int symbolId = GetSplitSymbolIdOnFace(neighborFaceId);
            if (symbolId == -1)
                return; // Not a split symbol, no topology split event could happen.
            TopologySplitEventData eventData = new TopologySplitEventData();

            eventData.splitSymbolId = symbolId;
            // It's always the left edge for true split symbols (as the right edge is
            // traversed first).
            eventData.splitEdge = (byte) EdgeFaceName.LeftFaceEdge;

            eventData.sourceSymbolId = srcSymbolId;
            eventData.sourceEdge = (byte) srcEdge;
            topologySplitEventData.Add(eventData);
        }

        /// <summary>
        /// Encodes connectivity of all attributes on a newly traversed face.
        /// </summary>
        public void EncodeAttributeConnectivitiesOnFace(int corner)
        {

            // Three corners of the face.
            int[] corners =
            {
                corner, cornerTable.Next(corner),
                cornerTable.Previous(corner)
            };

            int src_face_id = cornerTable.Face(corner);
            visitedFaces[src_face_id] = true;

            for (int c = 0; c < 3; ++c)
            {
                int oppCorner = cornerTable.Opposite(corners[c]);
                if (oppCorner < 0)
                    continue; // Don't encode attribute seams on boundary edges.

                int opp_face_id = cornerTable.Face(oppCorner);
                if (visitedFaces[opp_face_id])
                    continue;

                for (int i = 0; i < attributeData.Length; ++i)
                {
                    if (attributeData[i].connectivityData.IsCornerOppositeToSeamEdge(
                        corners[c]))
                    {
                        traversalEncoder.EncodeAttributeSeam(i, true);
                    }
                    else
                    {
                        traversalEncoder.EncodeAttributeSeam(i, false);
                    }
                }
            }
        }

        /// <summary>
        /// This function is used to to assign correct encoding order of attributes
        /// to unprocessed corners. The encoding order is equal to the order in which
        /// the attributes are going to be processed by the decoder and it is necessary
        /// for proper prediction of attribute values.
        ///public bool AssignPositionEncodingOrderToAllCorners();
        /// This function is used to generate encoding order for all non-position
        /// attributes.
        /// Returns false when one or more attributes failed to be processed.
        ///public bool GenerateEncodingOrderForAttributes();
        /// </summary>

        public MeshEdgeBreakerEncoderImpl(ITraversalEncoder encoder)
        {
            lastEncodedSymbolId = -1;
            this.traversalEncoder = encoder;
        }

        public void Init(MeshEdgeBreakerEncoder encoder)
        {

            this.encoder = encoder;
            mesh = encoder.Mesh;
            attributeEncoderToDataIdMap.Clear();
        }

        public MeshAttributeCornerTable GetAttributeCornerTable(int attId)
        {

            for (int i = 0; i < attributeData.Length; ++i)
            {
                if (attributeData[i].attributeIndex == attId)
                {
                    if (attributeData[i].isConnectivityUsed)
                        return attributeData[i].connectivityData;
                    return null;
                }
            }
            return null;
        }

        public MeshAttributeIndicesEncodingData GetAttributeEncodingData(int attId)
        {
            for (int i = 0; i < attributeData.Length; ++i)
            {
                if (attributeData[i].attributeIndex == attId)
                    return attributeData[i].encodingData;
            }
            return posEncodingData;
        }

        public void GenerateAttributesEncoder(int attId)
        {

            // For now, create one encoder for each attribute. Ideally we can share
            // the same encoder for attributes with the same connectivity (this is
            // especially true for per-vertex attributes).
            var elementType = Encoder.Mesh.GetAttributeElementType(attId);
            PointAttribute att = Encoder.PointCloud.Attribute(attId);
            int attDataId = -1;
            for (int i = 0; i < attributeData.Length; ++i)
            {
                if (attributeData[i].attributeIndex == attId)
                {
                    attDataId = i;
                    break;
                }
            }
            MeshTraversalMethod traversalMethod = MeshTraversalMethod.DepthFirst;
            PointsSequencer sequencer = null;
            if (att.AttributeType == AttributeType.Position || elementType == MeshAttributeElementType.Vertex ||
                (elementType == MeshAttributeElementType.Corner &&
                 attributeData[attDataId].connectivityData.NoInteriorSeams))
            {
                // Per-vertex attribute reached, use the basic corner table to traverse the
                // mesh.
                // Traverser that is used to generate the encoding order of each attribute.

                MeshAttributeIndicesEncodingData encodingData;
                if (att.AttributeType == AttributeType.Position)
                {
                    encodingData = posEncodingData;
                }
                else
                {
                    encodingData = attributeData[attDataId].encodingData;
                    attributeData[attDataId].isConnectivityUsed = false;
                }

                int speed = encoder.Options.GetSpeed();
                if (speed == 0 && att.AttributeType == AttributeType.Position)
                {
                    traversalMethod = MeshTraversalMethod.PredictionDegree;
                    if (mesh.NumAttributes > 1)
                    {
                        // Make sure we don't use the prediction degree traversal when we encode
                        // multiple attributes using the same connectivity.
                        // TODO(ostava): We should investigate this and see if the prediction
                        // degree can be actually used efficiently for non-position attributes.
                        traversalMethod = MeshTraversalMethod.DepthFirst;
                    }
                }

                // Defining sequencer via a traversal scheme.
                MeshTraversalSequencer<CornerTable> traversalSequencer = new MeshTraversalSequencer<CornerTable>(
                    mesh, encodingData);
                var attObserver = new MeshAttributeIndicesEncodingObserver<CornerTable>(cornerTable, mesh,
                    traversalSequencer, encodingData);
                ICornerTableTraverser<CornerTable> attTraverser = null;
                if (traversalMethod == MeshTraversalMethod.PredictionDegree) {
                  //typedef MeshAttributeIndicesEncodingObserver<CornerTable> AttObserver;
                  //typedef MaxPredictionDegreeTraverser<CornerTable, AttObserver> AttTraverser;
                    attTraverser = new EdgeBreakerTraverser<CornerTable>(cornerTable, attObserver);
                } else if (traversalMethod == MeshTraversalMethod.DepthFirst) {
                    //typedef MeshAttributeIndicesEncodingObserver<CornerTable> AttObserver;
                    //typedef DepthFirstTraverser<CornerTable, AttObserver> AttTraverser;
                    //sequencer = CreateVertexTraversalSequencer<AttTraverser>(encoding_data);
                    attTraverser = new DepthFirstTraverser(cornerTable, attObserver);
                }

                traversalSequencer.SetCornerOrder(processedConnectivityCorners);
                traversalSequencer.SetTraverser(attTraverser);
                sequencer = traversalSequencer;
            }
            else
            {
                // Else use a general per-corner encoder.
                // Traverser that is used to generate the encoding order of each attribute.

                var traversalSequencer = new MeshTraversalSequencer<MeshAttributeCornerTable>(mesh, attributeData[attDataId].encodingData) ;

                var attObserver =
                    new MeshAttributeIndicesEncodingObserver<MeshAttributeCornerTable>(
                        attributeData[attDataId].connectivityData,
                        mesh, traversalSequencer,
                        attributeData[attDataId].encodingData);

                var attTraverser = new EdgeBreakerTraverser<MeshAttributeCornerTable>(attributeData[attDataId].connectivityData, attObserver);

                traversalSequencer.SetCornerOrder(processedConnectivityCorners);
                traversalSequencer.SetTraverser(attTraverser);
                sequencer = traversalSequencer;
            }

            if (sequencer == null)
                throw DracoUtils.Failed();

            if (attDataId == -1)
                posTraversalMethod = traversalMethod;
            else
                attributeData[attDataId].traversalMethod = traversalMethod;

            SequentialAttributeEncodersController attController = new SequentialAttributeEncodersController(sequencer,
                attId);

            // Update the mapping between the encoder id and the attribute data id.
            // This will be used by the decoder to select the approperiate attribute
            // decoder and the correct connectivity.
            attributeEncoderToDataIdMap.Add(attDataId);
            Encoder.AddAttributesEncoder(attController);
        }

        public void EncodeAttributesEncoderIdentifier(int attEncoderId)
        {

            int attDataId = attributeEncoderToDataIdMap[attEncoderId];
            MeshTraversalMethod traversalMethod;
            encoder.Buffer.Encode((byte)attDataId);

            // Also encode the type of the encoder that we used.
            var elementType = MeshAttributeElementType.Vertex; 
            if (attDataId >= 0)
            {
                int attId = attributeData[attDataId].attributeIndex;
                elementType = Encoder.Mesh.GetAttributeElementType(attId);
                traversalMethod = attributeData[attDataId].traversalMethod;
            }
            else
            {
                traversalMethod = posTraversalMethod;
            }
            if (elementType == MeshAttributeElementType.Vertex ||
                (elementType == MeshAttributeElementType.Corner &&
                 attributeData[attDataId].connectivityData.NoInteriorSeams))
            {
                // Per-vertex encoder.
                encoder.Buffer.Encode((byte) (MeshAttributeElementType.Vertex));
            }
            else
            {
                // Per-corner encoder.
                encoder.Buffer.Encode((byte) (MeshAttributeElementType.Corner));
            }
            encoder.Buffer.Encode((byte) traversalMethod);
        }

        private CornerTable CreateCornerTableFromPositionAttribute(DracoMesh mesh)
        {
            PointAttribute att = mesh.GetNamedAttribute(AttributeType.Position);
            if (att == null)
                return null;
            int[,] faces = new int[mesh.NumFaces, 3];
            Span<int> face = stackalloc int[3];
            for (int i = 0; i < mesh.NumFaces; ++i)
            {
                mesh.ReadFace(i, face);
                for (int j = 0; j < 3; ++j)
                {
                    // Map general vertex indices to position indices.
                    faces[i, j] = att.MappedIndex(face[j]);
                }
            }

            // Build the corner table.
            var ret = new CornerTable();
            ret.Initialize(faces);
            return ret;
        }

        private CornerTable CreateCornerTableFromAllAttributes(DracoMesh mesh)
        {
            int[,] faces = new int[mesh.NumFaces, 3];
            Span<int> face = stackalloc int[3];
            for (int i = 0; i < mesh.NumFaces; ++i)
            {
                mesh.ReadFace(i, face);
                // Each face is identified by point indices that automatically split the
                // mesh along attribute seams.
                for (int j = 0; j < 3; ++j)
                {
                    faces[i, j] = face[j];
                }
            }
            // Build the corner table.
            var ret = new CornerTable();
            ret.Initialize(faces);
            return ret;
        }

        public void EncodeConnectivity()
        {

            // To encode the mesh, we need face connectivity data stored in a corner
            // table. To compute the connectivity we must use indices associated with
            // POSITION attribute, because they define which edges can be connected
            // together.

            if (encoder.Options.SplitMeshOnSeams)
                cornerTable = CreateCornerTableFromAllAttributes(mesh);
            else
                cornerTable = CreateCornerTableFromPositionAttribute(mesh);
            if (cornerTable == null || cornerTable.NumFaces == cornerTable.NumDegeneratedFaces)
            {
                // Failed to ruct the corner table.
                throw DracoUtils.Failed();
            }

            traversalEncoder.Init(this);


            // Also encode the total number of vertices that is going to be encoded.
            // This can be different from the mesh.numPoints() + numNewVertices,
            // because some of the vertices of the input mesh can be ignored (e.g.
            // vertices on degenerated faces or isolated vertices not attached to any
            // face).
            int numVerticesToBeEncoded =
                cornerTable.NumVertices - cornerTable.NumIsolatedVertices;
            Encoding.EncodeVarint((uint)numVerticesToBeEncoded, Encoder.Buffer);

            int numFaces =
                cornerTable.NumFaces - cornerTable.NumDegeneratedFaces;
            Encoding.EncodeVarint((uint)numFaces, Encoder.Buffer);

            // Reset encoder data that may have been initialized in previous runs.
            Assign(visitedFaces, mesh.NumFaces, false);
            DracoUtils.Fill(posEncodingData.vertexToEncodedAttributeValueIndexMap = new int[cornerTable.NumVertices], -1);
            posEncodingData.encodedAttributeValueIndexToCornerMap.Clear();
            posEncodingData.encodedAttributeValueIndexToCornerMap.Capacity = cornerTable.NumFaces * 3;
            //Assign(visitedVertexIds, cornerTable.NumVertices, false);
            visitedVertexIds = new bool[cornerTable.NumVertices];
            vertexTraversalLength.Clear();
            lastEncodedSymbolId = -1;
            numSplitSymbols = 0;
            topologySplitEventData.Clear();
            faceToSplitSymbolMap.Clear();
            visitedHoles.Clear();
            //Assign(vertexHoleId, cornerTable.NumVertices, -1);
            DracoUtils.Fill(vertexHoleId = new int[cornerTable.NumVertices], -1);
            holeEventData.Clear();
            processedConnectivityCorners.Clear();
            processedConnectivityCorners.Capacity = cornerTable.NumFaces;
            posEncodingData.numValues = 0;

            if (!FindHoles())
                throw DracoUtils.Failed();

            InitAttributeData();

            byte numAttributeData = (byte) attributeData.Length;
            encoder.Buffer.Encode(numAttributeData);

            int numCorners = cornerTable.NumCorners;

            traversalEncoder.Start();

            IntList initFaceConnectivityCorners = new IntList();
            // Traverse the surface starting from each unvisited corner.
            for (int cId = 0; cId < numCorners; ++cId)
            {
                int cornerIndex = cId;
                int faceId = cornerTable.Face(cornerIndex);
                if (visitedFaces[faceId])
                    continue; // Face has been already processed.
                if (cornerTable.IsDegenerated(faceId))
                    continue; // Ignore degenerated faces.

                int startCorner;
                bool interiorConfig =
                    FindInitFaceConfiguration(faceId, out startCorner);
                traversalEncoder.EncodeStartFaceConfiguration(interiorConfig);

                if (interiorConfig)
                {
                    // Select the correct vertex on the face as the root.
                    cornerIndex = startCorner;
                    int vertId = cornerTable.Vertex(cornerIndex);
                    // Mark all vertices of a given face as visited.
                    int nextVertId =
                        cornerTable.Vertex(cornerTable.Next(cornerIndex));
                    int prevVertId =
                        cornerTable.Vertex(cornerTable.Previous(cornerIndex));

                    visitedVertexIds[vertId] = true;
                    visitedVertexIds[nextVertId] = true;
                    visitedVertexIds[prevVertId] = true;
                    // New traversal started. Initiate it's length with the first vertex.
                    vertexTraversalLength.Add(1);

                    // Mark the face as visited.
                    visitedFaces[faceId] = true;
                    // Start compressing from the opposite face of the "next" corner. This way
                    // the first encoded corner corresponds to the tip corner of the regular
                    // edgebreaker traversal (essentially the initial face can be then viewed
                    // as a TOPOLOGYC face).
                    initFaceConnectivityCorners.Add(
                        cornerTable.Next(cornerIndex));
                    int oppId =
                        cornerTable.Opposite(cornerTable.Next(cornerIndex));
                    int oppFaceId = cornerTable.Face(oppId);
                    if (oppFaceId != -1 &&
                        !visitedFaces[oppFaceId])
                    {
                        EncodeConnectivityFromCorner(oppId);
                    }
                }
                else
                {
                    // Bounary configuration. We start on a boundary rather than on a face.
                    // First encode the hole that's opposite to the startCorner.
                    EncodeHole(cornerTable.Next(startCorner), true);
                    // Start processing the face opposite to the boundary edge (the face
                    // containing the startCorner).
                    EncodeConnectivityFromCorner(startCorner);
                }
            }
            // Reverse the order of connectivity corners to match the order in which
            // they are going to be decoded.
            processedConnectivityCorners.Reverse();
            // Append the init face connectivity corners (which are processed in order by
            // the decoder after the regular corners.
            processedConnectivityCorners.AddRange(initFaceConnectivityCorners);
            // Emcode connectivity for all non-position attributes.
            if (attributeData.Length > 0)
            {
                // Use the same order of corner that will be used by the decoder.
                for (int i = 0; i < visitedFaces.Count; i++)
                    visitedFaces[i] = false;
                for (int i = 0; i < processedConnectivityCorners.Count; i++)
                {
                    var ci = processedConnectivityCorners[i];
                    EncodeAttributeConnectivitiesOnFace(ci);
                }
            }
            traversalEncoder.Done();

            // Encode the number of symbols.
            Encoding.EncodeVarint((uint)traversalEncoder.NumEncodedSymbols, encoder.Buffer);

            // Encode the number of split symbols.
            Encoding.EncodeVarint((uint)numSplitSymbols, encoder.Buffer);

            // Append the traversal buffer.

            EncodeSplitData();

            encoder.Buffer.Encode(traversalEncoder.Buffer.Data, traversalEncoder.Buffer.Bytes);
        }

        void EncodeSplitData()
        {
            int numEvents = topologySplitEventData.Count;
            Encoding.EncodeVarint((uint)numEvents, encoder.Buffer);
            if (numEvents > 0)
            {
                // Encode split symbols using delta and varint coding. Split edges are
                // encoded using direct bit coding.
                int last_source_symbol_id = 0; // Used for delta coding.
                for (int i = 0; i < numEvents; ++i)
                {
                    TopologySplitEventData event_data = topologySplitEventData[i];
                    // Encode source symbol id as delta from the previous source symbol id.
                    // Source symbol ids are always stored in increasing order so the delta is
                    // going to be positive.
                    Encoding.EncodeVarint((uint)(event_data.sourceSymbolId - last_source_symbol_id), encoder.Buffer);
                    // Encode split symbol id as delta from the current source symbol id.
                    // Split symbol id is always smaller than source symbol id so the below
                    // delta is going to be positive.
                    Encoding.EncodeVarint((uint)(event_data.sourceSymbolId - event_data.splitSymbolId), encoder.Buffer);
                    last_source_symbol_id = event_data.sourceSymbolId;
                }

                encoder.Buffer.StartBitEncoding(numEvents, false);
                for (int i = 0; i < numEvents; ++i)
                {
                    TopologySplitEventData event_data = topologySplitEventData[i];
                    encoder.Buffer.EncodeLeastSignificantBits32(1, event_data.sourceEdge);
                }

                encoder.Buffer.EndBitEncoding();
            }

        }

        public CornerTable CornerTable
        {
            get { return cornerTable; }
        }

        public MeshEdgeBreakerEncoder Encoder
        {
            get { return encoder; }
        }

    }
}
