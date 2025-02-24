using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;
using FileFormat.Drako.Compression;

namespace FileFormat.Drako.Decoder
{

    interface ITraversalDecoder
    {


        /// <summary>
        /// Returns true if there is an attribute seam for the next processed pair
        /// of visited faces.
        /// |attribute| is used to mark the id of the non-position attribute (in range
        /// of &lt;0, numAttributes - 1&gt;).
        /// </summary>
        bool DecodeAttributeSeam(int attribute);
        void Init(IMeshEdgeBreakerDecoderImpl decoder);
        /// <summary>
        /// Used to tell the decoder what is the number of expected decoded vertices.
        /// Ignored by default.
        /// </summary>
        void SetNumEncodedVertices(int numVertices);

        /// <summary>
        /// Set the number of non-position attribute data for which we need to decode
        /// the connectivity.
        /// </summary>
        void SetNumAttributeData(int numData);

        /// <summary>
        /// Called before the traversal decoding is started.
        /// Returns a buffer decoder that points to data that was encoded after the
        /// traversal.
        /// </summary>
        DecoderBuffer Start();
        void Done();
        /// <summary>
        /// Returns the next edgebreaker symbol that was reached during the traversal.
        /// </summary>
        EdgeBreakerTopologyBitPattern DecodeSymbol();
        /// <summary>
        /// Called whenever |source| vertex is about to be merged into the |dest|
        /// vertex.
        /// </summary>
        void MergeVertices(int dest, int source);
        /// <summary>
        /// Called whenever a new active corner is set in the decoder.
        /// </summary>
        void NewActiveCornerReached(int corner);
        /// <summary>
        /// Returns the configuration of a new initial face.
        /// </summary>
        bool DecodeStartFaceConfiguration();
    }


    /// <summary>
    /// Implementation of the edgebreaker decoder that decodes data encoded with the
    /// MeshEdgeBreakerEncoderImpl class. The implementation of the decoder is based
    /// on the algorithm presented in Isenburg et al'02 "Spirale Reversi: Reverse
    /// decoding of the Edgebreaker encoding". Note that the encoding is still based
    /// on the standard edgebreaker method as presented in "3D Compression
    /// Made Simple: Edgebreaker on a Corner-Table" by Rossignac at al.'01.
    /// http://www.cc.gatech.edu/~jarek/papers/CornerTableSMI.pdf. One difference is
    /// caused by the properties of the spirale reversi algorithm that decodes the
    /// symbols from the last one to the first one. To make the decoding more
    /// efficient, we encode all symbols in the reverse order, therefore the decoder
    /// can process them one by one.
    /// The main advantage of the spirale reversi method is that the partially
    /// decoded mesh has valid connectivity data at any time during the decoding
    /// process (valid with respect to the decoded portion of the mesh). The standard
    /// Edgebreaker decoder used two passes (forward decoding + zipping) which not
    /// only prevented us from having a valid connectivity but it was also slower.
    /// The main benefit of having the valid connectivity is that we can use the
    /// known connectivity to predict encoded symbols that can improve the
    /// compression rate.
    /// </summary>
    class MeshEdgeBreakerDecoderImpl : IMeshEdgeBreakerDecoderImpl
    {
        private MeshEdgeBreakerDecoder decoder;

        private CornerTable cornerTable;

        /// <summary>
        /// Stack used for storing corners that need to be traversed when decoding
        /// mesh vertices. New corner is added for each initial face and a split
        /// symbol, and one corner is removed when the end symbol is reached.
        /// Stored as member variable to prevent frequent memory reallocations when
        /// handling meshes with lots of disjoint components.  Originally, we used
        /// recursive functions to handle this behavior, but that can cause stack
        /// memory overflow when compressing huge meshes.
        /// </summary>
        private IntList cornerTraversalStack = new IntList();

        /// <summary>
        /// Array stores the number of visited visited for each mesh traversal.
        /// </summary>
        private IntList vertexTraversalLength = new IntList();

        /// <summary>
        /// List of decoded topology split events.
        /// </summary>
        private List<TopologySplitEventData> topologySplitData = new List<TopologySplitEventData>();

        /// <summary>
        /// List of decoded hole events.
        /// </summary>
        private List<HoleEventData> holeEventData = new List<HoleEventData>();

#pragma warning disable 0414
        /// <summary>
        /// The number of processed hole events.
        /// </summary>
        private int numProcessedHoleEvents;

        /// <summary>
        /// Configuration of the initial face for each mesh component.
        /// </summary>
        private List<bool> initFaceConfigurations = new List<bool>();

        /// <summary>
        /// Initial corner for each traversal.
        /// </summary>
        private IntList initCorners = new IntList();

        /// <summary>
        /// Mapping between vertex ids assigned during connectivity decoding and vertex
        /// ids that were used during encoding.
        /// </summary>
        private IntList vertexIdMap = new IntList();

        /// <summary>
        /// Id of the last processed input symbol.
        /// </summary>
        private int lastSymbolId = -1;

        /// <summary>
        /// Id of the last decoded vertex.
        /// </summary>
        private int lastVertId = -1;

        /// <summary>
        /// Id of the last decoded face.
        /// </summary>
        private int lastFaceId = -1;

#pragma warning restore 0414

        /// <summary>
        /// Array for marking visited faces.
        /// </summary>
        private List<bool> visitedFaces = new List<bool>();
        /// <summary>
        /// Array for marking visited vertices.
        /// </summary>
        private List<bool> visitedVerts = new List<bool>();

        /// <summary>
        /// Array for marking vertices on open boundaries.
        /// </summary>
        private bool[] isVertHole;

        /// <summary>
        /// The number of new vertices added by the encoder (because of non-manifold
        /// vertices on the input mesh).
        /// If there are no non-manifold edges/vertices on the input mesh, this should
        /// be 0.
        /// </summary>
        private int numNewVertices;
        /// <summary>
        /// For every newly added vertex, this array stores it's mapping to the
        /// parent vertex id of the encoded mesh.
        /// </summary>
        private Dictionary<int, int> newToParentVertexMap = new Dictionary<int, int>();
        /// <summary>
        /// The number of vertices that were encoded (can be different from the number
        /// of vertices of the input mesh).
        /// </summary>
        private int numEncodedVertices;

        /// <summary>
        /// Array for storing the encoded corner ids in the order their associated
        /// vertices were decoded.
        /// </summary>
        private IntList processedCornerIds = new IntList();

        /// <summary>
        /// Array storing corners in the order they were visited during the
        /// connectivity decoding (always storing the tip corner of each newly visited
        /// face).
        /// </summary>
        private IntList processedConnectivityCorners = new IntList();

        private MeshAttributeIndicesEncodingData posEncodingData = new MeshAttributeIndicesEncodingData();


        /// <summary>
        /// Id of an attributes decoder that uses |pos_encoding_data_|.
        /// </summary>
        private int posDataDecoderId = -1;

        /// <summary>
        /// Data for non-position attributes used by the decoder.
        /// </summary>
        class AttributeData
        {
            /// <summary>
            /// Id of the attribute decoder that was used to decode this attribute data.
            /// </summary>
            internal int decoderId = -1;
            internal MeshAttributeCornerTable connectivityData;
            /// <summary>
            /// Flag that can mark the connectivityData invalid. In such case the base
            /// corner table of the mesh should be used instead.
            /// </summary>
            internal bool isConnectivityUsed = true;
            internal MeshAttributeIndicesEncodingData encodingData = new MeshAttributeIndicesEncodingData();
            /// <summary>
            /// Opposite corners to attribute seam edges.
            /// </summary>
            internal IntList attributeSeamCorners = new IntList();
        }

        private AttributeData[] attributeData;

        private ITraversalDecoder traversalDecoder;

        public MeshEdgeBreakerDecoderImpl(MeshEdgeBreakerDecoder decoder, ITraversalDecoder traversalDecoder)
        {
            Init(decoder);
            this.traversalDecoder = traversalDecoder;
        }

        public void Init(MeshEdgeBreakerDecoder decoder)
        {
            this.decoder = decoder;
        }

        public MeshAttributeCornerTable GetAttributeCornerTable(int attId)
        {

            for (int i = 0; i < attributeData.Length; ++i)
            {
                AttributesDecoder dec = decoder.AttributesDecoders[attributeData[i].decoderId];
                for (int j = 0; j < dec.NumAttributes; ++j)
                {
                    if (dec.GetAttributeId(j) == attId)
                    {
                        if (attributeData[i].isConnectivityUsed)
                            return attributeData[i].connectivityData;
                        return null;
                    }
                }
            }
            return null;
        }

        public MeshAttributeIndicesEncodingData GetAttributeEncodingData(int attId)
        {

            for (int i = 0; i < attributeData.Length; ++i)
            {
                AttributesDecoder dec = decoder.AttributesDecoders[attributeData[i].decoderId];
                for (int j = 0; j < dec.NumAttributes; ++j)
                {
                    if (dec.GetAttributeId(j) == attId)
                        return attributeData[i].encodingData;
                }
            }
            return posEncodingData;
        }

        public void CreateAttributesDecoder(int attDecoderId)
        {

            sbyte attDataId = decoder.Buffer.DecodeI8();
            byte decoderType = decoder.Buffer.DecodeU8();

            if (attDataId >= 0)
            {
                if (attDataId >= attributeData.Length)
                {
                    throw DracoUtils.Failed(); // Unexpected attribute data.
                }
                attributeData[attDataId].decoderId = attDecoderId;
            }
            else
            {
                if (posDataDecoderId >= 0)
                    throw DracoUtils.Failed();
                posDataDecoderId = attDecoderId;
            }


            MeshTraversalMethod traversalMethod = MeshTraversalMethod.DepthFirst;
            if (decoder.BitstreamVersion >= 12)
            {
                byte encoded = decoder.Buffer.DecodeU8();
                traversalMethod = (MeshTraversalMethod)encoded;
            }



            DracoMesh mesh = decoder.Mesh;
            PointsSequencer sequencer;

            if (decoderType == (byte) MeshAttributeElementType.Vertex)
            {
                // Per-vertex attribute decoder.
                // Traverser that is used to generate the encoding order of each attribute.

                MeshAttributeIndicesEncodingData encodingData = null;
                if (attDataId < 0)
                {
                    encodingData = posEncodingData;
                }
                else
                {
                    encodingData = attributeData[attDataId].encodingData;
                    // Mark the attribute connectivity data invalid to ensure it's not used
                    // later on.
                    attributeData[attDataId].isConnectivityUsed = false;
                }

                var traversalSequencer = new MeshTraversalSequencer<CornerTable>(mesh, encodingData);

                var attObserver = new MeshAttributeIndicesEncodingObserver<CornerTable>(cornerTable, mesh,
                    traversalSequencer, encodingData);
                ICornerTableTraverser<CornerTable> attTraverser;
                if (traversalMethod == MeshTraversalMethod.DepthFirst)
                    attTraverser = new EdgeBreakerTraverser<CornerTable>(attObserver);
                else if (traversalMethod == MeshTraversalMethod.PredictionDegree)
                    attTraverser = new PredictionDegreeTraverser(attObserver);
                else
                    throw DracoUtils.Failed();

                traversalSequencer.SetTraverser(attTraverser);
                sequencer = traversalSequencer;

            }
            else
            {
                if (attDataId < 0)
                    throw DracoUtils.Failed(); // Attribute data must be specified.

                // Per-corner attribute decoder.
                //typedef CornerTableTraversalProcessor<MeshAttributeCornerTable> AttProcessor;
                //typedef MeshAttributeIndicesEncodingObserver<MeshAttributeCornerTable> AttObserver;
                // Traverser that is used to generate the encoding order of each attribute.
                //typedef EdgeBreakerTraverser<AttProcessor, AttObserver> AttTraverser;

                MeshTraversalSequencer<MeshAttributeCornerTable> traversalSequencer =
                    new MeshTraversalSequencer<MeshAttributeCornerTable>(mesh,
                        attributeData[attDataId].encodingData);

                var attObserver =
                    new MeshAttributeIndicesEncodingObserver<MeshAttributeCornerTable>(
                        attributeData[attDataId].connectivityData,
                        mesh, traversalSequencer,
                        attributeData[attDataId].encodingData);
                var attTraverser = new EdgeBreakerTraverser<MeshAttributeCornerTable>(attObserver);

                traversalSequencer.SetTraverser(attTraverser);
                sequencer = traversalSequencer;
            }

            var attController = new SequentialAttributeDecodersController(sequencer);

            decoder.SetAttributesDecoder(attDecoderId, attController);
        }

        public void DecodeConnectivity()
        {
            numNewVertices = 0;
            newToParentVertexMap.Clear();
            if (decoder.BitstreamVersion < 22)
            {
                uint num_new_verts;
                if (decoder.BitstreamVersion < 20)
                {
                    num_new_verts = decoder.Buffer.DecodeU32();
                }
                else
                {
                    num_new_verts = Decoding.DecodeVarintU32(decoder.Buffer);
                }

                numNewVertices = (int) num_new_verts;
            }

            uint numEncodedVertices;
            if (decoder.BitstreamVersion < 20)
            {
                numEncodedVertices = decoder.Buffer.DecodeU32();
            }
            else
            {
                numEncodedVertices = Decoding.DecodeVarintU32(decoder.Buffer);
            }

            this.numEncodedVertices = (int)numEncodedVertices;

            uint numFaces;

            if (decoder.BitstreamVersion < 20)
            {
                numFaces = decoder.Buffer.DecodeU32();
            }
            else
            {
                numFaces = Decoding.DecodeVarintU32(decoder.Buffer);
            }

            if (numFaces > 805306367) //Upper limit of int32_t / 3
                throw DracoUtils.Failed();  // Draco cannot handle this many faces.

            if (numEncodedVertices > numFaces * 3)
            {
                throw DracoUtils.Failed(); // There cannot be more vertices than 3 * numFaces.
            }

            byte numAttributeData = decoder.Buffer.DecodeU8();

            uint numEncodedSymbols;
            if (decoder.BitstreamVersion < 20)
            {
                numEncodedSymbols = decoder.Buffer.DecodeU32();
            }
            else
            {
                numEncodedSymbols = Decoding.DecodeVarintU32(decoder.Buffer);
            }

            if (numFaces < numEncodedSymbols)
            {
                // Number of faces needs to be the same or greater than the number of
                // symbols (it can be greater because the initial face may not be encoded as
                // a symbol).
                throw DracoUtils.Failed();
            }

            uint max_encoded_faces = numEncodedSymbols + (numEncodedSymbols / 3);
            if (numFaces > max_encoded_faces)
            {
                // Faces can only be 1 1/3 times bigger than number of encoded symbols. This
                // could only happen if all new encoded components started with interior
                // triangles. E.g. A mesh with multiple tetrahedrons.
                throw DracoUtils.Failed();
            }

            uint numEncodedSplitSymbols;
            if (decoder.BitstreamVersion < 20)
            {
                numEncodedSplitSymbols = decoder.Buffer.DecodeU32();
            }
            else
            {
                numEncodedSplitSymbols = Decoding.DecodeVarintU32(decoder.Buffer);
            }

            if (numEncodedSplitSymbols > numEncodedSymbols)
            {
                throw DracoUtils.Failed(); // Split symbols are a sub-set of all symbols.
            }



            // Decode topology (connectivity).
            vertexTraversalLength.Clear();
            cornerTable = new CornerTable();
            processedCornerIds.Clear();
            processedCornerIds.Capacity = (int)numFaces;
            processedConnectivityCorners.Clear();
            processedConnectivityCorners.Capacity = (int)numFaces;
            topologySplitData.Clear();
            holeEventData.Clear();
            initFaceConfigurations.Clear();
            initCorners.Clear();

            numProcessedHoleEvents = 0;
            lastSymbolId = -1;

            lastFaceId = -1;
            lastVertId = -1;

            attributeData = new AttributeData[numAttributeData];


            cornerTable.Reset((int)numFaces, (int)(numEncodedVertices + numEncodedSplitSymbols));

            // Add one attribute data for each attribute decoder.
            for (int i = 0; i < attributeData.Length; i++)
            {
                attributeData[i] = new AttributeData();
            }

            // Start with all vertices marked as holes (boundaries).
            // Only vertices decoded with TOPOLOGYC symbol (and the initial face) will
            // be marked as non hole vertices. We need to allocate the array larger
            // because split symbols can create extra vertices during the decoding
            // process (these extra vertices are then eliminated during deduplication).
            isVertHole = new bool[numEncodedVertices + numEncodedSplitSymbols];
            for (int i = 0; i < isVertHole.Length; i++)
                isVertHole[i] = true;

            int topologySplitDecodedBytes = -1;
            if (decoder.BitstreamVersion < 22)
            {
                uint encodedConnectivitySize;
                if (decoder.BitstreamVersion < 20)
                {
                    encodedConnectivitySize = decoder.Buffer.DecodeU32();
                }
                else
                {
                    encodedConnectivitySize = Decoding.DecodeVarintU32(decoder.Buffer);
                }

                if (encodedConnectivitySize == 0 ||
                    encodedConnectivitySize > decoder.Buffer.RemainingSize)
                    throw DracoUtils.Failed();
                DecoderBuffer eventBuffer = decoder.Buffer.SubBuffer((int) encodedConnectivitySize); // new DecoderBuffer(buf, decoder.Buffer. decoder.Buffer.DecodedSize + encodedConnectivitySize, decoder.Buffer.BufferSize - decoder.Buffer.DecodedSize - encodedConnectivitySize);

                // Decode hole and topology split events.
                topologySplitDecodedBytes =
                    DecodeHoleAndTopologySplitEvents(eventBuffer);
                if (topologySplitDecodedBytes == -1)
                    throw DracoUtils.Failed();
            }
            else
            {
                if (DecodeHoleAndTopologySplitEvents(decoder.Buffer) == -1)
                    throw DracoUtils.Failed();
            }

            traversalDecoder.Init(this);
            traversalDecoder.SetNumEncodedVertices((int)(numEncodedVertices + numEncodedSplitSymbols));
            traversalDecoder.SetNumAttributeData(numAttributeData);

            DecoderBuffer traversalEndBuffer = traversalDecoder.Start();

            int numConnectivityVerts = DecodeConnectivity((int)numEncodedSymbols);
            if (numConnectivityVerts == -1)
                throw DracoUtils.Failed();

            // Set the main buffer to the end of the traversal.
            decoder.Buffer = traversalEndBuffer.SubBuffer(0);// .Initialize(traversalEndBuffer.GetBuffer(), traversalEndBuffer.  traversalEndBuffer.remainingCount);

            // Skip topology split data that was already decoded earlier.

            if (decoder.BitstreamVersion < 22)
            {
                // Skip topology split data that was already decoded earlier.
                decoder.Buffer.Advance(topologySplitDecodedBytes);
            }

            // Decode connectivity of non-position attributes.
            if (attributeData.Length > 0)
            {
                if (decoder.BitstreamVersion < 21)
                {
                    for (int ci = 0; ci < cornerTable.NumCorners; ci += 3)
                    {
                        DecodeAttributeConnectivitiesOnFaceLegacy(ci);
                    }
                }
                else
                {
                    for (int ci = 0; ci < cornerTable.NumCorners; ci += 3)
                    {
                        DecodeAttributeConnectivitiesOnFace(ci);
                    }
                }
            }
            traversalDecoder.Done();

            // Decode attribute connectivity.
            // Prepare data structure for decoding non-position attribute connectivites.
            for (int i = 0; i < attributeData.Length; ++i)
            {
                attributeData[i].connectivityData = new MeshAttributeCornerTable(cornerTable);
                // Add all seams.
                var corners = attributeData[i].attributeSeamCorners;
                for(int j = 0; j < corners.Count; j++)
                {
                    int c = corners[j];
                    attributeData[i].connectivityData.AddSeamEdge(c);
                }
                // Recompute vertices from the newly added seam edges.
                attributeData[i].connectivityData.RecomputeVertices(null, null);
            }

            //A3DUtils.Resize(posEncodingData.vertexToEncodedAttributeValueIndexMap, cornerTable.NumVertices);
            posEncodingData.vertexToEncodedAttributeValueIndexMap = new int[cornerTable.NumVertices];
            for (int i = 0; i < attributeData.Length; ++i)
            {
                // For non-position attributes, preallocate the vertex to value mapping
                // using the maximum number of vertices from the base corner table and the
                // attribute corner table (since the attribute decoder may use either of
                // it).
                int attConnectivityVerts =
                    attributeData[i].connectivityData.NumVertices;
                if (attConnectivityVerts < cornerTable.NumVertices)
                    attConnectivityVerts = cornerTable.NumVertices;
                //A3DUtils.Resize(attributeData[i].encodingData.vertexToEncodedAttributeValueIndexMap, attConnectivityVerts);
                attributeData[i].encodingData.vertexToEncodedAttributeValueIndexMap = new int[attConnectivityVerts];
            }
            if (!AssignPointsToCorners())
                throw DracoUtils.Failed();
        }

        public void OnAttributesDecoded()
        {
        }

        public MeshEdgeBreakerDecoder GetDecoder()
        {
            return decoder;
        }

        public CornerTable CornerTable
        {
            get { return cornerTable; }
        }

        private bool AssignPointsToCorners()
        {
            // Map between the existing and deduplicated point ids.
            // Note that at this point we have one point id for each corner of the
            // mesh so there is cornerTable.numCorners() point ids.
            decoder.Mesh.NumFaces = cornerTable.NumFaces;
            int[] face = new int[3];

            if (attributeData.Length == 0)
            {
                // We have position only. In this case we can simplify the deduplication
                // because the only thing we need to do is to remove isolated vertices that
                // were introduced during the decoding.

                int numPoints = 0;
                int[] vertexToPointMap;
                vertexToPointMap = new int[cornerTable.NumVertices];
                DracoUtils.Fill(vertexToPointMap, -1);
                // Add faces.
                for (int f = 0; f < decoder.Mesh.NumFaces; ++f)
                {
                    for (int c = 0; c < 3; ++c)
                    {
                        // Remap old points to the new ones.
                        int vertId =
                            cornerTable.Vertex(3*f + c);
                        int pointId = vertexToPointMap[vertId];
                        if (pointId == -1)
                            pointId = vertexToPointMap[vertId] = numPoints++;
                        face[c] = pointId;
                    }
                    decoder.Mesh.SetFace(f, face);
                }
                decoder.PointCloud.NumPoints = numPoints;
                return true;
            }
            // Else we need to deduplicate multiple attributes.

            // Map between point id and an associated corner id. Only one corner for
            // each point is stored. The corners are used to sample the attribute values
            // in the last stage of the deduplication.
            IntList pointToCornerMap = new IntList();
            // Map between every corner and their new point ids.
            int[] cornerToPointMap = new int[cornerTable.NumCorners];// A3DUtils.NewArray<int>(cornerTable.NumCorners, 0);
            for (int v = 0; v < cornerTable.NumVertices; ++v)
            {
                int c = cornerTable.LeftMostCorner(v);
                if (c < 0)
                    continue; // Isolated vertex.
                int deduplicationFirstCorner = c;
                if (isVertHole[v])
                {
                    // If the vertex is on a boundary, start deduplication from the left most
                    // corner that is guaranteed to lie on the boundary.
                    deduplicationFirstCorner = c;
                }
                else
                {
                    // If we are not on the boundary we need to find the first seam (of any
                    // attribute).
                    for (int i = 0; i < attributeData.Length; ++i)
                    {
                        if (!attributeData[i].connectivityData.IsCornerOnSeam(c))
                            continue; // No seam for this attribute, ignore it.
                        // Else there needs to be at least one seam edge.

                        // At this point, we use identity mapping between corners and point ids.
                        int vertId =
                            attributeData[i].connectivityData.Vertex(c);
                        int actC = cornerTable.SwingRight(c);
                        bool seamFound = false;
                        while (actC != c)
                        {
                            if (attributeData[i].connectivityData.Vertex(actC) != vertId)
                            {
                                // Attribute seam found. Stop.
                                deduplicationFirstCorner = actC;
                                seamFound = true;
                                break;
                            }
                            actC = cornerTable.SwingRight(actC);
                        }
                        if (seamFound)
                            break; // No reason to process other attributes if we found a seam.
                    }
                }

                // Do a deduplication pass over the corners on the processed vertex.
                // At this point each corner corresponds to one point id and our goal is to
                // merge similar points into a single point id.
                // We do one one pass in a clocwise direction over the corners and we add
                // a new point id whenever one of the attributes change.
                c = deduplicationFirstCorner;
                // Create a new point.
                cornerToPointMap[c] = pointToCornerMap.Count;
                pointToCornerMap.Add(c);
                // Traverse in CW direction.
                int prevC = c;
                c = cornerTable.SwingRight(c);
                while (c >= 0 && c != deduplicationFirstCorner)
                {
                    bool attributeSeam = false;
                    for (int i = 0; i < attributeData.Length; ++i)
                    {
                        if (attributeData[i].connectivityData.Vertex(c) !=
                            attributeData[i].connectivityData.Vertex(prevC))
                        {
                            // Attribute index changed from the previous corner. We need to add a
                            // new point here.
                            attributeSeam = true;
                            break;
                        }
                    }
                    if (attributeSeam)
                    {
                        cornerToPointMap[c] = pointToCornerMap.Count;
                        pointToCornerMap.Add(c);
                    }
                    else
                    {
                        cornerToPointMap[c] = cornerToPointMap[prevC];
                    }
                    prevC = c;
                    c = cornerTable.SwingRight(c);
                }
            }
            // Add faces.
            for (int f = 0; f < decoder.Mesh.NumFaces; ++f)
            {
                for (int c = 0; c < 3; ++c)
                {
                    // Remap old points to the new ones.
                    face[c] = cornerToPointMap[3*f + c];
                }
                decoder.Mesh.SetFace(f, face);
            }
            decoder.PointCloud.NumPoints = pointToCornerMap.Count;
            return true;
        }

        int DecodeHoleAndTopologySplitEvents(DecoderBuffer decoderBuffer)
        {
            // Prepare a new decoder from the provided buffer offset.
            uint numTopologySplits;
            if (decoder.BitstreamVersion < 20)
            {
                numTopologySplits = decoderBuffer.DecodeU32();
            }
            else
            {
                numTopologySplits = Decoding.DecodeVarintU32(decoderBuffer);
            }

            if (numTopologySplits > 0)
            {
                if (numTopologySplits > cornerTable.NumFaces)
                    return -1;

                if (decoder.BitstreamVersion < 12)
                {
                    for (int i = 0; i < numTopologySplits; ++i)
                    {
                        TopologySplitEventData eventData = new TopologySplitEventData();
                        eventData.splitSymbolId = decoderBuffer.DecodeI32();
                        eventData.sourceSymbolId = decoderBuffer.DecodeI32();
                        byte edgeData = decoderBuffer.DecodeU8();
                        eventData.sourceEdge = (byte)(edgeData & 1);
                        eventData.splitEdge = (byte)((edgeData >> 1) & 1);
                        topologySplitData.Add(eventData);
                    }
                }
                else
                {

                    // Decode source and split symbol ids using delta and varint coding. See
                    // description in mesh_edgebreaker_encoder_impl.cc for more details.
                    int last_source_symbol_id = 0;
                    for (int i = 0; i < numTopologySplits; ++i)
                    {
                        TopologySplitEventData event_data = new TopologySplitEventData();
                        uint delta = Decoding.DecodeVarintU32(decoderBuffer);
                        event_data.sourceSymbolId = (int)(delta + last_source_symbol_id);
                        delta = Decoding.DecodeVarintU32(decoderBuffer);
                        if (delta > event_data.sourceSymbolId)
                            return -1;
                        event_data.splitSymbolId = event_data.sourceSymbolId - (int)(delta);
                        last_source_symbol_id = event_data.sourceSymbolId;
                        topologySplitData.Add(event_data);
                    }

                    // Split edges are decoded from a direct bit decoder.
                    long tmp = decoderBuffer.StartBitDecoding(false);
                    for (int i = 0; i < numTopologySplits; ++i)
                    {
                        uint edge_data;
                        if (decoder.BitstreamVersion < 22)
                        {
                            decoderBuffer.DecodeLeastSignificantBits32(2, out edge_data);
                        }
                        else
                        {
                            decoderBuffer.DecodeLeastSignificantBits32(1, out edge_data);
                        }

                        topologySplitData[i].sourceEdge = (byte)(edge_data & 1);
                    }

                    decoderBuffer.EndBitDecoding();
                }
            }

            uint numHoleEvents = 0;
            if (decoder.BitstreamVersion < 20)
            {
                numHoleEvents = decoderBuffer.DecodeU32();
            }
            else if(decoder.BitstreamVersion < 21)
            {
                numHoleEvents = Decoding.DecodeVarintU32(decoderBuffer);
            }

            if (numHoleEvents > 0)
            {
                if (decoder.BitstreamVersion < 12)
                {
                    for (uint i = 0; i < numHoleEvents; ++i)
                    {
                        HoleEventData eventData = new HoleEventData();
                        eventData.symbolId = decoderBuffer.DecodeI32();
                        holeEventData.Add(eventData);
                    }
                }
                else
                {

      // Decode hole symbol ids using delta and varint coding.
                    int last_symbol_id = 0;
                    for (int i = 0; i < numHoleEvents; ++i)
                    {
                        HoleEventData event_data = new HoleEventData();
                        uint delta = Decoding.DecodeVarintU32(decoderBuffer);
                        event_data.symbolId = (int)(delta + last_symbol_id);
                        last_symbol_id = event_data.symbolId;
                        holeEventData.Add(event_data);
                    }
                }
            }

            return decoderBuffer.DecodedSize;
        }

        private int DecodeConnectivity(int numSymbols)
        {
            // Algorithm does the reverse decoding of the symbols encoded with the
            // edgebreaker method. The reverse decoding always keeps track of the active
            // edge identified by its opposite corner (active corner). New faces are
            // always added to this active edge. There may be multiple active corners at
            // one time that either correspond to separate mesh components or to
            // sub-components of one mesh that are going to be merged together using the
            // TOPOLOGYS symbol. We can store these active edges on a stack, because the
            // decoder always processes only the latest active edge. TOPOLOGYS then
            // removes the top edge from the stack and TOPOLOGYE adds a new edge to the
            // stack.
            IntList activeCornerStack = new IntList();

            // Additional active edges may be added as a result of topology split events.
            // They can be added in arbitrary order, but we always know the split symbol
            // id they belong to, so we can address them using this symbol id.
            Dictionary<int, int> topologySplitActiveCorners = new Dictionary<int, int>();
            bool removeInvalidVertices = attributeData.Length == 0;
            IntList invalidVertices = new IntList();

            //int numVertices = 0;
            int maxNumVertices = isVertHole.Length;
            int numFaces = 0;
            for (int symbolId = 0; symbolId < numSymbols; ++symbolId)
            {
                int face = numFaces++;
                // Used to flag cases where we need to look for topology split events.
                bool checkTopologySplit = false;
                EdgeBreakerTopologyBitPattern symbol = traversalDecoder.DecodeSymbol();
                if (symbol == EdgeBreakerTopologyBitPattern.C)
                {
                    // Create a new face between two edges on the open boundary.
                    // The first edge is opposite to the corner "a" from the image below.
                    // The other edge is opposite to the corner "b" that can be reached
                    // through a CCW traversal around the vertex "v".
                    // One new active boundary edge is created, opposite to the new corner
                    // "x".
                    //
                    //     *-------*
                    //    / \     / \
                    //   /   \   /   \
                    //  /     \ /     \
                    // *-------v-------*
                    //  \b    /x\    a/
                    //   \   /   \   /
                    //    \ /  C  \ /
                    //     *.......*
                    // Find the corner "b" from the corner "a" which is the corner on the
                    // top of the active stack.
                    if (activeCornerStack.Count == 0)
                        return -1;

                    int cornerA = activeCornerStack[activeCornerStack.Count - 1];
                    int vertexX = cornerTable.Vertex(cornerTable.Next(cornerA));
                    int cornerB =
                        cornerTable.Next(cornerTable.LeftMostCorner(vertexX));

                    // New tip corner.
                    int corner = 3 * face;
                    // Update opposite corner mappings.
                    SetOppositeCorners(cornerA, corner + 1);
                    SetOppositeCorners(cornerB, corner + 2);

                    // Update vertex mapping.
                    cornerTable.MapCornerToVertex(corner, vertexX);
                    cornerTable.MapCornerToVertex(
                        corner + 1, cornerTable.Vertex(cornerTable.Next(cornerB)));
                    int vert_a_prev =
                        cornerTable.Vertex(cornerTable.Previous(cornerA));
                    cornerTable.MapCornerToVertex(corner + 2, vert_a_prev);
                    cornerTable.SetLeftMostCorner(vert_a_prev, corner + 2);
                    // Mark the vertex |x| as interior.
                    isVertHole[vertexX] = false;
                    // Update the corner on the active stack.
                    activeCornerStack[activeCornerStack.Count - 1] = corner;
                }
                else if (symbol == EdgeBreakerTopologyBitPattern.R || symbol == EdgeBreakerTopologyBitPattern.L)
                {

                    // Create a new face extending from the open boundary edge opposite to the
                    // corner "a" from the image below. Two new boundary edges are created
                    // opposite to corners "r" and "l". New active corner is set to either "r"
                    // or "l" depending on the decoded symbol. One new vertex is created
                    // at the opposite corner to corner "a".
                    //     *-------*
                    //    /a\     / \
                    //   /   \   /   \
                    //  /     \ /     \
                    // *-------v-------*
                    //  .l   r.
                    //   .   .
                    //    . .
                    //     *
                    if (activeCornerStack.Count == 0)
                        return -1;
                    int corner_a = activeCornerStack[activeCornerStack.Count - 1];
                    // First corner on the new face is either corner "l" or "r".
                    int corner = 3 * face;
                    int opp_corner, corner_l, corner_r;
                    if (symbol == EdgeBreakerTopologyBitPattern.R)
                    {
                        // "r" is the new first corner.
                        opp_corner = corner + 2;
                        corner_l = corner + 1;
                        corner_r = corner;
                    }
                    else
                    {
                        // "l" is the new first corner.
                        opp_corner = corner + 1;
                        corner_l = corner;
                        corner_r = corner + 2;
                    }

                    SetOppositeCorners(opp_corner, corner_a);
                    // Update vertex mapping.
                    int new_vert_index = cornerTable.AddNewVertex();

                    if (cornerTable.NumVertices > maxNumVertices)
                        return -1; // Unexpected number of decoded vertices.

                    cornerTable.MapCornerToVertex(opp_corner, new_vert_index);
                    cornerTable.SetLeftMostCorner(new_vert_index, opp_corner);

                    int vertex_r =
                        cornerTable.Vertex(cornerTable.Previous(corner_a));
                    cornerTable.MapCornerToVertex(corner_r, vertex_r);
                    // Update left-most corner on the vertex on the |corner_r|.
                    cornerTable.SetLeftMostCorner(vertex_r, corner_r);

                    cornerTable.MapCornerToVertex(
                        corner_l, cornerTable.Vertex(cornerTable.Next(corner_a)));
                    activeCornerStack[activeCornerStack.Count - 1] = corner;
                    checkTopologySplit = true;
                }
                else if (symbol == EdgeBreakerTopologyBitPattern.S)
                {
                    // Create a new face that merges two last active edges from the active
                    // stack. No new vertex is created, but two vertices at corners "p" and
                    // "n" need to be merged into a single vertex.
                    //
                    // *-------v-------*
                    //  \a   p/x\n   b/
                    //   \   /   \   /
                    //    \ /  S  \ /
                    //     *.......*
                    //
                    if (activeCornerStack.Count == 0)
                        return -1;
                    int corner_b = activeCornerStack[activeCornerStack.Count - 1];
                    activeCornerStack.RemoveAt(activeCornerStack.Count - 1);

                    // Corner "a" can correspond either to a normal active edge, or to an edge
                    // created from the topology split event.
                    int tmp;
                    if (topologySplitActiveCorners.TryGetValue(symbolId, out tmp))
                    {
                        // Topology split event. Move the retrieved edge to the stack.
                        activeCornerStack.Add(tmp);
                    }

                    if (activeCornerStack.Count == 0)
                        return -1;
                    int corner_a = activeCornerStack[activeCornerStack.Count - 1];

                    if (cornerTable.Opposite(corner_a) != CornerTable.kInvalidCornerIndex ||
                        cornerTable.Opposite(corner_b) != CornerTable.kInvalidCornerIndex)
                    {
                        // One of the corners is already opposite to an existing face, which
                        // should not happen unless the input was tempered with.
                        return -1;
                    }

                    // First corner on the new face is corner "x" from the image above.
                    int corner = 3 * face;
                    // Update the opposite corner mapping.
                    SetOppositeCorners(corner_a, corner + 2);
                    SetOppositeCorners(corner_b, corner + 1);
                    // Update vertices. For the vertex at corner "x", use the vertex id from
                    // the corner "p".
                    int vertex_p =
                        cornerTable.Vertex(cornerTable.Previous(corner_a));
                    cornerTable.MapCornerToVertex(corner, vertex_p);
                    cornerTable.MapCornerToVertex(
                        corner + 1, cornerTable.Vertex(cornerTable.Next(corner_a)));
                    int vert_b_prev =
                        cornerTable.Vertex(cornerTable.Previous(corner_b));
                    cornerTable.MapCornerToVertex(corner + 2, vert_b_prev);
                    cornerTable.SetLeftMostCorner(vert_b_prev, corner + 2);
                    int corner_n = cornerTable.Next(corner_b);
                    int vertex_n = cornerTable.Vertex(corner_n);
                    traversalDecoder.MergeVertices(vertex_p, vertex_n);
                    // Update the left most corner on the newly merged vertex.
                    cornerTable.SetLeftMostCorner(vertex_p,
                        cornerTable.LeftMostCorner(vertex_n));

                    // Also update the vertex id at corner "n" and all corners that are
                    // connected to it in the CCW direction.
                    while (corner_n != CornerTable.kInvalidCornerIndex)
                    {
                        cornerTable.MapCornerToVertex(corner_n, vertex_p);
                        corner_n = cornerTable.SwingLeft(corner_n);
                    }

                    // Make sure the old vertex n is now mapped to an invalid corner (make it
                    // isolated).
                    cornerTable.MakeVertexIsolated(vertex_n);
                    if (removeInvalidVertices)
                        invalidVertices.Add(vertex_n);
                    activeCornerStack[activeCornerStack.Count - 1] = corner;
                }
                else if (symbol == EdgeBreakerTopologyBitPattern.E)
                {

                    int corner = 3 * face;
                    int firstVertIdx = cornerTable.AddNewVertex();
                    // Create three new vertices at the corners of the new face.
                    cornerTable.MapCornerToVertex(corner, firstVertIdx);
                    cornerTable.MapCornerToVertex(corner + 1, cornerTable.AddNewVertex());
                    cornerTable.MapCornerToVertex(corner + 2, cornerTable.AddNewVertex());

                    if (cornerTable.NumVertices > maxNumVertices)
                        return -1; // Unexpected number of decoded vertices.

                    cornerTable.SetLeftMostCorner(firstVertIdx, corner);
                    cornerTable.SetLeftMostCorner(firstVertIdx + 1, corner + 1);
                    cornerTable.SetLeftMostCorner(firstVertIdx + 2, corner + 2);
                    // Add the tip corner to the active stack.
                    activeCornerStack.Add(corner);
                    checkTopologySplit = true;
                }
                else
                {
                    // Error. Unknown symbol decoded.
                    return -1;
                }

                // Inform the traversal decoder that a new corner has been reached.
                traversalDecoder.NewActiveCornerReached(activeCornerStack[activeCornerStack.Count - 1]);

                if (checkTopologySplit)
                {
                    // Check for topology splits happens only for TOPOLOGY_L, TOPOLOGY_R and
                    // TOPOLOGY_E symbols because those are the symbols that correspond to
                    // faces that can be directly connected a TOPOLOGY_S face through the
                    // topology split event.
                    // If a topology split is detected, we need to add a new active edge
                    // onto the activeCornerStack because it will be used later when the
                    // corresponding TOPOLOGY_S event is decoded.

                    // Symbol id used by the encoder (reverse).
                    int encoder_symbol_id = numSymbols - symbolId - 1;
                    EdgeFaceName split_edge;
                    int encoderSplitSymbolId;
                    while (IsTopologySplit(encoder_symbol_id, out split_edge, out encoderSplitSymbolId))
                    {
                        if (encoderSplitSymbolId < 0)
                            return -1; // Wrong split symbol id.
                        // Symbol was part of a topology split. Now we need to determine which
                        // edge should be added to the active edges stack.
                        int act_top_corner = activeCornerStack.Back;
                        // The current symbol has one active edge (stored in act_top_corner) and
                        // two remaining inactive edges that are attached to it.
                        //              *
                        //             / \
                        //  left_edge /   \ right_edge
                        //           /     \
                        //          *.......*
                        //         active_edge

                        int new_active_corner;
                        if (split_edge == EdgeFaceName.RightFaceEdge)
                        {
                            new_active_corner = cornerTable.Next(act_top_corner);
                        }
                        else
                        {
                            new_active_corner = cornerTable.Previous(act_top_corner);
                        }

                        // Add the new active edge.
                        // Convert the encoder split symbol id to decoder symbol id.
                        int decoderSplitSymbolId =
                            numSymbols - encoderSplitSymbolId - 1;
                        topologySplitActiveCorners[decoderSplitSymbolId] =
                            new_active_corner;
                    }
                }
            }

            if (cornerTable.NumVertices > maxNumVertices)
                return -1; // Unexpected number of decoded vertices.
            // Decode start faces and connect them to the faces from the active stack.
            while (activeCornerStack.Count > 0)
            {
                int corner = activeCornerStack.Back;
                activeCornerStack.PopBack();
                bool interior_face =
                    traversalDecoder.DecodeStartFaceConfiguration();
                if (interior_face)
                {
                    // The start face is interior, we need to find three corners that are
                    // opposite to it. The first opposite corner "a" is the corner from the
                    // top of the active corner stack and the remaining two corners "b" and
                    // "c" are then the next corners from the left-most corners of vertices
                    // "n" and "x" respectively.
                    //
                    //           *-------*
                    //          / \     / \
                    //         /   \   /   \
                    //        /     \ /     \
                    //       *-------p-------*
                    //      / \a    . .    c/ \
                    //     /   \   .   .   /   \
                    //    /     \ .  I  . /     \
                    //   *-------n.......x------*
                    //    \     / \     / \     /
                    //     \   /   \   /   \   /
                    //      \ /     \b/     \ /
                    //       *-------*-------*
                    //

                    if (numFaces >= cornerTable.NumFaces)
                    {
                        return -1; // More faces than expected added to the mesh.
                    }

                    int corner_a = corner;
                    int vert_n = cornerTable.Vertex(cornerTable.Next(corner_a));
                    int corner_b = cornerTable.Next(cornerTable.LeftMostCorner(vert_n));

                    int vert_x = cornerTable.Vertex(cornerTable.Next(corner_b));
                    int corner_c = cornerTable.Next(cornerTable.LeftMostCorner(vert_x));

                    int vert_p = cornerTable.Vertex(cornerTable.Next(corner_c));

                    int face = numFaces ++;
                    // The first corner of the initial face is the corner opposite to "a".
                    int new_corner = 3 * face;
                    SetOppositeCorners(new_corner, corner);
                    SetOppositeCorners(new_corner + 1, corner_b);
                    SetOppositeCorners(new_corner + 2, corner_c);

                    // Map new corners to existing vertices.
                    cornerTable.MapCornerToVertex(new_corner, vert_x);
                    cornerTable.MapCornerToVertex(new_corner + 1, vert_p);
                    cornerTable.MapCornerToVertex(new_corner + 2, vert_n);

                    // Mark all three vertices as interior.
                    for (int ci = 0; ci < 3; ++ci)
                    {
                        isVertHole[cornerTable.Vertex(new_corner + ci)] = false;
                    }

                    initFaceConfigurations.Add(true);
                    initCorners.Add(new_corner);
                }
                else
                {
                    // The initial face wasn't interior and the traversal had to start from
                    // an open boundary. In this case no new face is added, but we need to
                    // keep record about the first opposite corner to this boundary.
                    initFaceConfigurations.Add(false);
                    initCorners.Add(corner);
                }
            }

            if (numFaces != cornerTable.NumFaces)
                return -1; // Unexpected number of decoded faces.

            int num_vertices = cornerTable.NumVertices;
            // If any vertex was marked as isolated, we want to remove it from the corner
            // table to ensure that all vertices in range <0, num_vertices> are valid.
            for(int i = 0; i < invalidVertices.Count; i++) {
                var invalidVert = invalidVertices[i];
                // Find the last valid vertex and swap it with the isolated vertex.
                int srcVert = num_vertices -1;
                while (cornerTable.LeftMostCorner(srcVert) == CornerTable.kInvalidCornerIndex)
                {
                    // The last vertex is invalid, proceed to the previous one.
                    srcVert = --num_vertices - 1;
                }

                if (srcVert < invalidVert)
                    continue; // No need to swap anything.

                // Remap all corners mapped to |src_vert| to |invalid_vert|.
                VertexCornersIterator vcit = VertexCornersIterator.FromVertex(cornerTable, srcVert);
                for (; !vcit.End; vcit.Next())
                {
                    int cid = vcit.Corner;
                    cornerTable.MapCornerToVertex(cid, invalidVert);
                }

                cornerTable.SetLeftMostCorner(invalidVert, cornerTable.LeftMostCorner(srcVert));

                // Make the |src_vert| invalid.
                cornerTable.MakeVertexIsolated(srcVert);
                isVertHole[invalidVert] = isVertHole[srcVert];
                isVertHole[srcVert] = false;

                // The last vertex is now invalid.
                num_vertices--;
            }
            return num_vertices;
        }

        private void SetOppositeCorners(int corner0, int corner1)
        {
            cornerTable.SetOppositeCorner(corner0, corner1);
            cornerTable.SetOppositeCorner(corner1, corner0);
        }

        /// <summary>
        /// Returns true if the current symbol was part of a topolgy split event. This
        /// means that the current face was connected to the left edge of a face
        /// encoded with the TOPOLOGYS symbol. |outSymbolEdge| can be used to
        /// identify which edge of the source symbol was connected to the TOPOLOGYS
        /// symbol.
        /// </summary>
        private bool IsTopologySplit(int encoderSymbolId, out EdgeFaceName outFaceEdge,
            out int outEncoderSplitSymbolId)
        {
            outFaceEdge = EdgeFaceName.LeftFaceEdge;
            outEncoderSplitSymbolId = 0;
            if (topologySplitData.Count == 0)
                return false;
            var back = topologySplitData[topologySplitData.Count - 1];
            if (back.sourceSymbolId > encoderSymbolId)
            {
                // Something is wrong; if the desired source symbol is greater than the
                // current encoderSymbolId, we missed it, or the input was tampered
                // (|encoderSymbolId| keeps decreasing).
                // Return invalid symbol id to notify the decoder that there was an
                // error.
                outEncoderSplitSymbolId = -1;
                return true;
            }
            if (back.sourceSymbolId != encoderSymbolId)
                return false;
            outFaceEdge = (EdgeFaceName)(back.sourceEdge);
            outEncoderSplitSymbolId = back.splitSymbolId;
            // Remove the latest split event.
            topologySplitData.RemoveAt(topologySplitData.Count - 1);
            return true;
        }

        /// <summary>
        /// Decodes all non-position attribute connectivities on the currently
        /// processed face.
        /// </summary>
        private void DecodeAttributeConnectivitiesOnFaceLegacy(int corner)
        {
            // Three corners of the face.
            DecodeAttributeConnectivitiesOnFaceLegacyImpl(corner);
            DecodeAttributeConnectivitiesOnFaceLegacyImpl(cornerTable.Next(corner));
            DecodeAttributeConnectivitiesOnFaceLegacyImpl(cornerTable.Previous(corner));
        }
        private void DecodeAttributeConnectivitiesOnFaceLegacyImpl(int corner)
        {
            // Three corners of the face.
            {
                int oppCorner = cornerTable.Opposite(corner);
                if (oppCorner < 0)
                {
                    // Don't decode attribute seams on boundary edges (every boundary edge
                    // is automatically an attribute seam).
                    for (int i = 0; i < attributeData.Length; ++i)
                    {
                        attributeData[i].attributeSeamCorners.Add(corner);
                    }
                    return;
                }

                for (int i = 0; i < attributeData.Length; ++i)
                {
                    bool isSeam = traversalDecoder.DecodeAttributeSeam(i);
                    if (isSeam)
                        attributeData[i].attributeSeamCorners.Add(corner);
                }
            }
        }

        private void DecodeAttributeConnectivitiesOnFace(int corner)
        {

            // Three corners of the face.
            int src_face_id = cornerTable.Face(corner);
            DecodeAttributeConnectivitiesOnFace(corner, src_face_id);
            DecodeAttributeConnectivitiesOnFace(cornerTable.Next(corner), src_face_id);
            DecodeAttributeConnectivitiesOnFace(cornerTable.Previous(corner), src_face_id);
        }

        private void DecodeAttributeConnectivitiesOnFace(int corner, int src_face_id)
        { 
            int opp_corner = cornerTable.Opposite(corner);
            if (opp_corner == CornerTable.kInvalidCornerIndex)
            {
                // Don't decode attribute seams on boundary edges (every boundary edge
                // is automatically an attribute seam).
                for (int i = 0; i < attributeData.Length; ++i)
                {
                    attributeData[i].attributeSeamCorners.Add(corner);
                }

                return;
            }

            int opp_face_id = cornerTable.Face(opp_corner);
            // Don't decode edges when the opposite face has been already processed.
            if (opp_face_id < src_face_id)
                return;

            for (int i = 0; i < attributeData.Length; ++i)
            {
                bool is_seam = traversalDecoder.DecodeAttributeSeam(i);
                if (is_seam)
                    attributeData[i].attributeSeamCorners.Add(corner);
            }
        }
    }
}
