using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Encoder
{
    /// <summary>
    /// Encoder that tries to predict the edgebreaker traversal symbols based on the
    /// vertex valences of the unencoded portion of the mesh. The current prediction
    /// scheme assumes that each vertex has valence 6 which can be used to predict
    /// the symbol preceeding the one that is currently encoded. Predictions are
    /// encoded using an arithmetic coding which can lead to less than 1 bit per
    /// triangle encoding for highly regular meshes.
    /// </summary>
    class MeshEdgeBreakerTraversalPredictiveEncoder : MeshEdgeBreakerTraversalEncoder
    {
        private CornerTable cornerTable;
        private IntList vertexValences = new IntList();
        private List<bool> predictions = new List<bool>();
        /// <summary>
        /// Previously encoded symbol.
        /// </summary>
        private EdgeBreakerTopologyBitPattern? prevSymbol;
        /// <summary>
        /// The total number of encoded split symbols.
        /// </summary>
        private int numSplitSymbols;
        private int lastCorner;
        /// <summary>
        /// Explicitly count the number of encoded symbols.
        /// </summary>
        private int numSymbols;

        private EdgeBreakerTopologyBitPattern ComputePredictedSymbol(int pivot)
        {
            int valence = vertexValences[pivot];
            if (valence < 0)
            {
                // This situation can happen only for split vertices. Returning
                // TOPOLOGYINVALID always cases misprediction.
                return EdgeBreakerTopologyBitPattern.Invalid;
            }
            if (valence < 6)
            {
                return EdgeBreakerTopologyBitPattern.R;
            }
            return EdgeBreakerTopologyBitPattern.C;
        }

        public override void EncodeSymbol(EdgeBreakerTopologyBitPattern symbol)
        {
            ++numSymbols;
            // Update valences on the mesh. And compute the predicted preceding symbol.
            // Note that the valences are computed for the so far unencoded part of the
            // mesh. Adding a new symbol either reduces valences on the vertices or
            // leaves the valence unchanged.
            EdgeBreakerTopologyBitPattern? predictedSymbol = null;
            int next = cornerTable.Next(lastCorner);
            int prev = cornerTable.Previous(lastCorner);
            switch (symbol)
            {
                case EdgeBreakerTopologyBitPattern.C:
                    // Compute prediction.
                    predictedSymbol = ComputePredictedSymbol(cornerTable.Vertex(next));
                    vertexValences[cornerTable.Vertex(next)] -= 1;
                    vertexValences[cornerTable.Vertex(prev)] -= 1;
                    break;
                case EdgeBreakerTopologyBitPattern.S:
                    // Update velences.
                    vertexValences[cornerTable.Vertex(next)] -= 1;
                    vertexValences[cornerTable.Vertex(prev)] -= 1;
                    // Whenever we reach a split symbol, mark its tip vertex as invalid by
                    // setting the valence to a negative value. Any prediction that will
                    // use this vertex will then cause a misprediction. This is currently
                    // necessary because the decodding works in the reverse direction and
                    // the decoder doesn't know about these vertices until the split
                    // symbol is decoded at which point two vertices are merged into one.
                    // This can be most likely solved on the encoder side by spliting the
                    // tip vertex into two, but since split symbols are relatively rare,
                    // it's probably not worth doing it.
                    vertexValences[cornerTable.Vertex(lastCorner)] = -1;
                    ++numSplitSymbols;
                    break;
                case EdgeBreakerTopologyBitPattern.R:
                    // Compute prediction.
                    predictedSymbol = ComputePredictedSymbol(cornerTable.Vertex(next));
                    // Update valences.
                    vertexValences[cornerTable.Vertex(lastCorner)] -= 1;
                    vertexValences[cornerTable.Vertex(next)] -= 1;
                    vertexValences[cornerTable.Vertex(prev)] -= 2;
                    break;
                case EdgeBreakerTopologyBitPattern.L:
                    vertexValences[cornerTable.Vertex(lastCorner)] -= 1;
                    vertexValences[cornerTable.Vertex(next)] -= 2;
                    vertexValences[cornerTable.Vertex(prev)] -= 1;
                    break;
                case EdgeBreakerTopologyBitPattern.E:
                    vertexValences[cornerTable.Vertex(lastCorner)] -= 2;
                    vertexValences[cornerTable.Vertex(next)] -= 2;
                    vertexValences[cornerTable.Vertex(prev)] -= 2;
                    break;
                default:
                    break;
            }
            // Flag used when it's necessary to explicitly store the previous symbol.
            bool storePrevSymbol = true;
            if (predictedSymbol != null)
            {
                if (prevSymbol != null)
                {
                    if (predictedSymbol.Value == prevSymbol)
                    {
                        predictions.Add(true);
                        storePrevSymbol = false;
                    }
                    else
                    {
                        predictions.Add(false);
                    }
                }
            }
            if (storePrevSymbol && prevSymbol != null)
            {
                base.EncodeSymbol((EdgeBreakerTopologyBitPattern) prevSymbol.Value);
            }
            prevSymbol = symbol;
        }

        public override void NewCornerReached(int corner)
        {
            lastCorner = corner;
        }



        public override void Done()
        {
            // We still need to store the last encoded symbol.
            if (prevSymbol != null)
            {
                base.EncodeSymbol(prevSymbol.Value);
            }
            // Store the init face configurations and the explicitly encoded symbols.
            base.Done();
            // Encode the number of split symbols.
            OutputBuffer.Encode(numSplitSymbols);
            // Store the predictions.
            RAnsBitEncoder predictionEncoder = new RAnsBitEncoder();
            predictionEncoder.StartEncoding();
            for (int i = predictions.Count - 1; i >= 0; --i)
            {
                predictionEncoder.EncodeBit(predictions[i]);
            }
            predictionEncoder.EndEncoding(OutputBuffer);
        }


        public override void Init(IMeshEdgeBreakerEncoder encoder)
        {
            base.Init(encoder);
            cornerTable = encoder.CornerTable;
            // Initialize valences of all vertices.
            vertexValences.Resize(cornerTable.NumVertices, 0);
            for (int i = 0; i < vertexValences.Count; ++i)
            {
                vertexValences[i] = cornerTable.Valence(i);
            }
        }

        public override int NumEncodedSymbols
        {
            get { return numSymbols; }
        }
    }
}
