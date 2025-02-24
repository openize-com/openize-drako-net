using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Decoder
{
    class MeshEdgeBreakerTraversalValenceDecoder : MeshEdgeBreakerTraversalDecoder
    {
        private static readonly EdgeBreakerTopologyBitPattern[] EdgeBreakerSymbolToTopologyId =
        {
            EdgeBreakerTopologyBitPattern.C, EdgeBreakerTopologyBitPattern.S, EdgeBreakerTopologyBitPattern.L, EdgeBreakerTopologyBitPattern.R, EdgeBreakerTopologyBitPattern.E
        };
        CornerTable cornerTable;
        int numVertices = 0;
        private int[] vertexValences;
        private EdgeBreakerTopologyBitPattern lastSymbol = EdgeBreakerTopologyBitPattern.Invalid;
        int activeContext = -1;

        int minValence = 2;
        int maxValence = 7;

        private int[][] contextSymbols;
        // Points to the active symbol in each context.
        private int[] contextCounters;

        public override void Init(IMeshEdgeBreakerDecoderImpl decoder)
        {
            base.Init(decoder);
            cornerTable = decoder.CornerTable;
        }

        public override void SetNumEncodedVertices(int numVertices)
        {
            this.numVertices = numVertices;
        }
        public override DecoderBuffer Start()
        {
            DecoderBuffer outBuffer = null;
            if (BitstreamVersion < 22)
            {
                base.DecodeTraversalSymbols();
            }

            base.DecodeStartFaces();
            base.DecodeAttributeSeams();
            outBuffer = buffer.Clone();

            if (BitstreamVersion < 22)
            {
                int numSplitSymbols;
                if (BitstreamVersion < 20)
                {
                    numSplitSymbols = outBuffer.DecodeI32();
                }
                else
                {
                    numSplitSymbols = (int)outBuffer.DecodeVarintU32();
                }

                if (numSplitSymbols >= numVertices)
                    throw DracoUtils.Failed();

                byte mode = outBuffer.DecodeU8();
                if (mode == 0)// Edgebreaker valence mode 2-7
                {
                    minValence = 2;
                    maxValence = 7;
                }
                else
                {
                    // Unsupported mode.
                    throw DracoUtils.Failed();
                }
            }
            else
            {
                minValence = 2;
                maxValence = 7;
            }

            if (numVertices < 0)
                throw DracoUtils.Failed();
            // Set the valences of all initial vertices to 0.
            vertexValences = new int[numVertices];

            int numUniqueValences = maxValence - minValence + 1;

            // Decode all symbols for all contexts.
            contextSymbols = new int[numUniqueValences][];
            contextCounters = new int[contextSymbols.Length];
            for (int i = 0; i < contextSymbols.Length; ++i)
            {
                uint numSymbols = outBuffer.DecodeVarintU32();
                if (numSymbols > 0)
                {

                    contextSymbols[i] = new int[(int)numSymbols];
                    Decoding.DecodeSymbols((int)numSymbols, 1, outBuffer, contextSymbols[i].AsSpan());
                    // All symbols are going to be processed from the back.
                    contextCounters[i] = (int)numSymbols;
                }
            }

            return outBuffer;
        }

        public override EdgeBreakerTopologyBitPattern DecodeSymbol()
        {
            // First check if we have a valid context.
            if (activeContext != -1)
            {
                int contextCounter = --contextCounters[activeContext];
                if (contextCounter < 0)
                    return EdgeBreakerTopologyBitPattern.Invalid;
                int symbol_id = contextSymbols[activeContext][contextCounter];
                lastSymbol = EdgeBreakerSymbolToTopologyId[symbol_id];
            }
            else
            {
                if (BitstreamVersion < 22)
                {
                    // We don't have a predicted symbol or the symbol was mis-predicted.
                    // Decode it directly.
                    lastSymbol = base.DecodeSymbol();
                }
                else
                {
                    // The first symbol must be E.
                    lastSymbol = EdgeBreakerTopologyBitPattern.E;
                }
            }

            return lastSymbol;
        }

        public override void NewActiveCornerReached(int corner)
        {
            int next = cornerTable.Next(corner);
            int prev = cornerTable.Previous(corner);
            // Update valences.
            switch (lastSymbol)
            {
                case EdgeBreakerTopologyBitPattern.S:
                case EdgeBreakerTopologyBitPattern.C:
                    vertexValences[cornerTable.Vertex(next)] += 1;
                    vertexValences[cornerTable.Vertex(prev)] += 1;
                    break;
                case EdgeBreakerTopologyBitPattern.R:
                    vertexValences[cornerTable.Vertex(corner)] += 1;
                    vertexValences[cornerTable.Vertex(next)] += 1;
                    vertexValences[cornerTable.Vertex(prev)] += 2;
                    break;
                case EdgeBreakerTopologyBitPattern.L:
                    vertexValences[cornerTable.Vertex(corner)] += 1;
                    vertexValences[cornerTable.Vertex(next)] += 2;
                    vertexValences[cornerTable.Vertex(prev)] += 1;
                    break;
                case EdgeBreakerTopologyBitPattern.E:
                    vertexValences[cornerTable.Vertex(corner)] += 2;
                    vertexValences[cornerTable.Vertex(next)] += 2;
                    vertexValences[cornerTable.Vertex(prev)] += 2;
                    break;
                default:
                    break;
            }

            // Compute the new context that is going to be used to decode the next
            // symbol.
            int activeValence = vertexValences[cornerTable.Vertex(next)];
            int clampedValence;
            if (activeValence < minValence)
            {
                clampedValence = minValence;
            }
            else if (activeValence > maxValence)
            {
                clampedValence = maxValence;
            }
            else
            {
                clampedValence = activeValence;
            }

            activeContext = (clampedValence - minValence);
        }

        public override void MergeVertices(int dest, int source)
        {
            // Update valences on the merged vertices.
            vertexValences[dest] += vertexValences[source];
        }
    }
}
