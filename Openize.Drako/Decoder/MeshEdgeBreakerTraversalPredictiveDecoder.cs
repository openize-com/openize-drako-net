using FileFormat.Drako.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Decoder
{

    /// <summary>
    /// Decoder for traversal encoded with the
    /// MeshEdgeBreakerTraversalPredictiveEncoder. The decoder maintains valences
    /// of the decoded portion of the traversed mesh and it uses them to predict
    /// symbols that are about to be decoded.
    /// </summary>
    class MeshEdgeBreakerTraversalPredictiveDecoder : MeshEdgeBreakerTraversalDecoder
    {
        private CornerTable corner_table_;
        private int num_vertices_;
        private int[] vertex_valences_;
        private RAnsBitDecoder prediction_decoder_ = new RAnsBitDecoder();
        private EdgeBreakerTopologyBitPattern last_symbol_ = EdgeBreakerTopologyBitPattern.Invalid;
        private EdgeBreakerTopologyBitPattern predicted_symbol_ = EdgeBreakerTopologyBitPattern.Invalid;

        public override void Init(IMeshEdgeBreakerDecoderImpl decoder)
        {
            base.Init(decoder);
            corner_table_ = decoder.CornerTable;
        }

        public override void SetNumEncodedVertices(int num_vertices)
        {
            num_vertices_ = num_vertices;
        }

        public override DecoderBuffer Start()
        {
            var out_buffer = base.Start();
            int num_split_symbols = out_buffer.DecodeI32();
            // Add one vertex for each split symbol.
            num_vertices_ += num_split_symbols;
            // Set the valences of all initial vertices to 0.
            Array.Resize(ref vertex_valences_, num_vertices_);
            prediction_decoder_.StartDecoding(out_buffer);
            return out_buffer;
        }

        public override EdgeBreakerTopologyBitPattern DecodeSymbol()
        {
            // First check if we have a predicted symbol.
            if (predicted_symbol_ != EdgeBreakerTopologyBitPattern.Invalid)
            {
                // Double check that the predicted symbol was predicted correctly.
                if (prediction_decoder_.DecodeNextBit())
                {
                    last_symbol_ = predicted_symbol_;
                    return predicted_symbol_;
                }
            }
            // We don't have a predicted symbol or the symbol was mis-predicted.
            // Decode it directly.
            last_symbol_ = base.DecodeSymbol();
            return last_symbol_;
        }

        public override void NewActiveCornerReached(int corner)
        {
            int next = corner_table_.Next(corner);
            int prev = corner_table_.Previous(corner);
            // Update valences.
            switch (last_symbol_)
            {
                case EdgeBreakerTopologyBitPattern.C:
                case EdgeBreakerTopologyBitPattern.S:
                    vertex_valences_[corner_table_.Vertex(next)] += 1;
                    vertex_valences_[corner_table_.Vertex(prev)] += 1;
                    break;
                case EdgeBreakerTopologyBitPattern.R:
                    vertex_valences_[corner_table_.Vertex(corner)] += 1;
                    vertex_valences_[corner_table_.Vertex(next)] += 1;
                    vertex_valences_[corner_table_.Vertex(prev)] += 2;
                    break;
                case EdgeBreakerTopologyBitPattern.L:
                    vertex_valences_[corner_table_.Vertex(corner)] += 1;
                    vertex_valences_[corner_table_.Vertex(next)] += 2;
                    vertex_valences_[corner_table_.Vertex(prev)] += 1;
                    break;
                case EdgeBreakerTopologyBitPattern.E:
                    vertex_valences_[corner_table_.Vertex(corner)] += 2;
                    vertex_valences_[corner_table_.Vertex(next)] += 2;
                    vertex_valences_[corner_table_.Vertex(prev)] += 2;
                    break;
                default:
                    break;
            }
            // Compute the new predicted symbol.
            if (last_symbol_ == EdgeBreakerTopologyBitPattern.C || last_symbol_ == EdgeBreakerTopologyBitPattern.R)
            {
                int pivot =
                    corner_table_.Vertex(corner_table_.Next(corner));
                if (vertex_valences_[pivot] < 6)
                {
                    predicted_symbol_ = EdgeBreakerTopologyBitPattern.R;
                }
                else
                {
                    predicted_symbol_ = EdgeBreakerTopologyBitPattern.C;
                }
            }
            else
            {
                predicted_symbol_ = EdgeBreakerTopologyBitPattern.Invalid;
            }
        }

        public override void MergeVertices(int dest, int source)
        {
            // Update valences on the merged vertices.
            vertex_valences_[dest] += vertex_valences_[source];
        }
    }
}
