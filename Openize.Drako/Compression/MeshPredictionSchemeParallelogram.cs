using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Compression
{
    /// <summary>
    /// Parallelogram prediction predicts an attribute value V from three vertices
    /// on the opposite face to the predicted vertex. The values on the three
    /// vertices are used to construct a parallelogram V' = O - A - B, where O is the
    /// value on the oppoiste vertex, and A, B are values on the shared vertices:
    ///     V
    ///    / \
    ///   /   \
    ///  /     \
    /// A-------B
    ///  \     /
    ///   \   /
    ///    \ /
    ///     O
    /// </summary>
    class MeshPredictionSchemeParallelogram : MeshPredictionScheme
    {

        public MeshPredictionSchemeParallelogram(PointAttribute attribute,
            PredictionSchemeTransform transform_,
            MeshPredictionSchemeData meshData)
            :base(attribute, transform_, meshData)
        {
            
        }
        public override PredictionSchemeMethod PredictionMethod { get {return PredictionSchemeMethod.Parallelogram;} }

        public override void ComputeCorrectionValues(Span<int> inData, Span<int> outCorr, int size, int numComponents, int[] entryToPointIdMap)
        {
            transform_.InitializeEncoding(inData, numComponents);
            Span<int> predVals = stackalloc int[numComponents];

            // We start processing from the end because this prediction uses data from
            // previous entries that could be overwritten when an entry is processed.
            ICornerTable table = this.meshData.CornerTable;
            var vertexToDataMap = meshData.vertexToDataMap;
            for (int p = meshData.dataToCornerMap.Count - 1; p > 0; --p)
            {
                int cornerId = meshData.dataToCornerMap[p];

                int dst_offset = p * numComponents;
                if (!ComputeParallelogramPrediction(p, cornerId, table,
                                                    vertexToDataMap, inData,
                                                    numComponents, predVals))
                {
                    // Parallelogram could not be computed, Possible because some of the
                    // vertices are not valid (not encoded yet).
                    // We use the last encoded point as a reference (delta coding).
                    int src_offset = (p - 1) * numComponents;
                    this.transform_.ComputeCorrection(
                        inData, dst_offset, inData, src_offset, outCorr, dst_offset, 0);
                }
                else
                {
                    // Apply the parallelogram prediction.
                    this.transform_.ComputeCorrection(inData, dst_offset, predVals, 0,
                                                        outCorr, dst_offset, 0);
                }
                /*
                // Initialize the vertex ids to "p" which ensures that if the opposite
                // corner does not exist we will not use the vertices to predict the
                // encoded value.
                int vertOpp = p, vertNext = p, vertPrev = p;
                int oppCorner = table.Opposite(cornerId);
                if (oppCorner >= 0)
                {
                    // Get vertices on the opposite face.
                    GetParallelogramEntries(oppCorner, table, vertexToDataMap, ref vertOpp, ref vertNext,
                        ref vertPrev);
                }
                int dstOffset = p * numComponents;
                if (vertOpp >= p || vertNext >= p || vertPrev >= p)
                {
                    // Some of the vertices are not valid (not encoded yet).
                    // We use the last encoded point as a reference.
                    int srcOffset = (p - 1) * numComponents;
                    transform_.ComputeCorrection(inData, dstOffset, inData, srcOffset, outCorr, dstOffset, 0);
                }
                else
                {
                    // Apply the parallelogram prediction.
                    int vOppOff = vertOpp * numComponents;
                    int vNextOff = vertNext * numComponents;
                    int vPrevOff = vertPrev * numComponents;
                    for (int c = 0; c < numComponents; ++c)
                    {
                        predVals[c] = (inData[vNextOff + c] + inData[vPrevOff + c]) -
                                       inData[vOppOff + c];
                    }
                    transform_.ComputeCorrection(inData, dstOffset, predVals, 0,
                        outCorr, dstOffset, 0);
                }
                */
            }
            // First element is always fixed because it cannot be predicted.
            for (int i = 0; i < numComponents; ++i)
            {
                predVals[i] = 0;
            }
            transform_.ComputeCorrection(inData, predVals, outCorr, 0);
        }

        public override void ComputeOriginalValues(Span<int> inCorr, Span<int> outData, int size, int numComponents, int[] entryToPointIdMap)
        {
            transform_.InitializeDecoding(numComponents);

            ICornerTable table = this.meshData.CornerTable;
            var vertexToDataMap = this.meshData.vertexToDataMap;

            Span<int> predVals = stackalloc int[numComponents];

            // Restore the first value.
            this.transform_.ComputeOriginalValue(predVals, inCorr, outData);

            int cornerMapSize = this.meshData.dataToCornerMap.Count;
            for (int p = 1; p < cornerMapSize; ++p)
            {
                int corner_id = this.meshData.dataToCornerMap[p];
                int dst_offset = p * numComponents;
                if (!ComputeParallelogramPrediction(p, corner_id, table,
                    vertexToDataMap, outData,
                    numComponents, predVals))
                {
                    // Parallelogram could not be computed, Possible because some of the
                    // vertices are not valid (not encoded yet).
                    // We use the last encoded point as a reference (delta coding).
                    int src_offset = (p - 1) * numComponents;
                    this.transform_.ComputeOriginalValue(
                        outData, src_offset, inCorr, dst_offset, outData, dst_offset);
                }
                else
                {
                    // Apply the parallelogram prediction.
                    this.transform_.ComputeOriginalValue(
                        predVals, 0, inCorr, dst_offset, outData, dst_offset);
                }
            }
        }

// Computes parallelogram prediction for a given corner and data entry id.
// The prediction is stored in |out_prediction|.
// Function returns false when the prediction couldn't be computed, e.g. because
// not all entry points were available.
        public static bool ComputeParallelogramPrediction(
            int data_entry_id, int ci, ICornerTable table,
            int[] vertex_to_data_map, Span<int> in_data,
            int num_components, Span<int> out_prediction)
        {
            int oci = table.Opposite(ci);
            if (oci == CornerTable.kInvalidCornerIndex)
                return false;
            int vert_opp = 0, vert_next = 0, vert_prev = 0;
            GetParallelogramEntries(oci, table, vertex_to_data_map,
                ref vert_opp, ref vert_next, ref vert_prev);
            if (vert_opp < data_entry_id && vert_next < data_entry_id &&
                vert_prev < data_entry_id)
            {
                // Apply the parallelogram prediction.
                int v_opp_off = vert_opp * num_components;
                int v_next_off = vert_next * num_components;
                int v_prev_off = vert_prev * num_components;
                for (int c = 0; c < num_components; ++c)
                {
                    out_prediction[c] = (in_data[v_next_off + c] + in_data[v_prev_off + c]) -
                                        in_data[v_opp_off + c];
                }

                return true;
            }

            return false; // Not all data is available for prediction
        }
    }
}
