using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Compression
{
    class MeshPredictionSchemeMultiParallelogram : MeshPredictionScheme
    {
        public MeshPredictionSchemeMultiParallelogram(PointAttribute attribute,
            PredictionSchemeTransform transform,
            MeshPredictionSchemeData meshData)
            :base(attribute, transform, meshData)
        {
            
        }

        public override PredictionSchemeMethod PredictionMethod { get { return PredictionSchemeMethod.MultiParallelogram;} }

        public override void ComputeCorrectionValues(Span<int> inData, Span<int> outCorr, int size, int numComponents, int[] entryToPointIdMap)
        {
            this.transform_.InitializeEncoding(inData, numComponents);
            ICornerTable table = this.meshData.CornerTable;
            var vertexToDataMap = this.meshData.vertexToDataMap;

            Span<int> predVals = stackalloc int[numComponents];

            // We start processing from the end because this prediction uses data from
            // previous entries that could be overwritten when an entry is processed.
            for (int p = this.meshData.dataToCornerMap.Count - 1; p > 0; --p)
            {
                int startCornerId = this.meshData.dataToCornerMap[p];

                // Go over all corners attached to the vertex and compute the predicted
                // value from the parallelograms defined by their opposite faces.
                int cornerId = startCornerId;
                int numParallelograms = 0;
                for (int i = 0; i < numComponents; ++i)
                {
                    predVals[i] = 0;
                }
                while (cornerId >= 0)
                {
                    // TODO(ostava): Move code shared between multiple predictors into a new
                    // file.
                    int vertOpp = p, vertNext = p, vertPrev = p;
                    int oppCorner = table.Opposite(cornerId);
                    if (oppCorner >= 0)
                    {
                        GetParallelogramEntries(oppCorner, table, vertexToDataMap, ref vertOpp, ref vertNext, ref vertPrev);
                    }
                    if (vertOpp < p && vertNext < p && vertPrev < p)
                    {
                        // Apply the parallelogram prediction.
                        int vOppOff = vertOpp * numComponents;
                        int vNextOff = vertNext * numComponents;
                        int vPrevOff = vertPrev * numComponents;
                        for (int c = 0; c < numComponents; ++c)
                        {
                            predVals[c] += (inData[vNextOff + c] + inData[vPrevOff + c]) -
                                            inData[vOppOff + c];
                        }
                        ++numParallelograms;
                    }

                    // Proceed to the next corner attached to the vertex.
                    // TODO(ostava): This will not go around the whole neighborhood on
                    // vertices on a mesh boundary. We need to SwingLeft from the start vertex
                    // again to get the full coverage.
                    cornerId = table.SwingRight(cornerId);
                    if (cornerId == startCornerId)
                    {
                        cornerId = -1;
                    }
                }
                int dstOffset = p * numComponents;
                if (numParallelograms == 0)
                {
                    // No parallelogram was valid.
                    // We use the last encoded point as a reference.
                    int srcOffset = (p - 1) * numComponents;
                    this.transform_.ComputeCorrection(inData, dstOffset, inData, srcOffset, outCorr, 0, dstOffset);
                }
                else
                {
                    // Compute the correction from the predicted value.
                    for (int c = 0; c < numComponents; ++c)
                    {
                        predVals[c] /= numParallelograms;
                    }
                    this.transform_.ComputeCorrection(inData, dstOffset, predVals, 0, outCorr, 0, dstOffset);

                }
            }
            // First element is always fixed because it cannot be predicted.
            for (int i = 0; i < numComponents; ++i)
            {
                predVals[i] = 0;
            }
            this.transform_.ComputeCorrection(inData, predVals, outCorr, 0);
        }

        public override void ComputeOriginalValues(Span<int> inCorr, Span<int> outData, int size, int numComponents, int[] entryToPointIdMap)
        {
            transform_.InitializeDecoding(numComponents);

            Span<int> predVals = stackalloc int[numComponents];
            Span<int> parallelogramPredVals = stackalloc int[numComponents];

            this.transform_.ComputeOriginalValue(predVals, inCorr, outData);

            ICornerTable table = this.meshData.CornerTable;
            var vertexToDataMap = this.meshData.vertexToDataMap;

            int cornerMapSize = this.meshData.dataToCornerMap.Count;
            for (int p = 1; p < cornerMapSize; ++p)
            {
                int startCornerId = this.meshData.dataToCornerMap[p];

                int cornerId = startCornerId;
                int numParallelograms = 0;
                for (int i = 0; i < numComponents; ++i)
                {
                    predVals[i] = 0;
                }
                while (cornerId != CornerTable.kInvalidCornerIndex) {
                  if (MeshPredictionSchemeParallelogram.ComputeParallelogramPrediction(
                          p, cornerId, table, vertexToDataMap, outData,
                          numComponents, parallelogramPredVals)) {
                    for (int c = 0; c < numComponents; ++c) {
                      predVals[c] += parallelogramPredVals[c];
                    }
                    ++numParallelograms;
                  }

                  // Proceed to the next corner attached to the vertex.
                  cornerId = table.SwingRight(cornerId);
                  if (cornerId == startCornerId) {
                    cornerId = CornerTable.kInvalidCornerIndex;
                  }
                }

                int dstOffset = p * numComponents;
                if (numParallelograms == 0)
                {
                    // No parallelogram was valid.
                    // We use the last decoded point as a reference.
                    int srcOffset = (p - 1) * numComponents;
                    this.transform_.ComputeOriginalValue(outData, srcOffset, inCorr, dstOffset,
                        outData, dstOffset);
                }
                else
                {
                    // Compute the correction from the predicted value.
                    for (int c = 0; c < numComponents; ++c)
                    {
                        predVals[c] /= numParallelograms;
                    }
                    this.transform_.ComputeOriginalValue(predVals, 0, inCorr, dstOffset,
                        outData, dstOffset);
                }
            }
        }
    }
}
