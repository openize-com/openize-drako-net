using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Decoder;
using FileFormat.Drako.Encoder;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Compression
{

    /// <summary>
    /// PredictionSchemeWrapTransform uses the min and max bounds of the original
    /// data to wrap stored correction values around these bounds centered at 0,
    /// i.e., when the range of the original values O is between &lt;MIN, MAX&gt; and
    /// N = MAX-MIN, we can then store any correction X = O - P, as:
    ///        X + N,   if X &lt; -N / 2
    ///        X - N,   if X &gt; N / 2
    ///        X        otherwise
    /// To unwrap this value, the decoder then simply checks whether the final
    /// corrected value F = P + X is out of the bounds of the input values.
    /// All out of bounds values are unwrapped using
    ///        F + N,   if F &lt; MIN
    ///        F - N,   if F &gt; MAX
    /// This wrapping can reduce the number of unique values, which translates to a
    /// better entropy of the stored values and better compression rates.
    /// </summary>
    class PredictionSchemeWrapTransform : PredictionSchemeTransform
    {
        private int minValue;
        private int maxValue;
        private int maxDif;
        private int maxCorrection;
        private int minCorrection;
        private int[] clampedValue = null;

        protected void InitCorrectionBounds()
        {
            maxDif = 1 + maxValue - minValue;
            maxCorrection = maxDif / 2;
            minCorrection = -maxCorrection;
            if ((maxDif & 1) == 0)
                maxCorrection -= 1;
        }

        public override PredictionSchemeTransformType Type
        {

            get { return PredictionSchemeTransformType.Wrap; }
        }

        public override void InitializeEncoding(Span<int> origData, int numComponents)
        {
            base.InitializeEncoding(origData, numComponents);
            // Go over the original values and compute the bounds.
            if (origData.Length == 0)
                return;
            minValue = maxValue = origData[0];
            for (int i = 1; i < origData.Length; ++i)
            {
                if (origData[i] < minValue)
                    minValue = origData[i];
                else if (origData[i] > maxValue)
                    maxValue = origData[i];
            }
            InitCorrectionBounds();
            Array.Resize(ref clampedValue, numComponents);
        }

        public override void InitializeDecoding(int numComponents)
        {
            base.InitializeDecoding(numComponents);
            Array.Resize(ref clampedValue, numComponents);
        }

        /// <summary>
        /// Computes the corrections based on the input original value and the
        /// predicted value. Out of bound correction values are wrapped around the max
        /// range of input values.
        /// </summary>
        public override void ComputeCorrection(Span<int> originalVals, int originalOffset,
            Span<int> predictedVals, int predictedOffset,
            Span<int> outCorrVals, int outOffset, int valId)
        {
            base.ComputeCorrection(originalVals, originalOffset, ClampPredictedValue(predictedVals, predictedOffset), 0,
                outCorrVals,
                outOffset, valId);
            // Wrap around if needed.
            for (int i = 0; i < numComponents; ++i)
            {
                int idx = outOffset + valId + i;
                int corrVal = outCorrVals[idx];
                if (corrVal < minCorrection)
                    outCorrVals[idx] = corrVal + maxDif;
                else if (corrVal > maxCorrection)
                    outCorrVals[idx] = corrVal - maxDif;
            }

            //////////////////
            /*

            for (int i = 0; i < numComponents; ++i)
            {
                ClampPredictedValue(predictedVals, predictedOffset);
                outCorrVals[i + outOffset] = originalVals[i + originalOffset] - predictedVals[i + predictedOffset];
                // Wrap around if needed.
                var corrVal = outCorrVals[i + originalOffset];
                if (corrVal < minCorrection)
                    corrVal += maxDif;
                else if (corrVal > maxCorrection)
                    corrVal -= maxDif;
                outCorrVals[i + originalOffset] = corrVal;
            }
            */

        }

        /// <summary>
        /// Computes the original value from the input predicted value and the decoded
        /// corrections. Values out of the bounds of the input values are unwrapped.
        /// </summary>
        public override void ComputeOriginalValue(Span<int> predictedVals,
            int predictedOffset,
            Span<int> corrVals,
            int corrOffset,
            Span<int> outOriginalVals, int outOffset)
        {
            //base.ComputeOriginalValue(ClampPredictedValue(predictedVals, predictedOffset), 0, corrVals, corrOffset, outOriginalVals, outOffset, valId);
            predictedVals = ClampPredictedValue(predictedVals, predictedOffset);
            for (int i = 0; i < numComponents; ++i)
            {
                int n = i + outOffset;
                outOriginalVals[n] = predictedVals[i] + corrVals[i + corrOffset];
                if (outOriginalVals[n] > maxValue)
                    outOriginalVals[n] -= maxDif;
                else if (outOriginalVals[n] < minValue)
                    outOriginalVals[n] += maxDif;
            }
        }

        Span<int> ClampPredictedValue(Span<int> predictedVal, int offset)
        {
            for (int i = 0; i < numComponents; ++i)
            {
                int v = predictedVal[i + offset];
                if (v > maxValue)
                    clampedValue[i] = maxValue;
                else if (v < minValue)
                    clampedValue[i] = minValue;
                else
                    clampedValue[i] = v;
            }
            return clampedValue.AsSpan();
        }

        public override void EncodeTransformData(EncoderBuffer buffer)
        {
            // Store the input value range as it is needed by the decoder.
            buffer.Encode(minValue);
            buffer.Encode(maxValue);
        }

        public override void DecodeTransformData(DecoderBuffer buffer)
        {
            minValue = buffer.DecodeI32();
            maxValue = buffer.DecodeI32();
            InitCorrectionBounds();
        }
    }
}
