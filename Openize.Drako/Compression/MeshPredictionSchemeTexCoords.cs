using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using Openize.Drako.Decoder;
using Openize.Drako.Encoder;
using Openize.Drako.Utils;

namespace Openize.Drako.Compression
{
    class MeshPredictionSchemeTexCoords : MeshPredictionScheme
    {
        private PointAttribute posAttribute;
        private int[] entryToPointIdMap;
        private int[] predictedValue;
        private int numComponents;
        /// <summary>
        /// Encoded / decoded array of UV flips.
        /// </summary>
        private List<bool> orientations = new List<bool>();
        public MeshPredictionSchemeTexCoords(PointAttribute attribute, PredictionSchemeTransform transform_, MeshPredictionSchemeData meshData)
            :base(attribute, transform_, meshData)
        {
        }
        public override PredictionSchemeMethod PredictionMethod { get {return PredictionSchemeMethod.TexCoordsDeprecated;} }

        public override int NumParentAttributes { get { return 1; } }

        public override AttributeType GetParentAttributeType(int i)
        {
            return AttributeType.Position;
        }

        public override void SetParentAttribute(PointAttribute att)
        {
            if (att.AttributeType != AttributeType.Position)
                throw new ArgumentException(); // Invalid attribute type.
            if (att.ComponentsCount != 3)
                throw new ArgumentException(); // Currently works only for 3 component positions.
            posAttribute = att;
        }

        public override void ComputeCorrectionValues(Span<int> inData, Span<int> outCorr, int size, int numComponents, int[] entryToPointIdMap)
        {
            this.numComponents = numComponents;
            this.entryToPointIdMap = entryToPointIdMap;
            predictedValue = new int[numComponents];
            var predictedValueSpan = predictedValue.AsSpan();
            this.transform_.InitializeEncoding(inData, numComponents);
            // We start processing from the end because this prediction uses data from
            // previous entries that could be overwritten when an entry is processed.
            for (int p = this.meshData.dataToCornerMap.Count - 1; p >= 0; --p)
            {
                int cornerId = this.meshData.dataToCornerMap[p];
                ComputePredictedValue(true, cornerId, inData, p);

                int dstOffset = p * numComponents;
                this.transform_.ComputeCorrection(inData, dstOffset, predictedValueSpan, 0, outCorr, 0, dstOffset);
            }
        }

        public override void ComputeOriginalValues(Span<int> inCorr, Span<int> outData, int size, int numComponents, int[] entryToPointIdMap)
        {
            this.numComponents = numComponents;
            this.entryToPointIdMap = entryToPointIdMap;
            predictedValue = new int[numComponents];
            var predictedValueSpan = predictedValue.AsSpan();
            this.transform_.InitializeDecoding(numComponents);

            int cornerMapSize = this.meshData.dataToCornerMap.Count;
            for (int p = 0; p < cornerMapSize; ++p)
            {
                int cornerId = this.meshData.dataToCornerMap[p];
                ComputePredictedValue(false, cornerId, outData, p);

                int dstOffset = p * numComponents;
                this.transform_.ComputeOriginalValue(predictedValueSpan, 0, inCorr, dstOffset, outData, dstOffset);
            }
        }

        public override void EncodePredictionData(EncoderBuffer buffer)
        {
            // Encode the delta-coded orientations using arithmetic coding.
            int numOrientations = orientations.Count;
            buffer.Encode(numOrientations);
            bool lastOrientation = true;
            RAnsBitEncoder encoder = new RAnsBitEncoder();
            encoder.StartEncoding();
            for (int i = 0; i < orientations.Count; i++)
            {
                bool orientation = this.orientations[i];
                encoder.EncodeBit(orientation == lastOrientation);
                lastOrientation = orientation;
            }
            encoder.EndEncoding(buffer);
            base.EncodePredictionData(buffer);
        }

        public override void DecodePredictionData(DecoderBuffer buffer)
        {
            // Decode the delta coded orientations.
            int numOrientations = buffer.DecodeI32();
            orientations.Clear();
            orientations.AddRange(new bool[numOrientations]);
            bool lastOrientation = true;
            RAnsBitDecoder decoder = new RAnsBitDecoder();
            decoder.StartDecoding(buffer);
            for (int i = 0; i < numOrientations; ++i)
            {
                if (!decoder.DecodeNextBit())
                    lastOrientation = !lastOrientation;
                orientations[i] = lastOrientation;
            }
            decoder.EndDecoding();
            base.DecodePredictionData(buffer);
        }

        private Vector3 GetPositionForEntryId(int entryId)
        {
            int pointId = entryToPointIdMap[entryId];
            Vector3 pos = posAttribute.GetValueAsVector3(posAttribute.MappedIndex(pointId));
            return pos;
        }

        private Vector2 GetTexCoordForEntryId(int entryId, Span<int> data)
        {
            int dataOffset = entryId * numComponents;
            return new Vector2(data[dataOffset], data[dataOffset + 1]);
        }

        private void ComputePredictedValue(bool isEncoder, int cornerId, Span<int> data,
            int dataId)
        {
            // Compute the predicted UV coordinate from the positions on all corners
            // of the processed triangle. For the best prediction, the UV coordinates
            // on the next/previous corners need to be already encoded/decoded.
            int nextCornerId =
                this.meshData.CornerTable.Next(cornerId);
            int prevCornerId =
                this.meshData.CornerTable.Previous(cornerId);
            // Get the encoded data ids from the next and previous corners.
            // The data id is the encoding order of the UV coordinates.
            int nextDataId, prevDataId;

            int nextVertId, prevVertId;
            nextVertId =
                this.meshData.CornerTable.Vertex(nextCornerId);
            prevVertId =
                this.meshData.CornerTable.Vertex(prevCornerId);

            nextDataId = this.meshData.vertexToDataMap[nextVertId];
            prevDataId = this.meshData.vertexToDataMap[prevVertId];

            if (prevDataId < dataId && nextDataId < dataId)
            {
                // Both other corners have available UV coordinates for prediction.
                Vector2 nUv = GetTexCoordForEntryId(nextDataId, data);
                Vector2 pUv = GetTexCoordForEntryId(prevDataId, data);
                if (pUv == nUv)
                {
                    // We cannot do a reliable prediction on degenerated UV triangles.
                    predictedValue[0] = (int) (pUv.X);
                    predictedValue[1] = (int) (pUv.Y);
                    return;
                }

                // Get positions at all corners.
                Vector3 tipPos = GetPositionForEntryId(dataId);
                Vector3 nextPos = GetPositionForEntryId(nextDataId);
                Vector3 prevPos = GetPositionForEntryId(prevDataId);
                // Use the positions of the above triangle to predict the texture coordinate
                // on the tip corner C.
                // Convert the triangle into a new coordinate system defined by orthoganal
                // bases vectors S, T, where S is vector prevPos - nextPos and T is an
                // perpendicular vector to S in the same plane as vector the
                // tipPos - nextPos.
                // The transformed triangle in the new coordinate system is then going to
                // be represented as:
                //
                //        1 ^
                //          |
                //          |
                //          |   C
                //          |  /  \
                //          | /      \
                //          |/          \
                //          N--------------P
                //          0              1
                //
                // Where nextPos point (N) is at position (0, 0), prevPos point (P) is
                // at (1, 0). Our goal is to compute the position of the tipPos point (C)
                // in this new coordinate space (s, t).
                //
                Vector3 pn = prevPos - nextPos;
                Vector3 cn = tipPos - nextPos;
                var pnNorm2Squared = Vector3.Dot(pn, pn);
                // Coordinate s of the tip corner C is simply the dot product of the
                // normalized vectors |pn| and |cn| (normalized by the length of |pn|).
                // Since both of these vectors are normalized, we don't need to perform the
                // normalization explicitly and instead we can just use the squared norm
                // of |pn| as a denominator of the resulting dot product of non normalized
                // vectors.
                var s = Vector3.Dot(pn, cn) / pnNorm2Squared;
                // To get the coordinate t, we can use formula:
                //      t = |C-N - (P-N) * s| / |P-N|
                var t = (float) Math.Sqrt((cn - pn * s).LengthSquared() / pnNorm2Squared);

                // Now we need to transform_ the point (s, t) to the texture coordinate space
                // UV. We know the UV coordinates on points N and P (NUV and PUV). Lets
                // denote PUV - NUV = PNUV. PNUV is then 2 dimensional vector that can
                // be used to define transformation from the normalized coordinate system
                // to the texture coordinate system using a 3x3 affine matrix M:
                //
                //  M = | PNUV[0]  -PNUV[1]  NUV[0] |
                //      | PNUV[1]   PNUV[0]  NUV[1] |
                //      | 0          0         1       |
                //
                // The predicted point CUV in the texture space is then equal to
                // CUV = M * (s, t, 1). Because the triangle in UV space may be flipped
                // around the PNUV axis, we also need to consider point CUV' = M * (s, -t)
                // as the prediction.
                Vector2 pnUv = pUv - nUv;
                var pnus = pnUv.X * s + nUv.X;
                var pnut = pnUv.X * t;
                var pnvs = pnUv.Y * s + nUv.Y;
                var pnvt = pnUv.Y * t;
                Vector2 predictedUv;
                if (isEncoder)
                {
                    // When encoding compute both possible vectors and determine which one
                    // results in a better prediction.
                    Vector2 predictedUv0 = new Vector2(pnus - pnvt, pnvs + pnut);
                    Vector2 predictedUv1 = new Vector2(pnus + pnvt, pnvs - pnut);
                    Vector2 cUv = GetTexCoordForEntryId(dataId, data);
                    if ((cUv - predictedUv0).LengthSquared() <
                        (cUv - predictedUv1).LengthSquared())
                    {
                        predictedUv = predictedUv0;
                        orientations.Add(true);
                    }
                    else
                    {
                        predictedUv = predictedUv1;
                        orientations.Add(false);
                    }
                }
                else
                {
                    // When decoding the data, we already know which orientation to use.
                    bool orientation = orientations[orientations.Count - 1];
                    orientations.RemoveAt(orientations.Count - 1);
                    if (orientation)
                        predictedUv = new Vector2(pnus - pnvt, pnvs + pnut);
                    else
                        predictedUv = new Vector2(pnus + pnvt, pnvs - pnut);
                }
                // Round the predicted value for integer types.
                predictedValue[0] = (int) (Math.Floor(predictedUv.X + 0.5f));
                predictedValue[1] = (int) (Math.Floor(predictedUv.Y + 0.5f));
                return;
            }
            // Else we don't have available textures on both corners. For such case we
            // can't use positions for predicting the uv value and we resort to delta
            // coding.
            int dataOffset = 0;
            if (prevDataId < dataId)
            {
                // Use the value on the previous corner as the prediction.
                dataOffset = prevDataId * numComponents;
            }
            if (nextDataId < dataId)
            {
                // Use the value on the next corner as the prediction.
                dataOffset = nextDataId * numComponents;
            }
            else
            {
                // None of the other corners have a valid value. Use the last encoded value
                // as the prediction if possible.
                if (dataId > 0)
                {
                    dataOffset = (dataId - 1) * numComponents;
                }
                else
                {
                    // We are encoding the first value. Predict 0.
                    for (int i = 0; i < numComponents; ++i)
                    {
                        predictedValue[i] = 0;
                    }
                    return;
                }
            }
            for (int i = 0; i < numComponents; ++i)
            {
                predictedValue[i] = data[dataOffset + i];
            }
        }
    }
}
