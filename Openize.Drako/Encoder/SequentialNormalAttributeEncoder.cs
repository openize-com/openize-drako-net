using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using Openize.Drako.Compression;
using Openize.Drako.Utils;

namespace Openize.Drako.Encoder
{

    /// <summary>
    /// Class for encoding normal vectors using an octahedral encoding (see Cigolle
    /// et al.'14 “A Survey of Efficient Representations for Independent Unit
    /// Vectors”. Compared to the basic quantization encoder, this encoder results
    /// in a better compression rate under the same accuracy settings. Note that this
    /// encoder doesn't preserve the lengths of input vectors, therefore it will not
    /// work correctly when the input values are not normalized.
    /// </summary>
    class SequentialNormalAttributeEncoder : SequentialIntegerAttributeEncoder
    {
        private AttributeOctahedronTransform attribute_octahedron_transform_;
        public override SequentialAttributeEncoderType GetUniqueId()
        {
            return SequentialAttributeEncoderType.Normals;
        }

        public override bool IsLossyEncoder()
        {
            return true;
        }

        public override void Initialize(PointCloudEncoder encoder, int attributeId)
        {
            base.Initialize(encoder, attributeId);
            // Currently this encoder works only for 3-component normal vectors.
            if (Attribute.ComponentsCount != 3)
                throw DracoUtils.Failed();
            var q = encoder.Options.GetQuantizationBits(attribute);
            if (q < 1)
                throw DracoUtils.Failed();
            attribute_octahedron_transform_ = new AttributeOctahedronTransform(q);
        }
        public override void EncodeDataNeededByPortableTransform(EncoderBuffer out_buffer)
        {
            attribute_octahedron_transform_.EncodeParameters(out_buffer);
        }

        protected override void PrepareValues(int[] pointIds, int numPoints)
        {
            this.portableAttribute = attribute_octahedron_transform_.GeneratePortableAttribute(Attribute, pointIds, numPoints);
        }

// Converts a unit vector into octahedral coordinates (0-1 range).
        static void UnitVectorToOctahedralCoords(ref Vector3 vector, out float outS, out float outT)
        {
            float absSum = (float)(Math.Abs(vector.X) + Math.Abs(vector.Y) + Math.Abs(vector.Z));
            Vector3 scaledVec = new Vector3();
            if (absSum > 1e-6)
            {
                // Scale needed to project the vector to the surface of an octahedron.
                float scale = 1.0f / absSum;
                scaledVec.X = vector.X * scale;
                scaledVec.Y = vector.Y * scale;
                scaledVec.Z = vector.Z * scale;
            }
            else
            {
                scaledVec.X = 1;
                scaledVec.Y = 0;
                scaledVec.Z = 0;
            }

            if (scaledVec.X >= 0.0f)
            {
                // Right hemisphere.
                outS = (float)((scaledVec.Y + 1.0f) * 0.5f);
                outT = (float)((scaledVec.Z + 1.0f) * 0.5f);
            }
            else
            {
                // Left hemisphere.
                if (scaledVec.Y < 0.0f)
                {
                    outS = (float)(0.5f * Math.Abs(scaledVec.Z));
                }
                else
                {
                    outS = (float)(0.5f * (2.0f - Math.Abs(scaledVec.Z)));
                }
                if (scaledVec.Z < 0.0f)
                {
                    outT = (float)(0.5f * Math.Abs(scaledVec.Y));
                }
                else
                {
                    outT = (float)(0.5f * (2.0f - Math.Abs(scaledVec.Y)));
                }
            }
        }

        static void UnitVectorToQuantizedOctahedralCoords(ref Vector3 vector,
            float maxQuantizedValue,
            out int outS, out int outT)
        {
            // In order to be able to represent the center normal we reduce the range
            // by one.
            float maxValue = maxQuantizedValue - 1;
            float ss, tt;
            UnitVectorToOctahedralCoords(ref vector, out ss, out tt);
            int s = (int) (Math.Floor(ss * maxValue + 0.5));
            int t = (int) (Math.Floor(tt * maxValue + 0.5));

            int centerValue = (int) (maxValue / 2);

            // Convert all edge points in the top left and bottom right quadrants to
            // their corresponding position in the bottom left and top right quadrants.
            // Convert all corner edge points to the top right corner. This is necessary
            // for the inversion to occur correctly.
            if ((s == 0 && t == 0) || (s == 0 && t == maxValue) ||
                (s == maxValue && t == 0))
            {
                s = (int) maxValue;
                t = (int) maxValue;
            }
            else if (s == 0 && t > centerValue)
            {
                t = centerValue - (t - centerValue);
            }
            else if (s == maxValue && t < centerValue)
            {
                t = centerValue + (centerValue - t);
            }
            else if (t == maxValue && s < centerValue)
            {
                s = centerValue + (centerValue - s);
            }
            else if (t == 0 && s > centerValue)
            {
                s = centerValue - (s - centerValue);
            }

            outS = s;
            outT = t;
        }

        protected override PredictionScheme CreateIntPredictionScheme(PredictionSchemeMethod method)
        {
            int quantizationBits = Encoder.Options.GetQuantizationBits(Attribute);
            int maxValue = (1 << quantizationBits) - 1;
            var transform = new PredictionSchemeNormalOctahedronCanonicalizedTransform(maxValue);

            PredictionSchemeMethod prediction_method = SelectPredictionMethod(AttributeId, Encoder);
            if(prediction_method == PredictionSchemeMethod.GeometricNormal ||  prediction_method == PredictionSchemeMethod.Difference)
                return PredictionScheme.Create(Encoder, prediction_method, AttributeId, transform);
            return null;
        }

    }
}
