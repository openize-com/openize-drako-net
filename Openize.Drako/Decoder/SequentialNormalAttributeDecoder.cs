using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using FileFormat.Drako.Compression;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Decoder
{
    class SequentialNormalAttributeDecoder : SequentialIntegerAttributeDecoder
    {
        private int quantizationBits = -1;


        protected override int GetNumValueComponents()
        {
            return 2;
        }

        public override void Initialize(PointCloudDecoder decoder, int attributeId)
        {
            base.Initialize(decoder, attributeId);
            // Currently, this encoder works only for 3-component normal vectors.
            if (Attribute.ComponentsCount != 3)
                throw DracoUtils.Failed();
            // Also the data type must be DTFLOAT32.
            if (Attribute.DataType != DataType.FLOAT32)
                throw DracoUtils.Failed();
        }

        public override void DecodeIntegerValues(int[] pointIds, DecoderBuffer inBuffer)
        {
            byte quantizationBits;
            if (decoder.BitstreamVersion < 20)
            {
                quantizationBits = inBuffer.DecodeU8();
                this.quantizationBits = quantizationBits;
            }
            base.DecodeIntegerValues(pointIds, inBuffer);
        }

        public override void DecodeDataNeededByPortableTransform(int[] pointIds, DecoderBuffer in_buffer)
        {
            if (decoder.BitstreamVersion >= 20)
            {

                // For newer file version, decode attribute transform data here.
                byte quantization_bits = in_buffer.DecodeU8();
                quantizationBits = quantization_bits;
            }

            // Store the decoded transform data in portable attribute.
            AttributeOctahedronTransform octahedral_transform = new AttributeOctahedronTransform(quantizationBits);
            octahedral_transform.TransferToAttribute(PortableAttribute);
        }

        protected override void StoreValues(int numPoints)
        {
            // Convert all quantized values back to floats.
            int maxQuantizedValue = (1 << quantizationBits) - 1;
            float maxQuantizedValueF = (float) (maxQuantizedValue);

            int numComponents = Attribute.ComponentsCount;
            int entrySize = sizeof(float) * numComponents;
            float[] attVal = new float[3];
            int quantValId = 0;
            int outBytePos = 0;
            var values = MemoryMarshal.Cast<byte, int>(PortableAttribute.Buffer.GetBuffer().AsSpan(0, numPoints * 2 * 4));
            for (int i = 0; i < numPoints; ++i)
            {
                int s = values[quantValId++];
                int t = values[quantValId++];
                QuantizedOctaherdalCoordsToUnitVector(s, t, maxQuantizedValueF, attVal);
                // Store the decoded floating point value into the attribute buffer.
                Attribute.Buffer.Write(outBytePos, attVal);
                outBytePos += entrySize;
            }
        }

        void OctaherdalCoordsToUnitVector(float inS, float inT, float[] outVector)
        {
            float s = inS;
            float t = inT;
            float spt = s + t;
            float smt = s - t;
            float xSign = 1.0f;
            if (spt >= 0.5f && spt <= 1.5f && smt >= -0.5f && smt <= 0.5f)
            {
                // Right hemisphere. Don't do anything.
            }
            else
            {
                // Left hemisphere.
                xSign = -1.0f;
                if (spt <= 0.5)
                {
                    s = 0.5f - inT;
                    t = 0.5f - inS;
                }
                else if (spt >= 1.5f)
                {
                    s = 1.5f - inT;
                    t = 1.5f - inS;
                }
                else if (smt <= -0.5f)
                {
                    s = inT - 0.5f;
                    t = inS + 0.5f;
                }
                else
                {
                    s = inT + 0.5f;
                    t = inS - 0.5f;
                }
                spt = s + t;
                smt = s - t;
            }
            float y = 2.0f * s - 1.0f;
            float z = 2.0f * t - 1.0f;
            float x = Math.Min(Math.Min(2.0f * spt - 1.0f, 3.0f - 2.0f * spt),
                          Math.Min(2.0f * smt + 1.0f, 1.0f - 2.0f * smt)) *
                      xSign;
            // Normalize the computed vector.
            float normSquared = x * x + y * y + z * z;
            if (normSquared < 1e-6)
            {
                outVector[0] = 0;
                outVector[1] = 0;
                outVector[2] = 0;
            }
            else
            {
                float d = 1.0f / (float)Math.Sqrt(normSquared);
                outVector[0] = x * d;
                outVector[1] = y * d;
                outVector[2] = z * d;
            }
        }

        void QuantizedOctaherdalCoordsToUnitVector(int inS, int inT,
            float maxQuantizedValue,
            float[] outVector)
        {
            // In order to be able to represent the center normal we reduce the range
            // by one. Also note that we can not simply identify the lower left and the
            // upper right edge of the tile, which forces us to use one value less.
            maxQuantizedValue -= 1;
            OctaherdalCoordsToUnitVector(inS / maxQuantizedValue,
                inT / maxQuantizedValue, outVector);
        }

        protected override PredictionScheme CreateIntPredictionScheme(PredictionSchemeMethod method,
            PredictionSchemeTransformType transformType)
        {

            // At this point the decoder has not read the quantization bits,
            // which is why we must construct the transform by default.
            // See Transform.DecodeTransformData for more details.
            if (transformType == PredictionSchemeTransformType.NormalOctahedron)
                return PredictionScheme.Create((MeshDecoder)Decoder, method, AttributeId,
                    new PredictionSchemeNormalOctahedronTransform());
            if (transformType == PredictionSchemeTransformType.NormalOctahedronCanonicalized)
            {
                return PredictionScheme.Create(Decoder, method, AttributeId,
                    new PredictionSchemeNormalOctahedronCanonicalizedTransform());
            }

            return null; // Currently, we support only octahedron transform and
            // octahedron transform canonicalized.
        }
    }
}
