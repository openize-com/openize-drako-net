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
    /// Base class containing shared functionality used by both encoding and decoding
    /// canonicalized normal octahedron prediction scheme transforms. See the
    /// encoding transform for more details about the method.
    /// </summary>
    class PredictionSchemeNormalOctahedronCanonicalizedTransformBase : PredictionSchemeNormalOctahedronTransformBase
    {
        public override PredictionSchemeTransformType Type => PredictionSchemeTransformType.NormalOctahedronCanonicalized;

        public PredictionSchemeNormalOctahedronCanonicalizedTransformBase()
        {

        }

        public PredictionSchemeNormalOctahedronCanonicalizedTransformBase(int max_quantized_value)
        {
            this.SetMaxQuantizedValue(max_quantized_value);
        }

        protected int GetRotationCount(IntVector pred)
        {
            return GetRotationCount(pred.x, pred.y);
        }
        protected int GetRotationCount(int sign_x, int sign_y)
        {

            int rotation_count = 0;
            if (sign_x == 0)
            {
                if (sign_y == 0)
                {
                    rotation_count = 0;
                }
                else if (sign_y > 0)
                {
                    rotation_count = 3;
                }
                else
                {
                    rotation_count = 1;
                }
            }
            else if (sign_x > 0)
            {
                if (sign_y >= 0)
                {
                    rotation_count = 2;
                }
                else
                {
                    rotation_count = 1;
                }
            }
            else
            {
                if (sign_y <= 0)
                {
                    rotation_count = 0;
                }
                else
                {
                    rotation_count = 3;
                }
            }

            return rotation_count;
        }

        protected void RotatePoint(ref IntVector p, int rotation_count)
        {
            int s = p.x;
            int t = p.y;
            switch (rotation_count)
            {
                case 1:
                    p.x = t;
                    p.y = -s;
                    break;// return new int[]{p[1], -p[0]};
                case 2:
                    p.x = -s;
                    p.y = -t;
                    break;// return new int[]{-p[0], -p[1]};
                case 3:
                    p.x = -t;
                    p.y = s;
                    break;// return new int[]{-p[1], p[0]};
                default:
                    break;//return p;
            }
        }

        protected bool IsInBottomLeft(IntVector p)
        {
            return IsInBottomLeft(p.x, p.y);
        }

        protected bool IsInBottomLeft(int s, int t)
        {
            if (s == 0 && t == 0)
                return true;
            return (s < 0 && t <= 0);
        }
    }
    class PredictionSchemeNormalOctahedronCanonicalizedTransform : PredictionSchemeNormalOctahedronCanonicalizedTransformBase
    {


        public PredictionSchemeNormalOctahedronCanonicalizedTransform()
        {

        }

        public PredictionSchemeNormalOctahedronCanonicalizedTransform(int mod_value)
            :base(mod_value)
        {

        }

        public override void InitializeDecoding(int numComponents)
        {
        }


        public override void DecodeTransformData(DecoderBuffer buffer)
        {
            int max_quantized_value = buffer.DecodeI32();
            int center_value = buffer.DecodeI32();

            SetMaxQuantizedValue(max_quantized_value);
            // Account for reading wrong values, e.g., due to fuzzing.
            if (QuantizationBits < 2)
                throw DracoUtils.Failed();
            if (QuantizationBits > 30)
                throw DracoUtils.Failed();
        }

        public override void ComputeOriginalValue(Span<int> predictedVals, int predictedOffset, Span<int> corrVals,
            int corrOffset, Span<int> outOriginalVals, int outOffset)
        {
            var pred = new IntVector(predictedVals[predictedOffset + 0], predictedVals[predictedOffset + 1]);
            IntVector corr = new IntVector(corrVals[corrOffset + 0], corrVals[corrOffset + 1]);
            IntVector orig = ComputeOriginalValue(pred, corr);

            outOriginalVals[outOffset] = orig.x;
            outOriginalVals[outOffset + 1] = orig.y;
        }

        private IntVector ComputeOriginalValue(IntVector pred, IntVector corr)
        {
            IntVector t = new IntVector(CenterValue, this.CenterValue);
            pred.x -= t.x;
            pred.y -= t.y;
            bool pred_is_in_diamond = IsInDiamond(pred.x, pred.y);
            if (!pred_is_in_diamond)
            {
                octahedronToolBox.InvertDiamond(ref pred);
            }

            bool pred_is_in_bottom_left = this.IsInBottomLeft(pred);
            int rotation_count = this.GetRotationCount(pred);
            if (!pred_is_in_bottom_left)
            {
                this.RotatePoint(ref pred, rotation_count);
            }

            var orig = new IntVector();
            orig.x = this.ModMax(pred.x + corr.x);
            orig.y = this.ModMax(pred.y + corr.y);
            if (!pred_is_in_bottom_left)
            {
                int reverse_rotation_count = (4 - rotation_count) % 4;
                this.RotatePoint(ref orig, reverse_rotation_count);
            }

            if (!pred_is_in_diamond)
            {
                octahedronToolBox.InvertDiamond(ref orig);
            }

            orig.x += t.x;
            orig.y += t.y;
            return orig;
        }

        public override void EncodeTransformData(EncoderBuffer buffer)
        {
            buffer.Encode(this.octahedronToolBox.MaxQuantizedValue);
            buffer.Encode(this.octahedronToolBox.CenterValue);
        }
        public override void ComputeCorrection(Span<int> originalVals, int originalOffset, Span<int> predictedVals,
            int predictedOffset, Span<int> outCorrVals, int outOffset, int valId)
        {

            var orig = new IntVector(originalVals[0] - CenterValue, originalVals[1] - CenterValue);
            var pred = new IntVector(predictedVals[0] - CenterValue, predictedVals[1] - CenterValue);
            if (!this.IsInDiamond(pred.x, pred.y))
            {
                octahedronToolBox.InvertDiamond(ref orig);
                octahedronToolBox.InvertDiamond(ref pred);
            }

            if (!this.IsInBottomLeft(pred))
            {
                int rotation_count = this.GetRotationCount(pred);
                this.RotatePoint(ref orig, rotation_count);
                this.RotatePoint(ref pred, rotation_count);
            }

            outCorrVals[0] = this.MakePositive(orig.x - pred.x);
            outCorrVals[1] = this.MakePositive(orig.y - pred.y);
        }
    }
}
