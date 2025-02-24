using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;
using Openize.Drako.Decoder;
using Openize.Drako.Encoder;

namespace Openize.Drako.Compression
{

    class PredictionSchemeNormalOctahedronTransformBase : PredictionSchemeTransform
    {

        protected OctahedronToolBox octahedronToolBox = new OctahedronToolBox();

        public int CenterValue => octahedronToolBox.CenterValue;
        public int QuantizationBits => octahedronToolBox.QuantizationBits;

        protected void SetMaxQuantizedValue(int value)
        {
            if (value % 2 == 0)
                throw new ArgumentException("Invalid quantized value");
            int q = DracoUtils.MostSignificantBit((uint)value) + 1;
            octahedronToolBox.SetQuantizationBits(q);
        }
        // For correction values.
        protected int MakePositive(int x)
        {
            return octahedronToolBox.MakePositive(x);
        }

        protected int ModMax(int x)
        {
            return octahedronToolBox.ModMax(x);
        }

        protected bool IsInDiamond(int s, int t)
        {
            return octahedronToolBox.IsInDiamond(s, t);
        }

        // We can return true as we keep correction values positive.
        public override bool AreCorrectionsPositive()
        {
            return true;
        }
    }


    /// <summary>
    /// The transform works on octahedral coordinates for normals. The square is
    /// subdivided into four inner triangles (diamond) and four outer triangles. The
    /// inner trianlges are associated with the upper part of the octahedron and the
    /// outer triangles are associated with the lower part.
    /// Given a preditiction value P and the actual value Q that should be encoded,
    /// this transform first checks if P is outside the diamond. If so, the outer
    /// triangles are flipped towards the inside and vice versa. The actual
    /// correction value is then based on the mapped P and Q values. This tends to
    /// result in shorter correction vectors.
    /// This is possible since the P value is also known by the decoder, see also
    /// ComputeCorrection and ComputeOriginalValue functions.
    /// Note that the tile is not periodic, which implies that the outer edges can
    /// not be identified, which requires us to use an odd number of values on each
    /// axis.
    /// DataTypeT is expected to be some integral type.
    ///
    /// This relates
    /// * IDF# 44535
    /// * Patent Application: GP-200957-00-US
    /// </summary>
    class PredictionSchemeNormalOctahedronTransform : PredictionSchemeNormalOctahedronTransformBase
    {

        public PredictionSchemeNormalOctahedronTransform()
        {
            
        }
        public PredictionSchemeNormalOctahedronTransform(int maxQuantizedValue)
        {
            SetMaxQuantizedValue(maxQuantizedValue);
        }

        public override PredictionSchemeTransformType Type
        {
            get { return PredictionSchemeTransformType.NormalOctahedron; }
        }


        public override void DecodeTransformData(DecoderBuffer buffer)
        {
            int maxQuantizedValue, centerValue;
            maxQuantizedValue = buffer.DecodeI32();
            if (buffer.BitstreamVersion < 22)
            {
                centerValue = buffer.DecodeI32(); 
            }

            SetMaxQuantizedValue(maxQuantizedValue);
        }

        public override void EncodeTransformData(EncoderBuffer buffer)
        {
            buffer.Encode(octahedronToolBox.MaxQuantizedValue);
        }

        public override void ComputeCorrection(Span<int> originalVals, int originalOffset, Span<int> predictedVals, int predictedOffset, Span<int> outCorrVals, int outOffset, int valId)
        {

            var orig = new IntVector(originalVals[originalOffset], originalVals[originalOffset + 1]);
            var pred = new IntVector(predictedVals[predictedOffset], predictedVals[predictedOffset + 1]);
            var corr = ComputeCorrection(orig, pred);

            outCorrVals[outOffset + valId] = corr.x;
            outCorrVals[outOffset + valId + 1] = corr.y;
        }

        public override void ComputeOriginalValue(Span<int> predictedVals, int predictedOffset, Span<int> corrVals,
            int corrOffset, Span<int> outOriginalVals,
            int outOffset)
        {
            var pred = new IntVector(predictedVals[predictedOffset + 0], predictedVals[predictedOffset + 1]);
            var corr = new IntVector(corrVals[corrOffset], corrVals[corrOffset + 1]);
            var orig = ComputeOriginalValue(pred, corr);

            outOriginalVals[outOffset + 0] = orig.x;
            outOriginalVals[outOffset + 1] = orig.y;
        }

        private IntVector ComputeCorrection(IntVector orig, IntVector pred)
        {
            var t = new IntVector(CenterValue, CenterValue);
            orig = orig - t;
            pred = pred - t;

            if (!IsInDiamond(pred.x, pred.y))
            {
                octahedronToolBox.InvertDiamond(ref orig.x, ref orig.y);
                octahedronToolBox.InvertDiamond(ref pred.x, ref pred.y);
            }

            var corr = orig - pred;
            corr.x = MakePositive(corr.x);
            corr.y = MakePositive(corr.y);
            return corr;
        }

        private IntVector ComputeOriginalValue(IntVector pred, IntVector corr)
        {
            var t = new IntVector(CenterValue, CenterValue);
            pred = pred - t;

            bool predIsInDiamond = IsInDiamond(pred.x, pred.y);
            if (!predIsInDiamond)
            {
                octahedronToolBox.InvertDiamond(ref pred);
            }
            var orig = pred + corr;
            orig.x = ModMax(orig.x);
            orig.y = ModMax(orig.y);
            if (!predIsInDiamond)
            {
                octahedronToolBox.InvertDiamond(ref orig);
            }
            orig = orig + t;
            return orig;
        }

    }
}
