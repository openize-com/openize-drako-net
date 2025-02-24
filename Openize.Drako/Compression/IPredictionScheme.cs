using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Decoder;
using Openize.Drako.Encoder;
using Openize.Drako.Utils;

namespace Openize.Drako.Compression
{
    abstract partial class PredictionScheme
    {
        protected PointAttribute attribute_;
        protected PredictionSchemeTransform transform_;

        protected PredictionScheme(PointAttribute attribute, PredictionSchemeTransform transform)
        {
            this.attribute_ = attribute;
            this.transform_ = transform;
        }


        /// <summary>
        /// Returns the encoded attribute.
        /// </summary>
        public PointAttribute Attribute => attribute_;

        public abstract bool Initialized { get; }

        /// <summary>
        /// Returns the number of parent attributes that are needed for the prediction.
        /// </summary>
        public virtual int NumParentAttributes
        {
            get { return 0; }
        }

        public virtual AttributeType GetParentAttributeType(int i)
        {
            return AttributeType.Invalid;
        }


      // Sets the required parent attribute.
        public virtual void SetParentAttribute(PointAttribute att)
        {
            throw DracoUtils.Failed();
        }

        public virtual bool AreCorrectionsPositive()
        {
            return transform_.AreCorrectionsPositive();
        }

        public virtual PredictionSchemeTransformType TransformType
        {
            get { return transform_.Type; }
        }
        public abstract PredictionSchemeMethod PredictionMethod { get; }


        public virtual void EncodePredictionData(EncoderBuffer buffer)
        {
            transform_.EncodeTransformData(buffer);
        }

        // Method that can be used to decode any prediction scheme specific data
        // from the input buffer.
        public virtual void DecodePredictionData(DecoderBuffer buffer)
        {
            transform_.DecodeTransformData(buffer);
        }

        public abstract void ComputeCorrectionValues(Span<int> in_data, Span<int> out_corr, int size, int num_components,
            int[] entry_to_point_id_map);

        // Reverts changes made by the prediction scheme during encoding.
        public abstract void ComputeOriginalValues(Span<int> in_corr, Span<int> out_data, int size, int num_components,
            int[] entry_to_point_id_map);

    }
}
