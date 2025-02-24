using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Decoder;
using Openize.Drako.Utils;

namespace Openize.Drako.Encoder
{

    /// <summary>
    /// Attribute encoder that quantizes floating point attribute values. The
    /// quantized values can be optionally compressed using an entropy coding.
    /// </summary>
    class SequentialQuantizationAttributeEncoder : SequentialIntegerAttributeEncoder
    {
        private AttributeQuantizationTransform attribute_quantization_transform_ = new AttributeQuantizationTransform();

        public override SequentialAttributeEncoderType GetUniqueId()
        {
            return SequentialAttributeEncoderType.Quantization;
        }

        public override bool IsLossyEncoder()
        {
            return true;
        }

        public override void Initialize(
            PointCloudEncoder encoder, int attributeId)
        {
            base.Initialize(encoder, attributeId);
            // This encoder currently works only for floating point attributes.
            PointAttribute attribute = Encoder.PointCloud.Attribute(attributeId);
            if (attribute.DataType != DataType.FLOAT32)
                throw DracoUtils.Failed();

            // Initialize AttributeQuantizationTransform.
            int quantization_bits = encoder.Options.GetQuantizationBits(attribute);
            if (quantization_bits < 1)
                throw DracoUtils.Failed();
            /*
            if (encoder.Options.IsAttributeOptionSet(attribute_id,
                    "quantization_origin") &&
                encoder.options().IsAttributeOptionSet(attribute_id,
                    "quantization_range"))
            {
                // Quantization settings are explicitly specified in the provided options.
                std::vector<float> quantization_origin(attribute.num_components());
                encoder.options().GetAttributeVector(attribute_id, "quantization_origin",
                    attribute.num_components(),
                    &quantization_origin[0]);
                const float range = encoder.options().GetAttributeFloat(
                    attribute_id, "quantization_range", 1.f);
                attribute_quantization_transform_.SetParameters(
                    quantization_bits, quantization_origin.data(),
                    attribute.ComponentsCount, range);
            }
            else
            */
            {
                // Compute quantization settings from the attribute values.
                attribute_quantization_transform_.ComputeParameters(attribute,
                    quantization_bits);
            }

        }

        public override void EncodeDataNeededByPortableTransform(EncoderBuffer out_buffer)
        {
            attribute_quantization_transform_.EncodeParameters(out_buffer);
        }

        protected override void PrepareValues(int[] pointIds, int numPoints)
        {
            var portable_attribute = attribute_quantization_transform_.InitTransformedAttribute(
                    attribute, pointIds.Length);
            attribute_quantization_transform_.TransformAttribute(
                    attribute, pointIds, portable_attribute);
            this.portableAttribute = portable_attribute;

        }
    }
}
