using FileFormat.Drako.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Decoder
{
    class SequentialAttributeDecodersController : AttributesDecoder
    {
        private SequentialAttributeDecoder[] sequentialDecoders;
        private int[] pointIds;
        private PointsSequencer sequencer;
        public SequentialAttributeDecodersController(PointsSequencer sequencer)
        {
            this.sequencer = sequencer;
        }

        public override void DecodeAttributesDecoderData(DecoderBuffer buffer)
        {

            base.DecodeAttributesDecoderData(buffer);
            // Decode unique ids of all sequential encoders and create them.
            sequentialDecoders = new SequentialAttributeDecoder[NumAttributes];
            for (int i = 0; i < NumAttributes; ++i)
            {
                byte decoderType = buffer.DecodeU8();
                // Create the decoder from the id.
                sequentialDecoders[i] = CreateSequentialDecoder((SequentialAttributeEncoderType)decoderType);
                if (sequentialDecoders[i] == null)
                    throw DracoUtils.Failed();
                sequentialDecoders[i].Initialize(Decoder, GetAttributeId(i));
            }
        }

        public override void DecodeAttributes(DecoderBuffer buffer)
        {

            if (sequencer == null)
                throw DracoUtils.Failed();
            pointIds = sequencer.GenerateSequence();
            // Initialize point to attribute value mapping for all decoded attributes.
            for (int i = 0; i < NumAttributes; ++i)
            {
                PointAttribute pa = Decoder.PointCloud.Attribute(GetAttributeId(i));
                sequencer.UpdatePointToAttributeIndexMapping(pa);
            }
            base.DecodeAttributes(buffer);
        }

        protected override void DecodePortableAttributes(DecoderBuffer buffer)
        {
            int num_attributes = NumAttributes;
            for (int i = 0; i < num_attributes; ++i)
            {
                sequentialDecoders[i].DecodePortableAttribute(pointIds, buffer);
            }

        }

        protected override void DecodeDataNeededByPortableTransforms(DecoderBuffer buffer)
        {
            int num_attributes = NumAttributes;
            for (int i = 0; i < num_attributes; ++i)
            {
                sequentialDecoders[i].DecodeDataNeededByPortableTransform(pointIds, buffer);
            }

        }

        protected override void TransformAttributesToOriginalFormat()
        {
            int num_attributes = NumAttributes;
            for (int i = 0; i < num_attributes; ++i)
            {
                // Check whether the attribute transform should be skipped.
                if (Decoder.options != null)
                {
                    PointAttribute attribute =
                    sequentialDecoders[i].Attribute;
                    if (Decoder.options.skipAttributeTransform)
                        // attribute.attribute_type(), "skip_attribute_transform", false))
                    {
                        // Attribute transform should not be performed. In this case, we replace
                        // the output geometry attribute with the portable attribute.
                        // TODO(ostava): We can potentially avoid this copy by introducing a new
                        // mechanism that would allow to use the final attributes as portable
                        // attributes for predictors that may need them.
                        sequentialDecoders[i].Attribute.CopyFrom(
                            sequentialDecoders[i].PortableAttribute);
                        continue;
                    }
                }

                sequentialDecoders[i].TransformAttributeToOriginalFormat(pointIds);
            }

        }

        protected virtual SequentialAttributeDecoder CreateSequentialDecoder(SequentialAttributeEncoderType decoderType)
        {
            switch (decoderType)
            {
                case SequentialAttributeEncoderType.Generic:
                    return new SequentialAttributeDecoder();
                case SequentialAttributeEncoderType.Integer:
                    return new SequentialIntegerAttributeDecoder();
                case SequentialAttributeEncoderType.Quantization:
                    return new SequentialQuantizationAttributeDecoder();
                case SequentialAttributeEncoderType.Normals:
                    return new SequentialNormalAttributeDecoder();
                default:
                    //unsupported decoder type
                    return null;
            }
        }

        public override PointAttribute GetPortableAttribute(int attId)
        {
            int loc_id = GetLocalIdForPointAttribute(attId);
            if (loc_id < 0)
                return null;
            return sequentialDecoders[loc_id].PortableAttribute;
        }
    }
}
