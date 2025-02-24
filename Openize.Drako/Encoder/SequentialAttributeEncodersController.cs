using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Encoder
{

    /// <summary>
    /// A basic implementation of an attribute encoder that can be used to encode
    /// an arbitrary set of attributes. The encoder creates a sequential attribute
    /// encoder for each encoded attribute (see sequentialAttributeEncoder.h) and
    /// then it encodes all attribute values in an order defined by a point sequence
    /// generated in the GeneratePointSequence() method. The default implementation
    /// generates a linear sequence of all points, but derived classes can generate
    /// any custom sequence.
    /// </summary>
    class SequentialAttributeEncodersController : AttributesEncoder
    {
        private int attId;
        private PointsSequencer sequencer;

        private SequentialAttributeEncoder[] sequentialEncoders;
        private int[] pointIds;

        public SequentialAttributeEncodersController(PointsSequencer sequencer, int attId)
            :base(attId)
        {
            this.sequencer = sequencer;
            this.attId = attId;
        }

        public override void Initialize(PointCloudEncoder encoder, DracoPointCloud pc)
        {
            base.Initialize(encoder, pc);
            CreateSequentialEncoders();
            // Initialize all value encoders.
            for (int i = 0; i < NumAttributes; ++i)
            {
                int attId = GetAttributeId(i);
                sequentialEncoders[i].Initialize(encoder, attId);
            }
        }

        public override bool MarkParentAttribute(int pointAttributeId)
        {
            int loc_id = GetLocalIdForPointAttribute(pointAttributeId);
            if (loc_id < 0)
                return false;
            sequentialEncoders[loc_id].MarkParentAttribute();
            return true;
        }

        /// <summary>
        /// Creates all sequential encoders (one for each attribute associated with the
        /// encoder).
        /// </summary>
        protected virtual void CreateSequentialEncoders()
        {
            sequentialEncoders = new SequentialAttributeEncoder[NumAttributes];
            for (int i = 0; i < NumAttributes; ++i)
            {
                sequentialEncoders[i] = CreateSequentialEncoder(i);
                if (sequentialEncoders[i] == null)
                    throw DracoUtils.Failed();
            }
        }

        /// <summary>
        /// Create a sequential encoder for a given attribute based on the attribute
        /// type
        /// and the provided encoder options.
        /// </summary>
        protected virtual SequentialAttributeEncoder CreateSequentialEncoder(int i)
        {

            int attId = GetAttributeId(i);
            PointAttribute att = Encoder.PointCloud.Attribute(attId);

            switch (att.DataType)
            {
                case DataType.UINT8:
                case DataType.INT8:
                case DataType.UINT16:
                case DataType.INT16:
                case DataType.UINT32:
                case DataType.INT32:
                    return new SequentialIntegerAttributeEncoder();
                case DataType.FLOAT32:
                {
                    int quantBits = 0;
                    var opts = Encoder.Options;
                        if(att.AttributeType == AttributeType.Normal)
                        quantBits = opts.NormalBits;
                    switch (att.AttributeType)
                    {
                        case AttributeType.Normal:
                            quantBits = opts.NormalBits;
                            break;
                        case AttributeType.Color:
                            quantBits = opts.ColorBits;
                            break;
                        case AttributeType.Position:
                            quantBits = opts.PositionBits;
                            break;
                        case AttributeType.TexCoord:
                            quantBits = opts.TextureCoordinateBits;
                            break;
                    }
                    if (quantBits > 0)
                    {
                        if (att.AttributeType == AttributeType.Normal)
                        {
                            // We currently only support normals with float coordinates
                            // and must be quantized.
                            return new SequentialNormalAttributeEncoder();
                        }
                        else
                        {
                            return new SequentialQuantizationAttributeEncoder();
                        }
                    }
                    break;
                }
            }
            // Return the default attribute encoder.
            return new SequentialAttributeEncoder();
        }

        public override void EncodeAttributesEncoderData(EncoderBuffer outBuffer)
        {
            base.EncodeAttributesEncoderData(outBuffer);
            // Encode a unique id of every sequential encoder.
            for (int i = 0; i < sequentialEncoders.Length; ++i)
            {
                outBuffer.Encode((byte)sequentialEncoders[i].GetUniqueId());
            }
        }

        public override void EncodeAttributes(EncoderBuffer outBuffer)
        {
            if (sequencer == null)
                throw DracoUtils.Failed();
            pointIds = sequencer.GenerateSequence();
            base.EncodeAttributes(outBuffer);
        }

        protected override void TransformAttributesToPortableFormat()
        {
            for (int i = 0; i < sequentialEncoders.Length; ++i)
            {
                sequentialEncoders[i].TransformAttributeToPortableFormat(pointIds);
            }
        }

        protected override void EncodePortableAttributes(EncoderBuffer out_buffer)
        {
            for (int i = 0; i < sequentialEncoders.Length; ++i)
            {
                sequentialEncoders[i].EncodePortableAttribute(pointIds, out_buffer);
            }
        }

        protected override void EncodeDataNeededByPortableTransforms(EncoderBuffer out_buffer)
        {
            for (int i = 0; i < sequentialEncoders.Length; ++i)
            {
                sequentialEncoders[i].EncodeDataNeededByPortableTransform(
                    out_buffer);
            }
        }

        public override PointAttribute GetPortableAttribute(int parentAttId)
        {
            int loc_id = GetLocalIdForPointAttribute(parentAttId);
            if (loc_id < 0)
                return null;
            return sequentialEncoders[loc_id].portableAttribute;
        }

        public override byte GetUniqueId()
        {
            return (byte)AttributeEncoderType.BASIC;
        }
    }
}
