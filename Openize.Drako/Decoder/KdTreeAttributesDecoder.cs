using Openize.Drako.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako.Decoder
{
    class KdTreeAttributesDecoder : AttributesDecoder
    {
        private List<AttributeQuantizationTransform> attribute_quantization_transforms_ =
            new List<AttributeQuantizationTransform>();

        private IntList min_signed_values_ = new IntList();
        private List<PointAttribute> quantized_portable_attributes_ = new List<PointAttribute>();

        protected override void DecodePortableAttributes(DecoderBuffer buffer)
        {
            if (buffer.BitstreamVersion < 23)
                return;

            byte compression_level = buffer.DecodeU8();
            int num_points = Decoder.PointCloud.NumPoints;

            // Decode data using the kd tree decoding into integer (portable) attributes.
            // We first need to go over all attributes and create a new portable storage
            // for those attributes that need it (floating point attributes that have to
            // be dequantized after decoding).

            int num_attributes = NumAttributes;
            int total_dimensionality = 0; // position is a required dimension
            var atts = new PointAttributeInfo[num_attributes];

            for (int i = 0; i < NumAttributes; ++i)
            {
                int att_id = GetAttributeId(i);
                PointAttribute att = Decoder.PointCloud.Attribute(att_id);
                // All attributes have the same number of values and identity mapping
                // between PointIndex and AttributeValueIndex.
                att.Reset(num_points);
                att.IdentityMapping = true;

                PointAttribute target_att = null;
                if (att.DataType == DataType.UINT32 || att.DataType == DataType.UINT16 ||
                    att.DataType == DataType.UINT8)
                {
                    // We can decode to these attributes directly.
                    target_att = att;
                }
                else if (att.DataType == DataType.INT32 || att.DataType == DataType.INT16 ||
                         att.DataType == DataType.INT8)
                {
                    // Prepare storage for data that is used to convert unsigned values back
                    // to the signed ones.
                    for (int c = 0; c < att.ComponentsCount; ++c)
                    {
                        min_signed_values_.Add(0);
                    }

                    target_att = att;
                }
                else if (att.DataType == DataType.FLOAT32)
                {
                    // Create a portable attribute that will hold the decoded data. We will
                    // dequantize the decoded data to the final attribute later on.
                    PointAttribute port_att = new PointAttribute(att.AttributeType, DataType.UINT32, att.ComponentsCount, false,att.ComponentsCount * DracoUtils.DataTypeLength(DataType.UINT32), 0);
                    port_att.IdentityMapping = true;
                    port_att.Reset(num_points);
                    quantized_portable_attributes_.Add(port_att);
                    target_att = port_att;
                }
                else
                {
                    // Unsupported type.
                    throw DracoUtils.Failed();
                }

                // Add attribute to the output iterator used by the core algorithm.
                DataType data_type = target_att.DataType;
                int data_size = Math.Max(0, DracoUtils.DataTypeLength(data_type));
                int num_components = target_att.ComponentsCount;
                atts[i] = new PointAttributeInfo(target_att, total_dimensionality, data_type,
                    data_size, num_components);
                total_dimensionality += num_components;
            }

            var out_it = new PointAttributeVectorOutputIterator(atts);

            var decoder = new DynamicIntegerPointsKdTreeDecoder (compression_level, total_dimensionality);
            decoder.DecodePoints(buffer, out_it);
        }

        protected override void DecodeDataNeededByPortableTransforms(DecoderBuffer in_buffer)
        {
            if (in_buffer.BitstreamVersion >= 23)
            {
                // Decode quantization data for each attribute that need it.
                // TODO(ostava): This should be moved to AttributeQuantizationTransform.
                float[] min_value;
                for (int i = 0; i < NumAttributes; ++i)
                {
                    int att_id = GetAttributeId(i);
                    PointAttribute att =
                        Decoder.PointCloud.Attribute(att_id);
                    if (att.DataType == DataType.FLOAT32)
                    {
                        int num_components = att.ComponentsCount;
                        min_value = new float[num_components];
                        if (!in_buffer.Decode(min_value))
                            throw DracoUtils.Failed();
                        float max_value_dif = in_buffer.DecodeF32();
                        byte quantization_bits = in_buffer.DecodeU8();
                        if (quantization_bits > 31)
                            throw DracoUtils.Failed();
                        AttributeQuantizationTransform transform = new AttributeQuantizationTransform();
                        transform.SetParameters(quantization_bits, min_value, num_components, max_value_dif);
                        int num_transforms = attribute_quantization_transforms_.Count;
                        transform.TransferToAttribute(
                            quantized_portable_attributes_[num_transforms]);
                        attribute_quantization_transforms_.Add(transform);
                    }
                }

                // Decode transform data for signed integer attributes.
                for (int i = 0; i < min_signed_values_.Count; ++i)
                {
                    uint val = in_buffer.DecodeVarintU32();
                    min_signed_values_[i] = (int)val;
                }

                return;
            }
#if DRACO_BACKWARDS_COMPATIBILITY_SUPPORTED
  // Handle old bitstream
  // Figure out the total dimensionality of the point cloud
  var attribute_count = NumAttributes;
  int total_dimensionality = 0;  // position is a required dimension
  PointAttributeInfo[] atts = new PointAttributeInfo[attribute_count];
  for (var attribute_index = 0;
       (uint)(attribute_index) < attribute_count;
       attribute_index += 1)  // increment the dimensionality as needed...
  {
    int att_id = GetAttributeId(attribute_index);
    PointAttribute att = Decoder.PointCloud.Attribute(att_id);
    DataType data_type = att.DataType;
    int data_size = Math.Max(0, DracoUtils.DataTypeLength(data_type));
    int num_components = att.ComponentsCount;
    atts[attribute_index] = new PointAttributeInfo(att, total_dimensionality, data_type, data_size, num_components);
    // everything is treated as 32bit in the encoder.
    total_dimensionality += num_components;
  }

  int att_id = GetAttributeId(0);
  PointAttribute att = Decoder.PointCloud.Attribute(att_id);
  att.IdentityMapping = true;
  // Decode method
  byte method;
  if (!in_buffer.Decode(&method))
    return false;
  if (method == KdTreeAttributesEncodingMethod::kKdTreeQuantizationEncoding) {
    byte compression_level = 0;
    if (!in_buffer.Decode(out compression_level))
      return false;
    uint num_points = 0;
    if (!in_buffer.Decode(out num_points))
      return false;
    att.Reset(num_points);
    FloatPointsTreeDecoder decoder;
    PointAttributeVectorOutputIterator<float> out_it(atts);
    if (!decoder.DecodePointCloud(in_buffer, out_it))
      return false;
  } else if (method == KdTreeAttributesEncodingMethod::kKdTreeIntegerEncoding) {
    byte compression_level = 0;
    if (!in_buffer.Decode(out compression_level))
      return false;
    if (6 < compression_level) {
      LOGE("KdTreeAttributesDecoder: compression level %i not supported.\n",
           compression_level);
      return false;
    }

    uint num_points;
    if (!in_buffer.Decode(out num_points))
      return false;

    for (int attribute_index = 0; attribute_index < attribute_count; attribute_index += 1) {
      int att_id = GetAttributeId(attribute_index);
      PointAttribute attr = Decoder.PointCloud.Attribute(att_id);
      attr.Reset(num_points);
                    attr.IdentityMapping = true;
    };

    PointAttributeVectorOutputIterator out_it = new PointAttributeVectorOutputIterator(atts);

        var decoder = new DynamicIntegerPointsKdTreeDecoder(compression_level, total_dimensionality);
        if (!decoder.DecodePoints(in_buffer, out_it))
          return false;
  } else {
    // Invalid method.
    return false;
  }
  return true;
#else
            throw DracoUtils.Failed();
#endif
        }

        void TransformAttributeBackToSignedType_short(PointAttribute att, int num_processed_signed_components)
        {
            ushort[] unsigned_val = new ushort[att.ComponentsCount];
            ushort[] signed_val = new ushort[att.ComponentsCount];

            for (int avi = 0; avi < att.NumUniqueEntries; ++avi)
            {
                att.GetValue(avi, unsigned_val);
                for (int c = 0; c < att.ComponentsCount; ++c)
                {
                    // Up-cast |unsigned_val| to int32_t to ensure we don't overflow it for
                    // smaller data types.
                    signed_val[c] = (ushort)((int)unsigned_val[c] + min_signed_values_[num_processed_signed_components + c]);
                }
                att.SetAttributeValue(avi, signed_val);
            }
        }
        void TransformAttributeBackToSignedType_sbyte(PointAttribute att, int num_processed_signed_components)
        {
            byte[] unsigned_val = new byte[att.ComponentsCount];
            byte[] signed_val = new byte[att.ComponentsCount];

            for (int avi = 0; avi < att.NumUniqueEntries; ++avi)
            {
                att.GetValue(avi, unsigned_val);
                for (int c = 0; c < att.ComponentsCount; ++c)
                {
                    // Up-cast |unsigned_val| to int32_t to ensure we don't overflow it for
                    // smaller data types.
                    signed_val[c] = (byte)((int)unsigned_val[c] + min_signed_values_[num_processed_signed_components + c]);
                }
                att.SetAttributeValue(avi, signed_val, 0);
            }
        }
        void TransformAttributeBackToSignedType_int(PointAttribute att, int num_processed_signed_components)
        {
            uint[] unsigned_val = new uint[att.ComponentsCount];
            uint[] signed_val = new uint[att.ComponentsCount];

            for (int avi = 0; avi < att.NumUniqueEntries; ++avi)
            {
                att.GetValue(avi, unsigned_val);
                for (int c = 0; c < att.ComponentsCount; ++c)
                {
                    // Up-cast |unsigned_val| to int32_t to ensure we don't overflow it for
                    // smaller data types.
                    signed_val[c] = (uint)((int)unsigned_val[c] + min_signed_values_[num_processed_signed_components + c]);
                }
                att.SetAttributeValue(avi, signed_val);
            }
        }

        protected override void TransformAttributesToOriginalFormat()
        {

            if (quantized_portable_attributes_.Count == 0 && min_signed_values_.Count == 0)
            {
                return;
            }

            int num_processed_quantized_attributes = 0;
            int num_processed_signed_components = 0;
            // Dequantize attributes that needed it.
            for (int i = 0; i < NumAttributes; ++i)
            {
                int att_id = GetAttributeId(i);
                PointAttribute att = Decoder.PointCloud.Attribute(att_id);
                if (att.DataType == DataType.INT32 || att.DataType == DataType.INT16 ||
                    att.DataType == DataType.INT8)
                {
                    uint[] unsigned_val = new uint[att.ComponentsCount];
                    int[] signed_val = new int[att.ComponentsCount];
                    // Values are stored as unsigned in the attribute, make them signed again.
                    if (att.DataType == DataType.INT32)
                    {
                        TransformAttributeBackToSignedType_int(att, num_processed_signed_components);
                    }
                    else if (att.DataType == DataType.INT16)
                    {
                        TransformAttributeBackToSignedType_short(att, num_processed_signed_components);
                    }
                    else if (att.DataType == DataType.INT8)
                    {
                        TransformAttributeBackToSignedType_sbyte(att, num_processed_signed_components);
                    }
                    num_processed_signed_components += att.ComponentsCount;
                }
                else if (att.DataType == DataType.FLOAT32)
                {
                    // TODO(ostava): This code should be probably moved out to attribute
                    // transform and shared with the SequentialQuantizationAttributeDecoder.

                    PointAttribute src_att = quantized_portable_attributes_[num_processed_quantized_attributes];

                    AttributeQuantizationTransform transform = attribute_quantization_transforms_[num_processed_quantized_attributes];

                    num_processed_quantized_attributes++;

                    if (Decoder.options.skipAttributeTransform)
                    {
                        // Attribute transform should not be performed. In this case, we replace
                        // the output geometry attribute with the portable attribute.
                        // TODO(ostava): We can potentially avoid this copy by introducing a new
                        // mechanism that would allow to use the final attributes as portable
                        // attributes for predictors that may need them.
                        att.CopyFrom(src_att);
                        continue;
                    }

                    // Convert all quantized values back to floats.
                    int max_quantized_value = (1 << transform.quantization_bits_) - 1;
                    int num_components = att.ComponentsCount;
                    int entry_size = sizeof(float) * num_components;
                    float[] att_val = new float[num_components];
                    int quant_val_id = 0;
                    int out_byte_pos = 0;
                    Dequantizer dequantizer = new Dequantizer(transform.range_, max_quantized_value);
                    var portable_attribute_data = src_att.Buffer.AsIntArray();
                    for (int j = 0; j < src_att.NumUniqueEntries; ++j)
                    {
                        for (int c = 0; c < num_components; ++c)
                        {
                            float value = dequantizer.DequantizeFloat(portable_attribute_data[quant_val_id++]);
                            value = value + transform.min_values_[c];
                            att_val[c] = value;
                        }

                        // Store the floating point value into the attribute buffer.
                        att.Buffer.Write(out_byte_pos, att_val);
                        out_byte_pos += entry_size;
                    }
                }
            }

        }

        public override PointAttribute GetPortableAttribute(int attId)
        {
            return null;
        }
    }
}
