using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Encoder
{
    /// <summary>
    /// Encodes all attributes of a given PointCloud using one of the available
    /// Kd-tree compression methods.
    /// </summary>
    class KdTreeAttributesEncoder : AttributesEncoder
    {
        private List<AttributeQuantizationTransform> attribute_quantization_transforms_ =
            new List<AttributeQuantizationTransform>();

        // Min signed values are used to transform signed integers into unsigned ones
        // (by subtracting the min signed value for each component).
        private IntList min_signed_values_ = new IntList();
        private List<PointAttribute> quantized_portable_attributes_ = new List<PointAttribute>();
        private int num_components_;

        public KdTreeAttributesEncoder(int attId)
            : base(attId)
        {

        }

        public override byte GetUniqueId()
        {
            return (byte)AttributeEncoderType.KdTree;
        }

        protected override void TransformAttributesToPortableFormat()
        {
            // Convert any of the input attributes into a format that can be processed by
            // the kd tree encoder (quantization of floating attributes for now).
            int num_points = Encoder.PointCloud.NumPoints;
            int num_components = 0;
            for (int i = 0; i < NumAttributes; ++i)
            {
                int att_id = GetAttributeId(i);
                PointAttribute att = Encoder.PointCloud.Attribute(att_id);
                num_components += att.ComponentsCount;
            }

            num_components_ = num_components;

            // Go over all attributes and quantize them if needed.
            for (int i = 0; i < NumAttributes; ++i)
            {
                int att_id = GetAttributeId(i);
                PointAttribute att = Encoder.PointCloud.Attribute(att_id);
                if (att.DataType == DataType.FLOAT32)
                {
                    // Quantization path.
                    AttributeQuantizationTransform attribute_quantization_transform =
                        new AttributeQuantizationTransform();
                    int quantization_bits = Encoder.Options.GetQuantizationBits(att);
                    if (quantization_bits < 1)
                        throw DracoUtils.Failed();
                    /*
                    if (Encoder.Options().IsAttributeOptionSet(att_id,
                            "quantization_origin") &&
                        Encoder.Options().IsAttributeOptionSet(att_id,
                            "quantization_range"))
                    {
                        // Quantization settings are explicitly specified in the provided
                        // options.

                        float[] quantization_origin = new float[att.ComponentsCount];
                        Encoder.Options.GetAttributeVector(att_id, "quantization_origin",
                            att.ComponentsCount,
                            &quantization_origin[0]);
                        float range = 1.0f; // Encoder.Options.GetAttributeFloat( att_id, "quantization_range", 1.f);
                        attribute_quantization_transform.SetParameters(
                            quantization_bits, quantization_origin,
                            att.ComponentsCount, range);
                    }
                    else
                    */
                    {
                        // Compute quantization settings from the attribute values.
                        attribute_quantization_transform.ComputeParameters(att, quantization_bits);
                    }

                    attribute_quantization_transforms_.Add(attribute_quantization_transform);
                    // Store the quantized attribute in an array that will be used when we do
                    // the actual encoding of the data.
                    var portable_att = attribute_quantization_transform.InitTransformedAttribute(att, num_points);
                        attribute_quantization_transform.GeneratePortableAttribute(att, (int) num_points, portable_att);
                    quantized_portable_attributes_.Add(portable_att);
                }
                else if (att.DataType == DataType.INT32 || att.DataType == DataType.INT16 ||
                         att.DataType == DataType.INT8)
                {
                    // For signed types, find the minimum value for each component. These
                    // values are going to be used to transform the attribute values to
                    // unsigned integers that can be processed by the core kd tree algorithm.
                    int[] min_value = new int[att.ComponentsCount];
                    for (int j = 0; j < min_value.Length; j++)
                        min_value[j] = int.MaxValue;
                    int[] act_value = new int[att.ComponentsCount];
                    for (int avi = 0; avi < att.NumUniqueEntries; ++avi)
                    {
                        att.ConvertValue(avi, out act_value[0]);
                        for (int c = 0; c < att.ComponentsCount; ++c)
                        {
                            if (min_value[c] > act_value[c])
                                min_value[c] = act_value[c];
                        }
                    }

                    for (int c = 0; c < att.ComponentsCount; ++c)
                    {
                        min_signed_values_.Add(min_value[c]);
                    }
                }
            }

        }

        protected override void EncodePortableAttributes(EncoderBuffer out_buffer)
        {

            // Encode the data using the kd tree encoder algorithm. The data is first
            // copied to a PointDVector that provides all the API expected by the core
            // encoding algorithm.

            // We limit the maximum value of compression_level to 6 as we don't currently
            // have viable algorithms for higher compression levels.
            byte compression_level =
                (byte)Math.Min(10 - Encoder.Options.GetSpeed(), 6);

            if (compression_level == 6 && num_components_ > 15)
            {
                // Don't use compression level for CL >= 6. Axis selection is currently
                // encoded using 4 bits.
                compression_level = 5;
            }

            out_buffer.Encode(compression_level);

            // Init PointDVector. The number of dimensions is equal to the total number
            // of dimensions across all attributes.
            int num_points = Encoder.PointCloud.NumPoints;
            int[][] point_vector = new int[num_points][];
            for (int i = 0; i < num_points; i++)
                point_vector[i] = new int[num_components_];

            int num_processed_components = 0;
            int num_processed_quantized_attributes = 0;
            int num_processed_signed_components = 0;
            // Copy data to the point vector.
            for (int i = 0; i < NumAttributes; ++i)
            {
                int att_id = GetAttributeId(i);
                PointAttribute att = Encoder.PointCloud.Attribute(att_id);
                PointAttribute source_att = null;
                if (att.DataType == DataType.UINT32 || att.DataType == DataType.UINT16 ||
                    att.DataType == DataType.UINT8 || att.DataType == DataType.INT32 ||
                    att.DataType == DataType.INT16 || att.DataType == DataType.INT8)
                {
                    // Use the original attribute.
                    source_att = att;
                }
                else if (att.DataType == DataType.FLOAT32)
                {
                    // Use the portable (quantized) attribute instead.
                    source_att = quantized_portable_attributes_[num_processed_quantized_attributes];
                    num_processed_quantized_attributes++;
                }
                else
                {
                    // Unsupported data type.
                    throw DracoUtils.Failed();
                }

                if (source_att == null)
                    throw DracoUtils.Failed();

                // Copy source_att to the vector.
                if (source_att.DataType == DataType.UINT32)
                {
                    // If the data type is the same as the one used by the point vector, we
                    // can directly copy individual elements.
                    for (int pi = 0; pi < num_points; ++pi)
                    {
                        int avi = source_att.MappedIndex(pi);
                        int offset = source_att.GetBytePos(avi);
                        CopyAttribute(point_vector, source_att.ComponentsCount,
                            num_processed_components, pi,
                            source_att.Buffer.GetBuffer(), offset);
                    }
                }
                else if (source_att.DataType == DataType.INT32 ||
                         source_att.DataType == DataType.INT16 ||
                         source_att.DataType == DataType.INT8)
                {
                    // Signed values need to be converted to unsigned before they are stored
                    // in the point vector.
                    int[] signed_point = new int[source_att.ComponentsCount];
                    int[] unsigned_point = new int[source_att.ComponentsCount];
                    for (int pi = 0; pi < num_points; ++pi)
                    {
                        int avi = source_att.MappedIndex(pi);
                        source_att.ConvertValue(avi, out signed_point[0]);
                        for (int c = 0; c < source_att.ComponentsCount; ++c)
                        {
                            unsigned_point[c] =
                                signed_point[c] -
                                min_signed_values_[num_processed_signed_components + c];
                        }

                        CopyAttribute(point_vector, num_processed_components, pi, unsigned_point);
                    }

                    num_processed_signed_components += source_att.ComponentsCount;
                }
                else
                {
                    // If the data type of the attribute is different, we have to convert the
                    // value before we put it to the point vector.
                    int[] point = new int[source_att.ComponentsCount];
                    for (int pi = 0; pi < num_points; ++pi)
                    {
                        int avi = source_att.MappedIndex(pi);
                        source_att.ConvertValue(avi, out point[0]);
                        CopyAttribute(point_vector, num_processed_components, pi, point);
                    }
                }

                num_processed_components += source_att.ComponentsCount;
            }

            // Compute the maximum bit length needed for the kd tree encoding.
            int num_bits = 0;
            for (int r = 0; r < num_points; r++)
            {
                for (int c = 0; c < num_components_; c++)
                {
                    if (point_vector[r][c] > 0)
                    {
                        int msb = DracoUtils.MostSignificantBit((uint) point_vector[r][c]) + 1;
                        if (msb > num_bits)
                        {
                            num_bits = msb;
                        }
                    }
                }
            }

            var points_encoder = new DynamicIntegerPointsKdTreeEncoder(compression_level, num_components_);
            points_encoder.EncodePoints(point_vector, num_bits, out_buffer);
        }

        protected override void EncodeDataNeededByPortableTransforms(EncoderBuffer out_buffer)
        {
            // Store quantization settings for all attributes that need it.
            for (int i = 0; i < attribute_quantization_transforms_.Count; ++i)
            {
                attribute_quantization_transforms_[i].EncodeParameters(out_buffer);
            }

            // Encode data needed for transforming signed integers to unsigned ones.
            for (int i = 0; i < min_signed_values_.Count; ++i)
            {
                Encoding.EncodeVarint(min_signed_values_[i], out_buffer);
            }

        }

        // Copy data directly off of an attribute buffer interleaved into internal
        // memory.
        void CopyAttribute(int[][] attribute,
            // The offset in dimensions to insert this attribute.
            int offset_dimensionality, int index,
            // The direct pointer to the data
            int[] attribute_item_data)
        {
            //int attribute_dimensionality = attribute_item_data.Length;
            // chunk copy
            int copy_size = attribute_item_data.Length;

            // a multiply and add can be optimized away with an iterator
            // std::memcpy(data0_ + index * dimensionality_ + offset_dimensionality, attribute_item_data, copy_size);
            var face = attribute[index];
            for (int j = offset_dimensionality, i = 0; i < copy_size; i++, j++)
            {
                face[j] = attribute_item_data[i];
            }
        }
        void CopyAttribute(int[][] attribute,
            // The dimensionality of the attribute being integrated
            int attribute_dimensionality,
            // The offset in dimensions to insert this attribute.
            int offset_dimensionality, int index,
            // The direct pointer to the data
            byte[] attribute_item_data, int offset)
        {
            //int copy_size = sizeof(int) * attribute_dimensionality;
            int copy_size = attribute_dimensionality;

            // a multiply and add can be optimized away with an iterator
            // std::memcpy(data0_ + index * dimensionality_ + offset_dimensionality, attribute_item_data, copy_size);
            var face = attribute[index];
            for (int j = offset_dimensionality, i = 0; i < copy_size; i++, j++)
            {
                face[j] = (int)Unsafe.GetLE32(attribute_item_data, offset);
                offset += 4;
            }
        }
    }
}
