using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using Openize.Drako.Encoder;
using Openize.Drako.Utils;

namespace Openize.Drako
{
    class AttributeOctahedronTransform : AttributeTransform
    {
        private int quantizationBits;

        public AttributeOctahedronTransform(int quantizationBits)
        {
            this.quantizationBits = quantizationBits;
        }

        protected override DataType GetTransformedDataType(PointAttribute attribute)
        {
            return DataType.UINT32;
        }
        protected override int GetTransformedNumComponents(PointAttribute attribute)
        {
            return 2;
        }

        public override void CopyToAttributeTransformData(AttributeTransformData outData)
        {
            outData.transformType = AttributeTransformType.OctahedronTransform;
            outData.AppendValue(quantizationBits);
        }

        public void EncodeParameters(EncoderBuffer encoder_buffer)
        {
            if (quantizationBits != -1)
            {
                encoder_buffer.Encode((byte)quantizationBits);
            }
            else
                throw DracoUtils.Failed();
        }

        public PointAttribute GeneratePortableAttribute(PointAttribute attribute, int[] point_ids, int num_points)
        {

            // Allocate portable attribute.
            int num_entries = point_ids.Length;
            PointAttribute portable_attribute = InitPortableAttribute(num_entries, 2, num_points, attribute, true);

            //IntArray portable_attribute_data = IntArray.Wrap(portable_attribute.Buffer.GetBuffer(), 0, num_entries * 2 * 4);
            //IntArray portable_attribute_data2 = IntArray.Wrap(portable_attribute.Buffer.GetBuffer(), 0, num_entries * 4);
            Span<byte> buffer = new Span<byte>(portable_attribute.Buffer.GetBuffer(), 0, num_entries * 2 * 4);
            Span<int> portable_attribute_data = MemoryMarshal.Cast<byte, int>(buffer);


            // Quantize all values in the order given by point_ids into portable
            // attribute.
            Span<float> att_val = stackalloc float[3];
            int dst_index = 0;
            OctahedronToolBox converter = new OctahedronToolBox();
            converter.SetQuantizationBits(quantizationBits);
            for (int i = 0; i < point_ids.Length; ++i)
            {
                int att_val_id = attribute.MappedIndex(point_ids[i]);
                attribute.GetValue(att_val_id, att_val);
                // Encode the vector into a s and t octahedral coordinates.
                int s, t;
                converter.FloatVectorToQuantizedOctahedralCoords(att_val, out s, out t);
                portable_attribute_data[dst_index++] = s;
                portable_attribute_data[dst_index++] = t;
            }
            return portable_attribute;
        }
    }
}
