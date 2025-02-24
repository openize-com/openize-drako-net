using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Decoder
{

        class PointAttributeInfo
        {
            public PointAttribute attribute;
            public int offset_dimensionality;
            public DataType data_type;
            public int data_size;
            public int num_components;

        public PointAttributeInfo(PointAttribute target_att, int total_dimensionality, DataType data_type, int data_size, int num_components)
        {
            this.attribute = target_att;
            this.offset_dimensionality = total_dimensionality;
            this.data_type = data_type;
            this.data_size = data_size;
            this.num_components = num_components;
        }
    }
// Output iterator that is used to decode values directly into the data buffer
// of the modified PointAttribute.
// The extension of this iterator beyond the DT_UINT32 concerns itself only with
// the size of the data for efficiency, not the type.  DataType is conveyed in
// but is an unused field populated for any future logic/special casing.
// DT_UINT32 and all other 4-byte types are naturally supported from the size of
// data in the kd tree encoder.  DT_UINT16 and DT_UINT8 are supported by way
// of byte copies into a temporary memory buffer.
    class PointAttributeVectorOutputIterator
    {

        // preallocated memory for buffering different data sizes.  Never reallocated.
        byte[] data_;

        private PointAttributeInfo[] attributes_;
        int point_id_;

        public PointAttributeVectorOutputIterator(PointAttributeInfo[] atts)
        {
            attributes_ = atts;


            //DRACO_DCHECK_GE(atts.size(), 1);
            int required_decode_bytes = 0;
            for (int index = 0; index < attributes_.Length; index++)
            {
                var att = attributes_[index];
                required_decode_bytes = Math.Max(required_decode_bytes, att.data_size * att.num_components);
            }

            data_ = new byte[required_decode_bytes];
        }

        public void Next()
        {
            ++point_id_;
        }

        // We do not want to do ANY copying of this constructor so this particular
        // operator is disabled for performance reasons.
        // Self operator++(int) {
        //   Self copy = *this;
        //   ++point_id_;
        //   return copy;
        // }

        // Still needed in some cases.
        // TODO(hemmer): remove.
        // hardcoded to 3 based on legacy usage.
        public void SetTriple(byte[] val)
        {
            //DRACO_DCHECK_EQ(attributes_.size(), 1);  // Expect only ONE attribute.
            var att = attributes_[0];
            PointAttribute attribute = att.attribute;
            attribute.SetAttributeValue(attribute.MappedIndex(point_id_), val, att.offset_dimensionality);
        }

        // Additional operator taking std::vector as argument.
        private byte[] tmp;
        public void Set(int[] val)
        {
            byte[] bytes = tmp == null || tmp.Length != val.Length ? new byte[val.Length * 4] : tmp;
            tmp = bytes;
            int offset = 0;
            for (int i = 0; i < val.Length; i++)
            {
                Unsafe.PutLE32(bytes, offset, (uint)val[i]);
                offset += 4;
            }
            Set(bytes);
        }
        public void Set(byte[] val)
        {
            for (int index = 0; index < attributes_.Length; index++)
            {
                var att = attributes_[index];
                PointAttribute attribute = att.attribute;
                //const uint32_t &offset = std::get<1>(att);
                //const uint32_t &data_size = std::get<3>(att);
                //const uint32_t &num_components = std::get<4>(att);
                //const uint32_t *data_source = val.data() + offset;
                byte[] src = val;
                int src_offset = att.offset_dimensionality;
                if (att.data_size != 4)
                {
                    // handle uint16_t, uint8_t
                    // selectively copy data bytes
                    int dst_offset = 0;
                    for (int i = 0; i < att.num_components; i += 1, dst_offset += att.data_size)
                    {
                        Array.Copy(src, src_offset + i, data_, dst_offset, att.data_size);
                    }

                    // redirect to copied data
                    src = data_;
                    src_offset = 0;
                }

                var avi = attribute.MappedIndex(point_id_);
                if (avi >= attribute.NumUniqueEntries)
                    return;
                attribute.SetAttributeValue(avi, src, 0);
            }
        }
    }
}
