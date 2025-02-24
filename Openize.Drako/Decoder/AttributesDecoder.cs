using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Decoder
{


    /// <summary>
    /// Base class for decoding one or more attributes that were encoded with a
    /// matching AttributesEncoder. This base class provides only the basic interface
    /// that is used by the PointCloudDecoder. The actual encoding must be
    /// implemented in derived classes using the DecodeAttributes() method.
    /// </summary>
    abstract class AttributesDecoder
    {
        private int[] pointAttributeIds = new int[0];

        // Map between point attribute id and the local id (i.e., the inverse of the
        // |point_attribute_ids_|.
        private IntList point_attribute_to_local_id_map_ = new IntList();
        private PointCloudDecoder pointCloudDecoder;
        private DracoPointCloud pointCloud;


        public virtual void Initialize(PointCloudDecoder decoder, DracoPointCloud pointCloud)
        {
            this.pointCloud = pointCloud;
            this.pointCloudDecoder = decoder;
        }

        public virtual void DecodeAttributesDecoderData(DecoderBuffer inBuffer)
        {
            // Decode and create attributes.
            int numAttributes;
            if (pointCloudDecoder.BitstreamVersion < 20)
            {
                numAttributes = inBuffer.DecodeI32();
            }
            else
            {
                uint n = inBuffer.DecodeVarintU32();
                numAttributes = (int) n;
            }

            if (numAttributes <= 0)
                throw DracoUtils.Failed();

            pointAttributeIds = new int[numAttributes];
            DracoPointCloud pc = pointCloud;
            int version = pointCloudDecoder.BitstreamVersion;
            for (int i = 0; i < numAttributes; ++i)
            {
                // Decode attribute descriptor data.
                byte attType, dataType, componentsCount, normalized;
                attType = inBuffer.DecodeU8();
                dataType = inBuffer.DecodeU8();
                componentsCount = inBuffer.DecodeU8();
                normalized = inBuffer.DecodeU8();
                DataType dracoDt = (DataType) (dataType);

                // Add the attribute to the point cloud
                PointAttribute ga = new PointAttribute();
                ga.AttributeType = (AttributeType) attType;
                ga.ComponentsCount = componentsCount;
                ga.DataType = dracoDt;
                ga.Normalized = normalized > 0;
                ga.ByteStride = DracoUtils.DataTypeLength(dracoDt) * componentsCount;

                ushort customId;
                if (version < 13)
                {
                    customId = inBuffer.DecodeU16();
                }
                else
                {
                    customId = inBuffer.DecodeVarintU16();
                }



                int attId = pc.AddAttribute(ga);
                ga.UniqueId = customId;
                pointAttributeIds[i] = attId;

                // Update the inverse map.
                if (attId >= point_attribute_to_local_id_map_.Count)
                    point_attribute_to_local_id_map_.Resize(attId + 1, -1);
                point_attribute_to_local_id_map_[attId] = i;
            }

        }

        public virtual void DecodeAttributes(DecoderBuffer buffer)
        {
            DecodePortableAttributes(buffer);
            DecodeDataNeededByPortableTransforms(buffer);
            TransformAttributesToOriginalFormat();
        }

        protected abstract void DecodePortableAttributes(DecoderBuffer buffer);

        protected virtual void DecodeDataNeededByPortableTransforms(DecoderBuffer buffer)
        {
        }

        protected virtual void TransformAttributesToOriginalFormat()
        {
        }

        public int GetAttributeId(int i)
        {
            return pointAttributeIds[i];
        }

        public int NumAttributes
        {
            get { return pointAttributeIds.Length; }
        }

        public PointCloudDecoder Decoder
        {
            get { return pointCloudDecoder; }
        }

        public abstract PointAttribute GetPortableAttribute(int attId);

        protected int GetLocalIdForPointAttribute(int point_attribute_id)
        {
            int id_map_size = point_attribute_to_local_id_map_.Count;
            if (point_attribute_id >= id_map_size)
                return -1;
            return point_attribute_to_local_id_map_[point_attribute_id];
        }
    }
}
