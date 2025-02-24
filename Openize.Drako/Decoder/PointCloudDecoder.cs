using Openize.Drako.Utils;
using System;
using System.Collections.Generic;
using System.Text;

namespace Openize.Drako.Decoder
{
    class DracoLoadOptions
    {
        public bool skipAttributeTransform = false;
        public bool flipTextureCoordinateV = true;

    }
    abstract class PointCloudDecoder
    {
        protected DecoderBuffer buffer;
        protected DracoPointCloud pointCloud;  
        protected AttributesDecoder[] attributesDecoders = new AttributesDecoder[0];
        protected EncodedGeometryType geometryType;
        internal DracoLoadOptions options = new DracoLoadOptions();


  // Map between attribute id and decoder id.
        private int[] attributeToDecoderMap;
        public int BitstreamVersion { get; set; }

        protected PointCloudDecoder(EncodedGeometryType geometryType)
        {
            this.geometryType = geometryType;
        }

        public DracoPointCloud PointCloud
        {
            get { return pointCloud; }
        }

        public EncodedGeometryType GeometryType
        {
            get { return geometryType; }
        }

        public virtual void Decode(DracoHeader header, DecoderBuffer buffer, DracoPointCloud result, bool decodeData)
        {
            this.buffer = buffer;
            this.pointCloud = result;
            BitstreamVersion = header.version;
            if (header.version >= 13 && (header.flags & DracoHeader.MetadataFlagMask) == DracoHeader.MetadataFlagMask)
            {
                DecodeMetadata();
            }
            InitializeDecoder();
            DecodeGeometryData();
            DecodePointAttributes(decodeData);
        }

        private void DecodeMetadata()
        {
            MetadataDecoder decoder = new MetadataDecoder();
            var metadata = decoder.Decode(buffer);
            pointCloud.Metadatas.Add(metadata);
        }

        public void SetAttributesDecoder(int attDecoderId, AttributesDecoder decoder)
        {
            if (attDecoderId < 0)
                throw DracoUtils.Failed();
            if (attDecoderId >= attributesDecoders.Length)
            {
                Array.Resize(ref attributesDecoders, attDecoderId + 1);
            }
            attributesDecoders[attDecoderId] = decoder;
        }

        public AttributesDecoder[] AttributesDecoders
        {
            get { return attributesDecoders;}
            
        }

        protected virtual void InitializeDecoder()
        {
        }

        protected virtual void DecodeGeometryData()
        {
        }

        protected void DecodePointAttributes(bool decodeAttributeData)
        {
            byte numAttributesDecoders = buffer.DecodeU8();
            //create attributes decoders
            for (int i = 0; i < numAttributesDecoders; i++)
            {
                CreateAttributesDecoder(i);
            }
            //initialize all decoders
            foreach(AttributesDecoder dec in attributesDecoders)
            {
                dec.Initialize(this, pointCloud);
            }
            //decode data
            foreach (AttributesDecoder dec in attributesDecoders)
            {
                dec.DecodeAttributesDecoderData(buffer);
            }

  // Create map between attribute and decoder ids.
            int maxAttrId = -1;

            for (int i = 0; i < attributesDecoders.Length; ++i)
            {
                int numAttributes = attributesDecoders[i].NumAttributes;
                for (int j = 0; j < numAttributes; ++j)
                {
                    int attId = attributesDecoders[i].GetAttributeId(j);
                    maxAttrId = Math.Max(attId, maxAttrId);
                }
            }
            attributeToDecoderMap = new int[maxAttrId + 1];
            for (int i = 0; i < attributesDecoders.Length; ++i)
            {
                int numAttributes = attributesDecoders[i].NumAttributes;
                for (int j = 0; j < numAttributes; ++j)
                {
                    int attId = attributesDecoders[i].GetAttributeId(j);
                    attributeToDecoderMap[attId] = i;
                }
            }

            //decode attributes
            if (decodeAttributeData)
            {
                DecodeAllAttributes();
            }
            OnAttributesDecoded();
        }

        protected abstract void CreateAttributesDecoder(int attrDecoderId);
        

        protected virtual void OnAttributesDecoded()
        {
        }
        protected virtual void DecodeAllAttributes()
        {
            foreach (AttributesDecoder dec in attributesDecoders)
            {
                dec.DecodeAttributes(buffer);
            }
        }

        public DecoderBuffer Buffer
        {
            get { return buffer; }
            set { buffer = value; }
        }

        public PointAttribute GetPortableAttribute(int attId)
        {
            if (attId < 0 || attId >= pointCloud.NumAttributes)
                return null;
            int parentAttDecoderId = attributeToDecoderMap[attId];
            return attributesDecoders[parentAttDecoderId].GetPortableAttribute(attId);
        }
    }
}
