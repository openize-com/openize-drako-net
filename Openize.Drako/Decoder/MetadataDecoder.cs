using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Decoder
{
    class MetadataDecoder
    {
        class MetadataTuple
        {
            public Metadata parent;
            public Metadata decoded;
            public int level;
            public MetadataTuple(Metadata parent, Metadata decoded, int level)
            {
                this.parent = parent;
                this.decoded = decoded;
                this.level = level;
            }
        }

        const int kMaxSubmetadataLevel = 1000;

        public GeometryMetadata Decode(DecoderBuffer buffer)
        {
            uint numAttrMetadata = buffer.DecodeVarintU32();
            var metadata = new GeometryMetadata();
            for(int i = 0; i < numAttrMetadata; i++)
            {
                uint attUniqueId = buffer.DecodeVarintU32();
                var attMetadata = new Metadata();
                DecodeMetadata(buffer, attMetadata);
                metadata.AttributeMetadata[(int)attUniqueId] = attMetadata;
            }

            DecodeMetadata(buffer, metadata);
            return metadata;
        }
        private void DecodeMetadata(DecoderBuffer buffer, Metadata metadata)
        {
            var stack = new List<MetadataTuple>();
            stack.Add(new MetadataTuple(null, metadata, 0));
            while(stack.Count > 0)
            {
                var mp = stack[stack.Count - 1];
                stack.RemoveAt(stack.Count - 1);
                metadata = mp.decoded;
                if(mp.parent != null)
                {
                    if (mp.level > kMaxSubmetadataLevel)
                        throw DracoUtils.Failed();
                    var subMetadataName = DecodeName(buffer);
                    if (subMetadataName == null)
                        throw DracoUtils.Failed();
                    var subMetadata = new Metadata();
                    metadata = subMetadata;
                    mp.parent.SubMetadata[subMetadataName] = subMetadata;
                }
                if (metadata == null)
                    throw DracoUtils.Failed();
                uint numEntries = buffer.DecodeVarintU32();
                for(int i = 0; i < numEntries; i++)
                {
                    DecodeEntry(buffer, metadata);
                }
                uint numSubMetadata = buffer.DecodeVarintU32();
                if (numSubMetadata > buffer.RemainingSize)
                    throw DracoUtils.Failed();
                for(var i = 0; i < numSubMetadata; i++)
                {
                    stack.Add(new MetadataTuple(metadata, null, mp.parent != null ? mp.level + 1 : mp.level));
                }
            }
        }
        private void DecodeEntry(DecoderBuffer buffer, Metadata metadata)
        {
            var entryName = DecodeName(buffer);
            if (entryName == null)
                throw DracoUtils.Failed();
            uint dataSize = buffer.DecodeVarintU32();
            if (dataSize == 0 || dataSize > buffer.RemainingSize)
                throw DracoUtils.Failed();
            var entryValue = new byte[dataSize];
            if (!buffer.Decode(entryValue, (int)dataSize))
                throw DracoUtils.Failed();
            metadata.Entries[entryName] = entryValue;
        }
        private string DecodeName(DecoderBuffer buffer)
        {
            uint nameLen = buffer.DecodeVarintU32();
            var bytes = new byte[nameLen];
            if (!buffer.Decode(bytes, bytes.Length))
                return null;
            return Encoding.UTF8.GetString(bytes);
        }
    }
}
