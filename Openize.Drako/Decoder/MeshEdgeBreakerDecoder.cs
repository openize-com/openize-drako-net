using FileFormat.Drako.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Decoder
{

    /// <summary>
    /// Class for decoding data encoded by MeshEdgeBreakerEncoder.
    /// </summary>
    class MeshEdgeBreakerDecoder : MeshDecoder
    {
        private IMeshEdgeBreakerDecoderImpl impl;

        public override CornerTable GetCornerTable()
        {
            return impl.CornerTable;
        }

        public override MeshAttributeCornerTable GetAttributeCornerTable(int attId)
        {
            return impl.GetAttributeCornerTable(attId);
        }

        public override MeshAttributeIndicesEncodingData GetAttributeEncodingData(int attId)
        {
            return impl.GetAttributeEncodingData(attId);
        }

        protected override void InitializeDecoder()
        {

            byte traversalDecoderType = buffer.DecodeU8();
            impl = null;
            if (traversalDecoderType == 0)
            {
                impl = new MeshEdgeBreakerDecoderImpl(this, new MeshEdgeBreakerTraversalDecoder());
            }
            else if (traversalDecoderType == 1)
            {
                impl = new MeshEdgeBreakerDecoderImpl(this, new MeshEdgeBreakerTraversalPredictiveDecoder());
            }
            else if(traversalDecoderType == 2)
                impl = new MeshEdgeBreakerDecoderImpl(this, new MeshEdgeBreakerTraversalValenceDecoder());
            else
                throw DracoUtils.Failed();
        }

        protected override void CreateAttributesDecoder(int attDecoderId)
        {
          impl.CreateAttributesDecoder(attDecoderId);
        }

        protected override void DecodeConnectivity()
        {
          impl.DecodeConnectivity();
        }

        protected override void OnAttributesDecoded()
        {
          impl.OnAttributesDecoded();
        }
    }
}
