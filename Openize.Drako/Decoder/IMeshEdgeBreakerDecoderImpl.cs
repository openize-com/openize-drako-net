using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako.Decoder
{
    interface IMeshEdgeBreakerDecoderImpl
    {

        void Init(MeshEdgeBreakerDecoder decoder);

        MeshAttributeCornerTable GetAttributeCornerTable(int attId);
        MeshAttributeIndicesEncodingData GetAttributeEncodingData(int attId);
        void CreateAttributesDecoder(int attDecoderId);
        void DecodeConnectivity();
        void OnAttributesDecoded();

        MeshEdgeBreakerDecoder GetDecoder();
        CornerTable CornerTable { get; }
    }
}
