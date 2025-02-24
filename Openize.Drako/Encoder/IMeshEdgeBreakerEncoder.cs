using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako.Encoder
{
    interface IMeshEdgeBreakerEncoder
    {

        void Init(MeshEdgeBreakerEncoder encoder);

        MeshAttributeCornerTable GetAttributeCornerTable(int attId);
        MeshAttributeIndicesEncodingData GetAttributeEncodingData(int attId);
        void GenerateAttributesEncoder(int attId);
        void EncodeAttributesEncoderIdentifier(int attEncoderId);
        void EncodeConnectivity();

        /// <summary>
        /// Returns corner table of the encoded mesh.
        /// </summary>
        CornerTable CornerTable { get; }
        MeshEdgeBreakerEncoder Encoder { get; }
    }
}
