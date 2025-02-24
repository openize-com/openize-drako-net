using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Encoder
{
    /// <summary>
    /// Class implements the edge breaker geometry compression method as described
    /// in "3D Compression Made Simple: Edgebreaker on a Corner-Table" by Rossignac
    /// at al.'01. http://www.cc.gatech.edu/~jarek/papers/CornerTableSMI.pdf
    /// </summary>
    class MeshEdgeBreakerEncoder : MeshEncoder
    {
        private IMeshEdgeBreakerEncoder impl;

        public override CornerTable CornerTable
        {
            get { return impl.CornerTable; }
        }

        public override MeshAttributeCornerTable GetAttributeCornerTable(int attId)
        {
            return impl.GetAttributeCornerTable(attId);
        }

        public override MeshAttributeIndicesEncodingData GetAttributeEncodingData(int attId)
        {
            return impl.GetAttributeEncodingData(attId);
        }

        public override DracoEncodingMethod GetEncodingMethod()
        {
            return DracoEncodingMethod.EdgeBreaker;
        }

        protected override void InitializeEncoder()
        {
            impl = null;
            if (options.CompressionLevel == DracoCompressionLevel.Optimal)
            {
                Buffer.Encode((byte) (1));
                impl = new MeshEdgeBreakerEncoderImpl(new MeshEdgeBreakerTraversalPredictiveEncoder());
            }
            else
            {
                Buffer.Encode((byte) (0));
                impl = new MeshEdgeBreakerEncoderImpl(new MeshEdgeBreakerTraversalEncoder());
            }
            impl.Init(this);
        }

        protected override void EncodeConnectivity()
        {
            impl.EncodeConnectivity();
        }

        protected override void GenerateAttributesEncoder(int attId)
        {
            impl.GenerateAttributesEncoder(attId);
        }

        protected override void EncodeAttributesEncoderIdentifier(int attEncoderId)
        {
            impl.EncodeAttributesEncoderIdentifier(attEncoderId);

        }

    }
}
