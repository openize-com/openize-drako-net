using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako.Encoder
{

    /// <summary>
    /// Default implementation of the edgebreaker traversal encoder. Face
    /// configurations are stored directly into the output buffer and the symbols
    /// are first collected and then encoded in the reverse order to make the
    /// decoding faster.
    /// </summary>
    class MeshEdgeBreakerTraversalEncoder : ITraversalEncoder
    {

        /// <summary>
        /// Bit-length of symbols in the EdgeBreakerTopologyBitPattern stored as a
        /// look up table for faster indexing.
        /// </summary>
        private static readonly int[] edgeBreakerTopologyBitPatternLength =
        {
            1, 3, 0, 3,
            0, 3, 0, 3
        };
        /// <summary>
        /// Buffers for storing encoded data.
        /// </summary>
        private RAnsBitEncoder start_face_encoder_ = new RAnsBitEncoder();
        private EncoderBuffer traversalBuffer = new EncoderBuffer();
        IMeshEdgeBreakerEncoder encoderImpl;
        /// <summary>
        /// Symbols collected during the traversal.
        /// </summary>
        private List<EdgeBreakerTopologyBitPattern> symbols = new List<EdgeBreakerTopologyBitPattern>();
        /// <summary>
        /// Arithmetic encoder for encoding attribute seams.
        /// One context for each non-position attribute.
        /// </summary>
        private RAnsBitEncoder[] attributeConnectivityEncoders;

        public virtual void Init(IMeshEdgeBreakerEncoder encoder)
        {
            encoderImpl = encoder;
        }

        /// <summary>
        /// Called before the traversal encoding is started.
        /// </summary>
        public void Start()
        {
            DracoMesh mesh = encoderImpl.Encoder.Mesh;
            // Allocate enough storage to store initial face configurations. This can
            // consume at most 1 bit per face if all faces are isolated.
            start_face_encoder_.StartEncoding();
            if (mesh.NumAttributes > 1)
            {
                // Init and start arithemtic encoders for storing configuration types
                // of non-position attributes.
                attributeConnectivityEncoders = new RAnsBitEncoder[mesh.NumAttributes - 1];
                for (int i = 0; i < mesh.NumAttributes - 1; ++i)
                {
                    attributeConnectivityEncoders[i] = new RAnsBitEncoder();
                    attributeConnectivityEncoders[i].StartEncoding();
                }
            }
        }

        /// <summary>
        /// Called when a traversal starts from a new initial face.
        /// </summary>
        public void EncodeStartFaceConfiguration(bool interior)
        {
            start_face_encoder_.EncodeBit(interior);
        }

        protected void EncodeStartFaces()
        {
            start_face_encoder_.EndEncoding(traversalBuffer);
        }

        protected void EncodeTraversalSymbols()
        {
            // Bit encode the collected symbols.
            // Allocate enough storage for the bit encoder.
            // It's guaranteed that each face will need only up to 3 bits.
            traversalBuffer.StartBitEncoding(encoderImpl.Encoder.Mesh.NumFaces * 3, true);
            for (int i = symbols.Count - 1; i >= 0; --i)
            {
                traversalBuffer.EncodeLeastSignificantBits32(edgeBreakerTopologyBitPatternLength[(int)symbols[i]], (uint)symbols[i]);
            }
            traversalBuffer.EndBitEncoding();
        }

        protected void EncodeAttributeSeams()
        {
            if (attributeConnectivityEncoders != null)
            {
                for (int i = 0; i < attributeConnectivityEncoders.Length; ++i)
                {
                    attributeConnectivityEncoders[i].EndEncoding(traversalBuffer);
                }
            }
        }

        /// <summary>
        /// Called when a new corner is reached during the traversal. No-op for the
        /// default encoder.
        /// </summary>
        public virtual void NewCornerReached(int corner)
        {
        }

        /// <summary>
        /// Called whenever a new symbol is reached during the edgebreaker traversal.
        /// </summary>
        public virtual void EncodeSymbol(EdgeBreakerTopologyBitPattern symbol)
        {
            // Store the symbol. It will be encoded after all symbols are processed.
            symbols.Add(symbol);
        }

        /// <summary>
        /// Called for every pair of connected and visited faces. |isSeam| specifies
        /// whether there is an attribute seam between the two faces.
        /// </summary>
        public void EncodeAttributeSeam(int attribute, bool isSeam)
        {
            attributeConnectivityEncoders[attribute].EncodeBit(isSeam);
        }

        /// <summary>
        /// Called when the traversal is finished.
        /// </summary>
        public virtual void Done()
        {
            EncodeTraversalSymbols();
            EncodeStartFaces();
            EncodeAttributeSeams();
        }

        /// <summary>
        /// Returns the number of encoded symbols.
        /// </summary>
        public virtual int NumEncodedSymbols
        {
            get { return symbols.Count; }
        }

        public EncoderBuffer Buffer
        {
            get { return traversalBuffer; }
        }

        protected EncoderBuffer OutputBuffer
        {
            get { return traversalBuffer; }
        }

    }
}
