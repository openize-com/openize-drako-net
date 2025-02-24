using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Encoder
{
    interface ITraversalEncoder
    {
        void EncodeSymbol(EdgeBreakerTopologyBitPattern sym);

        /// <summary>
        /// Called when a new corner is reached during the traversal. No-op for the
        /// default encoder.
        /// </summary>
        void NewCornerReached(int corner);

        /// <summary>
        /// Called for every pair of connected and visited faces. |isSeam| specifies
        /// whether there is an attribute seam between the two faces.
        /// </summary>
        void EncodeAttributeSeam(int attribute, bool isSeam);

        /// <summary>
        /// Called before the traversal encoding is started.
        /// </summary>
        void Start();

        /// <summary>
        /// Called when the traversal is finished.
        /// </summary>
        void Done();

        /// <summary>
        /// Called when a traversal starts from a new initial face.
        /// </summary>
        void EncodeStartFaceConfiguration(bool interior);

        void Init(IMeshEdgeBreakerEncoder encoder);

        EncoderBuffer Buffer { get; }
        int NumEncodedSymbols { get; }

    }
}
