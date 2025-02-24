using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako
{

    /// <summary>
    /// Data used for encoding and decoding of mesh attributes.
    /// </summary>
    class MeshAttributeIndicesEncodingData
    {

        /// <summary>
        /// Array for storing the corner ids in the order their associated attribute
        /// entries were encoded/decoded. For every encoded attrbiute value entry we
        /// store exactly one corner. I.e., this is the mapping between an encoded
        /// attribute entry ids and corner ids. This map is needed for example by
        /// prediction schemes. Note that not all corners are included in this map,
        /// e.g., if multiple corners share the same attribute value, only one of these
        /// corners will be usually included.
        /// </summary>
        internal IntList encodedAttributeValueIndexToCornerMap = new IntList();

        /// <summary>
        /// Map for storing encoding order of attribute entries for each vertex.
        /// i.e. Mapping between vertices and their corresponding attribute entry ids
        /// that are going to be used by the decoder.
        /// -1 if an attribute entry hasn't been encoded/decoded yet.
        /// </summary>
        internal int[] vertexToEncodedAttributeValueIndexMap = new int[0];

        /// <summary>
        /// Total number of encoded/decoded attribute entries.
        /// </summary>
        internal int numValues;
    }
}
