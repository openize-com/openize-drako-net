using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako
{
    enum EncodedGeometryType
    {
        Invalid= -1,
        PointCloud = 0,
        TriangularMesh,
    }

    /// <summary>
    /// A variable length encoding for storing all possible topology configurations
    /// during traversal of mesh's surface. The configurations are based on visited
    /// state of neighboring triangles around a currently processed face corner.
    /// Note that about half of the encountered configurations is expected to be of
    /// type TOPOLOGYC. It's guaranteed that the encoding will use at most 2 bits
    /// per triangle for meshes with no holes and up to 6 bits per triangle for
    /// general meshes. In addition, the encoding will take up to 4 bits per triangle
    /// for each non-position attribute attached to the mesh.
    ///
    ///     *-------*          *-------*          *-------*
    ///    / \     / \        / \     / \        / \     / \
    ///   /   \   /   \      /   \   /   \      /   \   /   \
    ///  /     \ /     \    /     \ /     \    /     \ /     \
    /// *-------v-------*  *-------v-------*  *-------v-------*
    ///  \     /x\     /          /x\     /    \     /x\
    ///   \   /   \   /          /   \   /      \   /   \
    ///    \ /  C  \ /          /  L  \ /        \ /  R  \
    ///     *-------*          *-------*          *-------*
    ///
    ///     *       *
    ///    / \     / \
    ///   /   \   /   \
    ///  /     \ /     \
    /// *-------v-------*          v
    ///  \     /x\     /          /x\
    ///   \   /   \   /          /   \
    ///    \ /  S  \ /          /  E  \
    ///     *-------*          *-------*
    ///
    /// </summary>
    enum EdgeBreakerTopologyBitPattern {
      C = 0x0,  // 0
      S = 0x1,  // 1 0 0
      L = 0x3,  // 1 1 0
      R = 0x5,  // 1 0 1
      E = 0x7,  // 1 1 1
      // A special symbol that's not actually encoded, but it can be used to mark
      // the initial face that triggers the mesh encoding of a single connected
      // component.
      InitFace,
      // A special value used to indicate an invalid symbol.
      Invalid
    };

    /// <summary>
    /// Types of edges used during mesh traversal relative to the tip vertex of a
    /// visited triangle.
    /// </summary>
    enum EdgeFaceName : byte
    {
        LeftFaceEdge = 0,
        RightFaceEdge = 1
    }

#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    enum DataType
    {
        // Not a legal value for DataType. Used to indicate a field has not been set.
        INVALID = 0,
        INT8,
        UINT8,
        INT16,
        UINT16,
        INT32,
        UINT32,
        INT64,
        UINT64,
        FLOAT32,
        FLOAT64,
        BOOL,
        TYPESCOUNT
    }

    ///<summary>
    /// List of encoding methods for point clouds.
    ///</summary>
    enum DracoEncodingMethod : int
    {
        Sequential = 0,
        KdTree = 1,
        EdgeBreaker = 1,
    }
    ///<summary>
    /// List of various attribute encoders supported by our framework. The entries
    /// are used as unique identifiers of the encoders and their values should not
    /// be changed!
    ///</summary>
    enum AttributeEncoderType
    {
        BASIC = 0,
        MeshTraversal,
        KdTree,
    }

    ///<summary>
    /// List of various sequential attribute encoder/decoders that can be used in our
    /// pipeline. The values represent unique identifiers used by the decoder and
    /// they should not be changed.
    ///</summary>
    enum SequentialAttributeEncoderType
    {
        Generic = 0,
        Integer,
        Quantization,
        Normals,
    }

    ///<summary>
    /// List of all prediction methods currently supported by our framework.
    ///</summary>
    enum PredictionSchemeMethod
    {
        ///<summary>
        /// Special value indicating that no prediction scheme was used.
        ///</summary>
        None = -2,
        ///<summary>
        /// Used when no specific prediction scheme is required.
        ///</summary>
        Undefined = -1,
        Difference = 0,
        Parallelogram,
        MultiParallelogram,
        TexCoordsDeprecated,
        ConstrainedMultiParallelogram,
        TexCoordsPortable,
        GeometricNormal
        //NUMPREDICTIONSCHEMES
    }

    ///<summary>
    /// List of all prediction scheme transforms used by our framework.
    ///</summary>
    enum PredictionSchemeTransformType
    {
        None = -1,
        ///<summary>
        /// Basic delta transform where the prediction is computed as difference the
        /// predicted and original value.
        ///</summary>
        Delta = 0,
        ///<summary>
        /// An improved delta transform where all computed delta values are wrapped
        /// around a fixed interval which lowers the entropy.
        ///</summary>
        Wrap = 1,
        ///<summary>
        /// Specialized transform for normal coordinates using inverted tiles.
        ///</summary>
        NormalOctahedron = 2,
        ///<summary>
        /// Specialized transform for normal coordinates using canonicalized inverted tiles.
        ///</summary>
        NormalOctahedronCanonicalized = 3,
    }

    enum MeshTraversalMethod : byte
    {
        DepthFirst,
        PredictionDegree
    }
}
