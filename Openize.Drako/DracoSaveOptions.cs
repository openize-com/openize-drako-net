using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako
{

    /// <summary>
    /// Compression level for draco file
    /// </summary>
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    enum DracoCompressionLevel
    {
        /// <summary>
        /// No compression, this will result in the minimum encoding time.
        /// </summary>
        NoCompression,
        /// <summary>
        /// Encoder will perform a compression as quickly as possible.
        /// </summary>
        Fast,
        /// <summary>
        /// Standard mode, with good compression and speed.
        /// </summary>
        Standard,
        /// <summary>
        /// Encoder will compress the scene optimally, which may takes longer time to finish.
        /// </summary>
        Optimal
    }

    /// <summary>
    /// Save options for Google draco files
    /// </summary>
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    class DracoEncodeOptions 
    {
        /// <summary>
        /// Quantization bits for position, default value is 14
        /// </summary>
        public int PositionBits { get; set; }
        /// <summary>
        /// Quantization bits for texture coordinate, default value is 12
        /// </summary>
        public int TextureCoordinateBits { get; set; }
        /// <summary>
        /// Quantization bits for vertex color, default value is 10
        /// </summary>
        public int ColorBits { get; set; }
        /// <summary>
        /// Quantization bits for normal vectors, default value is 10
        /// </summary>
        public int NormalBits { get; set; }
        /// <summary>
        /// Compression level, default value is <see cref="DracoCompressionLevel.Standard"/>
        /// </summary>
        public DracoCompressionLevel CompressionLevel { get; set; }

        /// <summary>
        /// Export the scene as point cloud, default value is false.
        /// </summary>
        public bool PointCloud { get; set; } 

        internal bool SplitMeshOnSeams { get; set; }
        //value is defined in Encoding Tagged/Raw
        internal int? SymbolEncodingMethod { get; set; }

        /*
             Compression Level/    Encoder/                 Predictive Scheme
             NoCompression         Sequential               None
             Fast                  Edgebreaker              Difference
             Normal                EdgeBreaker              Parallelogram
             Best                  Predictive Edgebreaker   MultiParallelogram
         */


        internal bool UseBuiltinAttributeCompression = true;

        /// <summary>
        /// Construct a default configuration for saving draco files.
        /// </summary>
        public DracoEncodeOptions()
        {

            PositionBits = 11;
            TextureCoordinateBits = 12;
            NormalBits = 10;
            ColorBits = 10;
            CompressionLevel = DracoCompressionLevel.Standard;
        }

        internal PredictionSchemeMethod GetPredictionMethod(EncodedGeometryType geometryType, PointAttribute attr)
        {
            /*
             * 
             * 
             **/
            if (CompressionLevel == DracoCompressionLevel.NoCompression)
                return PredictionSchemeMethod.None; //No prediction is required when fastest speed is requested.
            if (geometryType == EncodedGeometryType.TriangularMesh)
            {
                if (attr.AttributeType == AttributeType.TexCoord)
                {
                    if(CompressionLevel != DracoCompressionLevel.Fast && CompressionLevel != DracoCompressionLevel.NoCompression)
                        return PredictionSchemeMethod.TexCoordsPortable;
                }
                // Use speed setting to select the best encoding method.
                if (CompressionLevel== DracoCompressionLevel.Fast)
                    return PredictionSchemeMethod.Difference;
                if (CompressionLevel == DracoCompressionLevel.Standard)
                    return PredictionSchemeMethod.Parallelogram;
                return PredictionSchemeMethod.MultiParallelogram;
            }
            return PredictionSchemeMethod.Undefined;
        }
        internal int GetQuantizationBits(PointAttribute attribute)
        {
            switch (attribute.AttributeType)
            {
                case AttributeType.Color:
                    return ColorBits;
                case AttributeType.Normal:
                    return NormalBits;
                case AttributeType.Position:
                    return PositionBits;
                case AttributeType.TexCoord:
                    return TextureCoordinateBits;
                default:
                    throw new Exception("Not supported quantization bits option for the specified attribute type");
            }
        }

        internal PredictionSchemeMethod GetAttributePredictionScheme(PointAttribute attribute)
        {
            return PredictionSchemeMethod.Undefined;
        }

        internal int GetCompressionLevel()
        {
            switch (CompressionLevel)
            {
                case DracoCompressionLevel.NoCompression:
                    return 2;
                case DracoCompressionLevel.Fast:
                    return 6;
                case DracoCompressionLevel.Standard:
                    return 7;
                case DracoCompressionLevel.Optimal:
                    return 10;
            }

            return 7;
        }

        internal int GetSpeed()
        {
            return 3;
        }
    }
}
