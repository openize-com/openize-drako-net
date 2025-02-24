using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Decoder;
using FileFormat.Drako.Encoder;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Compression
{
    partial class PredictionScheme 
    {
        public static PredictionScheme Create(PredictionSchemeMethod method, PointAttribute att,
            PredictionSchemeTransform transform)
        {
            if (method == PredictionSchemeMethod.None)
                return null;
            return new PredictionSchemeDifference(att, transform);
        }


        public static PredictionScheme Create(PointCloudEncoder encoder, PredictionSchemeMethod method, int attId, PredictionSchemeTransform transform)
        {
            PointAttribute attr = encoder.PointCloud.Attribute(attId);
            if (method == PredictionSchemeMethod.Undefined)
                method = encoder.Options.GetPredictionMethod(encoder.GeometryType, attr);
            if (method == PredictionSchemeMethod.None)
                return null; //No prediction is required when fastest speed is requested.
            if (encoder.GeometryType == EncodedGeometryType.TriangularMesh)
            {

                // Cast the encoder to mesh encoder. This is not necessarily safe if there
                // is some other encoder decides to use TRIANGULARMESH as the return type,
                // but unfortunately there is not nice work around for this without using
                // RTTI (double dispatch and similar conecepts will not work because of the
                // template nature of the prediction schemes).
                MeshEncoder meshEncoder = (MeshEncoder) encoder;
                var ret = CreateMeshPredictionScheme( meshEncoder, method, attId, transform);
                if (ret != null)
                    return ret;
                // Otherwise try to create another prediction scheme.
            }
            return new PredictionSchemeDifference(attr, transform);
        }

        static PredictionScheme CreateMeshPredictionScheme(MeshEncoder source,
            PredictionSchemeMethod method, int attId, PredictionSchemeTransform transform)
        {
            PointAttribute att = source.PointCloud.Attribute(attId);
            if (source.GeometryType == EncodedGeometryType.TriangularMesh &&
                (method == PredictionSchemeMethod.Parallelogram ||
                 method == PredictionSchemeMethod.MultiParallelogram ||
                 method == PredictionSchemeMethod.GeometricNormal ||
                 method == PredictionSchemeMethod.TexCoordsPortable ||
                 method == PredictionSchemeMethod.TexCoordsDeprecated))
            {
                ICornerTable ct = source.CornerTable;
                MeshAttributeIndicesEncodingData encodingData = source.GetAttributeEncodingData(attId);
                if (ct == null || encodingData == null)
                {
                    // No connectivity data found.
                    return null;
                }
                // Connectivity data exists.
                ICornerTable attCt = source.GetAttributeCornerTable(attId);
                var md = new MeshPredictionSchemeData(source.Mesh, attCt ?? ct,
                    encodingData.encodedAttributeValueIndexToCornerMap,
                    encodingData.vertexToEncodedAttributeValueIndexMap);
                return CreateMeshPredictionSchemeInternal(method, att, transform, md);
            }

            return null;
        }

        public static PredictionScheme Create(PointCloudDecoder source, PredictionSchemeMethod method, int attId,
            PredictionSchemeTransform transform)
        {
            if (method == PredictionSchemeMethod.None)
                return null;
            PredictionScheme ret = null;
            if (source.GeometryType == EncodedGeometryType.TriangularMesh)
            {
                ret = CreateMeshPredictionScheme((MeshDecoder)source, method, attId, transform);
            }
            if (ret != null)
                return ret;
            // Create delta decoder.
            PointAttribute att = source.PointCloud.Attribute(attId);
            return new PredictionSchemeDeltaDecoder(att, transform);
        }
        public static PredictionScheme CreateMeshPredictionScheme(MeshDecoder source, PredictionSchemeMethod method, int attId,
            PredictionSchemeTransform transform)
        {

            PointAttribute att = source.PointCloud.Attribute(attId);
            if (method == PredictionSchemeMethod.Parallelogram ||
                method == PredictionSchemeMethod.MultiParallelogram ||
                method == PredictionSchemeMethod.TexCoordsDeprecated ||
                method == PredictionSchemeMethod.TexCoordsPortable ||
                method == PredictionSchemeMethod.GeometricNormal)
            {
                CornerTable ct = source.GetCornerTable();
                MeshAttributeIndicesEncodingData encodingData = source.GetAttributeEncodingData(attId);
                if (ct == null || encodingData == null)
                {
                    // No connectivity data found.
                    return null;
                }
                // Connectivity data exists.
                MeshAttributeCornerTable attCt = source.GetAttributeCornerTable(attId);
                if (attCt != null)
                {
                    MeshPredictionSchemeData md =
                        new MeshPredictionSchemeData(source.Mesh, attCt,
                            encodingData.encodedAttributeValueIndexToCornerMap,
                            encodingData.vertexToEncodedAttributeValueIndexMap);
                    var ret = CreateMeshPredictionSchemeInternal(method, att, transform, md);
                    if (ret != null)
                        return ret;
                }
                else
                {
                    MeshPredictionSchemeData md = new MeshPredictionSchemeData(source.Mesh, ct,
                        encodingData.encodedAttributeValueIndexToCornerMap,
                        encodingData.vertexToEncodedAttributeValueIndexMap);
                    var ret = CreateMeshPredictionSchemeInternal(
                        method, att, transform, md);
                    if (ret != null)
                        return ret;
                }
            }
            return null;
        }

        static PredictionScheme CreateMeshPredictionSchemeInternal(PredictionSchemeMethod method,
            PointAttribute attribute,
            PredictionSchemeTransform transform,
            MeshPredictionSchemeData meshData)
        {
            if (method == PredictionSchemeMethod.Parallelogram)
            {
                return new MeshPredictionSchemeParallelogram(attribute, transform, meshData);
            }
            else if (method == PredictionSchemeMethod.MultiParallelogram)
            {
                return
                    new MeshPredictionSchemeMultiParallelogram(attribute, transform, meshData);
            }
            else if (method == PredictionSchemeMethod.TexCoordsDeprecated)
            {
                return new MeshPredictionSchemeTexCoords(attribute, transform, meshData);
            }
            else if(method == PredictionSchemeMethod.TexCoordsPortable)
                return new MeshPredictionSchemeTexCoordsPortableDecoder(attribute, transform, meshData);
            else if(method == PredictionSchemeMethod.GeometricNormal)
                return new MeshPredictionSchemeGeometricNormal(attribute, transform, meshData);
            return null;
        }


    }
}
