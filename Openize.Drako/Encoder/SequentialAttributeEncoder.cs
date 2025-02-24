using Openize.Drako.Compression;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Encoder
{
    /// <summary>
    /// A base class for encoding attribute values of a single attribute using a
    /// given sequence of point ids. The default implementation encodes all attribute
    /// values directly to the buffer but derived classes can perform any custom
    /// encoding (such as quantization) by overriding the EncodeValues() method.
    /// </summary>
    class SequentialAttributeEncoder
    {
        private PointCloudEncoder encoder;
        protected PointAttribute attribute;
        internal PointAttribute portableAttribute;
        private int attributeId = -1;

        /// <summary>
        /// List of attribute encoders that need to be encoded before this attribute.
        /// E.g. The parent attributes may be used to predict values used by this
        /// attribute encoder.
        /// </summary>
        private IntList parentAttributes = new IntList();

        private bool isParentEncoder;

        /// <summary>
        /// Method that can be used for custom initialization of an attribute encoder,
        /// such as creation of prediction schemes and initialization of attribute
        /// encoder dependencies.
        /// |encoder| is the parent PointCloudEncoder,
        /// |attributeId| is the id of the attribute that is being encoded by this
        /// encoder.
        /// This method is automatically called by the PointCloudEncoder after all
        /// attribute encoders are created and it should not be called explicitly from
        /// other places.
        /// </summary>
        public virtual void Initialize(PointCloudEncoder encoder, int attributeId)
        {

            this.encoder = encoder;
            attribute = encoder.PointCloud.Attribute(attributeId);
            this.attributeId = attributeId;
        }

        /// <summary>
        /// Intialization for a specific attribute. This can be used mostly for
        /// standalone encoding of an attribute without an PointCloudEncoder.
        /// </summary>
        public virtual void InitializeStandalone(PointAttribute attribute)
        {

            this.attribute = attribute;
            attributeId = -1;
        }

        public virtual void TransformAttributeToPortableFormat(int[] point_ids)
        {
            // Default implementation doesn't transform the input data.
        }

        public virtual void EncodePortableAttribute(int[] point_ids, EncoderBuffer out_buffer)
        {
            // Lossless encoding of the input values.
            EncodeValues(point_ids, out_buffer);
        }

        public virtual void EncodeDataNeededByPortableTransform(EncoderBuffer out_buffer)
        {
            // Default implementation doesn't transform the input data.
        }

        protected void SetPredictionSchemeParentAttributes(PredictionScheme ps)
        {
            for (int i = 0; i < ps.NumParentAttributes; ++i)
            {
                int att_id = encoder.PointCloud.GetNamedAttributeId(ps.GetParentAttributeType(i));
                if (att_id == -1)
                    throw DracoUtils.Failed(); // Requested attribute does not exist.
                ps.SetParentAttribute(encoder.GetPortableAttribute(att_id));
            }

        }

        /// <summary>
        /// Encode all attribute values in the order of the provided points.
        /// The actual implementation of the encoding is done in the EncodeValues()
        /// method.
        /// </summary>
        public bool Encode(int[] pointIds, EncoderBuffer outBuffer)
        {
            EncodeValues(pointIds, outBuffer);
            if (isParentEncoder && IsLossyEncoder())
            {
                if (!PrepareLossyAttributeData())
                    return false;
            }
            return true;
        }

        public virtual bool IsLossyEncoder()
        {
            return false;
        }

        public int NumParentAttributes
        {
            get { return parentAttributes.Count; }
        }

        public int GetParentAttributeId(int i)
        {
            return parentAttributes[i];
        }

        /// <summary>
        /// Called when this attribute encoder becomes a parent encoder of another
        /// encoder.
        /// </summary>
        public void MarkParentAttribute()
        {

            isParentEncoder = true;
        }

        public virtual SequentialAttributeEncoderType GetUniqueId()
        {
            return SequentialAttributeEncoderType.Generic;
        }

        public PointAttribute Attribute
        {
            get { return attribute; }
        }

        public int AttributeId
        {
            get { return attributeId; }
        }

        public PointCloudEncoder Encoder
        {
            get { return encoder; }
        }

        /// <summary>
        /// Should be used to initialize newly created prediction scheme.
        /// Returns false when the initialization failed (in which case the scheme
        /// cannot be used).
        /// </summary>
        protected virtual void InitPredictionScheme(PredictionScheme ps)
        {

            for (int i = 0; i < ps.NumParentAttributes; ++i)
            {
                int attId = encoder.PointCloud.GetNamedAttributeId(
                    ps.GetParentAttributeType(i));
                if (attId == -1)
                    throw DracoUtils.Failed(); // Requested attribute does not exist.
                parentAttributes.Add(attId);
                encoder.MarkParentAttribute(attId);
            }
        }

        /// <summary>
        /// Encodes all attribute values in the specified order. Should be overriden
        /// for specialized  encoders.
        /// </summary>
        protected virtual void EncodeValues(int[] pointIds,
            EncoderBuffer outBuffer)
        {

            int entrySize = attribute.ByteStride;
            byte[] valueData = new byte[entrySize];
            // Encode all attribute values in their native raw format.
            for (int i = 0; i < pointIds.Length; ++i)
            {
                int entryId = attribute.MappedIndex(pointIds[i]);
                attribute.GetValue(entryId, valueData);
                outBuffer.Encode(valueData, entrySize);
            }
        }

        /// <summary>
        /// Method that can be used by lossy encoders to compute encoded lossy
        /// attribute data.
        /// If the return value is true, the caller can call either
        /// GetLossyAttributeData() or encodedLossyAttributeData() to get a new
        /// attribute that is filled with lossy version of the original data (i.e.,
        /// the same data that is going to be used by the decoder).
        /// </summary>
        protected virtual bool PrepareLossyAttributeData()
        {
            return false;
        }

        protected bool IsParentEncoder()
        {
            return isParentEncoder;
        }

        protected static PredictionSchemeMethod SelectPredictionMethod(int att_id, PointCloudEncoder encoder)
        {
            int speed = encoder.Options.GetSpeed();
            if (speed >= 10)
            {
                // Selected fastest, though still doing some compression.
                return PredictionSchemeMethod.Difference;
            }

            if (encoder.GeometryType == EncodedGeometryType.TriangularMesh)
            {
                // Use speed setting to select the best encoding method.
                PointAttribute att = encoder.PointCloud.Attribute(att_id);
                if (att.AttributeType == AttributeType.TexCoord)
                {
                    if (speed < 4)
                    {
                        // Use texture coordinate prediction for speeds 0, 1, 2, 3.
                        return PredictionSchemeMethod.TexCoordsPortable;
                    }
                }

                if (att.AttributeType == AttributeType.Normal)
                {
                    if (speed < 4)
                    {
                        // Use geometric normal prediction for speeds 0, 1, 2, 3.
                        return PredictionSchemeMethod.GeometricNormal;
                    }
                    return PredictionSchemeMethod.Difference; // default
                }

                // Handle other attribute types.
                if (speed >= 8)
                {
                    return PredictionSchemeMethod.Difference;
                }
                if (speed >= 2 || encoder.PointCloud.NumPoints < 40)
                {
                    // Parallelogram prediction is used for speeds 2 - 7 or when the overhead
                    // of using constrained multi-parallelogram would be too high.
                    return PredictionSchemeMethod.Parallelogram;
                }
                // Multi-parallelogram is used for speeds 0, 1.
                return PredictionSchemeMethod.ConstrainedMultiParallelogram;
            }
            // Default option is delta coding.
            return PredictionSchemeMethod.Difference;
        }

    }
}
