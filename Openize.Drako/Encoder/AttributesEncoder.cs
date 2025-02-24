using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Encoder
{

    /// <summary>
    /// Base class for encoding one or more attributes of a PointCloud (or other
    /// geometry). This base class provides only the basic interface that is used
    /// by the PointCloudEncoder. The actual encoding must be implemented in derived
    /// classes using the EncodeAttributes() method.
    /// </summary>
    abstract class AttributesEncoder
    {

        /// <summary>
        /// List of attribute ids that need to be encoded with this encoder.
        /// </summary>
        private int[] pointAttributeIds;

        /// <summary>
        /// Map between point attribute id and the local id (i.e., the inverse of the
        /// |pointAttributeIds|.
        /// </summary>
        private int[] pointAttributeToLocalIdMap;

        PointCloudEncoder pointCloudEncoder;
        DracoPointCloud pointCloud;

        public AttributesEncoder()
        {

        }

        /// <summary>
        /// Constructs an attribute encoder assosciated with a given point attribute.
        /// </summary>
        public AttributesEncoder(int pointAttribId)
        {
            pointAttributeIds = new int[] {pointAttribId};
            pointAttributeToLocalIdMap = new int[] {0};
        }

        /// <summary>
        /// Called after all attribute encoders are created. It can be used to perform
        /// any custom initialization, including setting up attribute dependencies.
        /// Note: no data should be encoded in this function, because the decoder may
        /// process encoders in a different order from the decoder.
        /// </summary>
        public virtual void Initialize(PointCloudEncoder encoder, DracoPointCloud pc)
        {

            pointCloudEncoder = encoder;
            pointCloud = pc;
        }

        /// <summary>
        /// Encodes data needed by the target attribute decoder.
        /// </summary>
        public virtual void EncodeAttributesEncoderData(EncoderBuffer outBuffer)
        {

            // Encode data about all attributes.
            Encoding.EncodeVarint((uint) NumAttributes, outBuffer);
            for (int i = 0; i < NumAttributes; ++i)
            {
                int attId = pointAttributeIds[i];
                PointAttribute pa = pointCloud.Attribute(attId);
                outBuffer.Encode((byte) (pa.AttributeType));
                outBuffer.Encode((byte) (pa.DataType));
                outBuffer.Encode((byte) (pa.ComponentsCount));
                outBuffer.Encode((byte) (pa.Normalized ? 1 : 0));
                Encoding.EncodeVarint((ushort) (pa.UniqueId), outBuffer);
            }
        }

        /// <summary>
        /// Returns a unique identifier of the given encoder type, that is used during
        /// decoding to ruct the corresponding attribute decoder.
        /// </summary>
        public abstract byte GetUniqueId();

        /// <summary>
        /// Encode attribute data to the target buffer. Needs to be implmented by the
        /// derived classes.
        /// </summary>
        public virtual void EncodeAttributes(EncoderBuffer out_buffer)
        {

            TransformAttributesToPortableFormat();
            EncodePortableAttributes(out_buffer);
            // Encode data needed by portable transforms after the attribute is encoded.
            // This corresponds to the order in which the data is going to be decoded by
            // the decoder.
            EncodeDataNeededByPortableTransforms(out_buffer);
        }

  // Transforms the input attribute data into a form that should be losslessly
  // encoded (transform itself can be lossy).
        protected virtual void TransformAttributesToPortableFormat()
        {
        }

        // Losslessly encodes data of all portable attributes.
  // Precondition: All attributes must have been transformed into portable
  // format at this point (see TransformAttributesToPortableFormat() method).
        protected abstract void EncodePortableAttributes(EncoderBuffer out_buffer);

  // Encodes any data needed to revert the transform to portable format for each
  // attribute (e.g. data needed for dequantization of quantized values).
        protected virtual void EncodeDataNeededByPortableTransforms(EncoderBuffer out_buffer)
        {
        }

        /// <summary>
        /// Returns the number of attributes that need to be encoded before the
        /// specified attribute is encoded.
        /// Note that the attribute is specified by its point attribute id.
        /// </summary>
        public virtual int NumParentAttributes(int pointAttributeId)
        {
            return 0;
        }

        public virtual int GetParentAttributeId(int pointAttributeId,
            int parentI)
        {
            return -1;
        }

        /// <summary>
        /// Marks a given attribute as a parent of another attribute.
        /// </summary>
        public virtual bool MarkParentAttribute(int pointAttributeId)
        {
            return false;
        }


        public void AddAttributeId(int id)
        {
            int[] ids = pointAttributeIds;
            if (ids == null || ids.Length == 0)
                ids = new int[] {id};
            else
            {
                ids = new int[pointAttributeIds.Length + 1];
                Array.Copy(pointAttributeIds, ids, pointAttributeIds.Length);
                ids[ids.Length - 1] = id;
            }
            SetAttributeIds(ids);
        }

        /// <summary>
        /// Sets new attribute point ids (replacing the existing ones).
        /// </summary>
        public void SetAttributeIds(int[] pointAttributeIds)
        {
            this.pointAttributeIds = new int[pointAttributeIds.Length];
            pointAttributeToLocalIdMap = new int[pointAttributeIds.Length];
            for (int i = 0; i < pointAttributeIds.Length; i++)
            {
                this.pointAttributeIds[i] = pointAttributeIds[i];
                pointAttributeToLocalIdMap[i] = this.pointAttributeIds.Length - 1;
            }
        }

        public int GetAttributeId(int i)
        {
            return pointAttributeIds[i];
        }

        public int NumAttributes
        {
            get { return pointAttributeIds == null ? 0 : pointAttributeIds.Length; }
        }

        public PointCloudEncoder Encoder
        {
            get { return pointCloudEncoder; }
        }

        protected int GetLocalIdForPointAttribute(int pointAttributeId)
        {
            if (pointAttributeToLocalIdMap == null)
                return -1;
            int idMapSize = pointAttributeToLocalIdMap.Length;
            if (pointAttributeId >= idMapSize)
                return -1;
            return pointAttributeToLocalIdMap[pointAttributeId];
        }

        public virtual PointAttribute GetPortableAttribute(int parentAttId)
        {
            return null;
        }
    }
}
