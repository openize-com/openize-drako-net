using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Encoder
{
    /// <summary>
    /// Abstract base class for all point cloud and mesh encoders. It provides a
    /// basic funcionality that's shared between different encoders.
    /// </summary>
    abstract class PointCloudEncoder
    {
        private DracoPointCloud pointCloud;
        private List<AttributesEncoder> attributesEncoders = new List<AttributesEncoder>();

        /// <summary>
        /// Map between attribute id and encoder id.
        /// </summary>
        private int[] attributeToEncoderMap;

        /// <summary>
        /// Encoding order of individual attribute encoders (i.e., the order in which
        /// they are processed during encoding that may be different from the order
        /// in which they were created because of attribute dependencies.
        /// </summary>
        private int[] attributesEncoderIdsOrder;

        /// <summary>
        /// This buffer holds the final encoded data.
        /// </summary>
        protected EncoderBuffer buffer;


        protected DracoEncodeOptions options;

        public PointCloudEncoder()
        {

        }

        /// <summary>
        /// The main entry point that encodes provided point cloud.
        /// </summary>
        /// <exception cref="DrakoException">Raised when failed to encode the point cloud.</exception>
        public void Encode(DracoEncodeOptions options, EncoderBuffer outBuffer)
        {
            this.options = options;
            buffer = outBuffer;

            // Cleanup from previous runs.
            attributesEncoders.Clear();
            attributeToEncoderMap = null;
            attributesEncoderIdsOrder = null;

            if (pointCloud == null)
                throw DracoUtils.Failed();
            InitializeEncoder();
            EncodeEncoderData();
            EncodeGeometryData();
            EncodePointAttributes();
        }

        public virtual EncodedGeometryType GeometryType
        {
            get { return EncodedGeometryType.PointCloud; }
        }

        /// <summary>
        /// Returns the unique identifier of the encoding method (such as Edgebreaker
        /// for mesh compression).
        /// </summary>
        public abstract DracoEncodingMethod GetEncodingMethod();

        public int NumAttributesEncoders
        {
            get { return attributesEncoders.Count; }
        }

        public AttributesEncoder AttributesEncoder(int i)
        {
            return attributesEncoders[i];
        }

        /// <summary>
        /// Adds a new attribute encoder, returning its id.
        /// </summary>
        public int AddAttributesEncoder(AttributesEncoder attEnc)
        {
            attributesEncoders.Add(attEnc);
            return attributesEncoders.Count - 1;
        }

        /// <summary>
        /// Marks one attribute as a parent of another attribute. Must be called after
        /// all attribute encoders are created (usually in the
        /// AttributeEncoder::Initialize() method).
        /// </summary>
        public bool MarkParentAttribute(int parentAttId)
        {

            if (parentAttId < 0 || parentAttId >= pointCloud.NumAttributes)
                return false;
            int parentAttEncoderId =
                attributeToEncoderMap[parentAttId];
            if (!attributesEncoders[parentAttEncoderId].MarkParentAttribute(
                parentAttId))
                return false;
            return true;
        }


        public EncoderBuffer Buffer
        {
            get { return buffer; }
        }

        public DracoEncodeOptions Options
        {
            get { return options; }
        }

        public DracoPointCloud PointCloud
        {
            get { return pointCloud; }
            set { pointCloud = value; }
        }

        /// <summary>
        /// Can be implemented by derived classes to perform any custom initialization
        /// of the encoder. Called in the Encode() method.
        /// </summary>
        protected virtual void InitializeEncoder()
        {
        }

        /// <summary>
        /// Should be used to encode any encoder-specific data.
        /// </summary>
        protected virtual void EncodeEncoderData()
        {
        }

        /// <summary>
        /// Encodes any global geometry data (such as the number of points).
        /// </summary>
        protected virtual void EncodeGeometryData()
        {
        }

        /// <summary>
        /// encode all attribute values. The attribute encoders are sorted to resolve
        /// any attribute dependencies and all the encoded data is stored into the
        /// |buffer|.
        /// Returns false if the encoding failed.
        /// </summary>
        protected virtual void EncodePointAttributes()
        {

            GenerateAttributesEncoders();

            // Encode the number of attribute encoders.
            buffer.Encode((byte) (attributesEncoders.Count));

            // Initialize all the encoders (this is used for example to init attribute
            // dependencies, no data is encoded in this step).
            foreach (var attEnc in attributesEncoders)
            {
                attEnc.Initialize(this, pointCloud);
            }

            // Rearrange attributes to respect dependencies between individual attributes.
            RearrangeAttributesEncoders();

            // Encode any data that is necessary to create the corresponding attribute
            // decoder.
            for (int i = 0; i < attributesEncoderIdsOrder.Length; i++)
            {
                int attEncoderId = attributesEncoderIdsOrder[i];
                EncodeAttributesEncoderIdentifier(attEncoderId);
            }

            // Also encode any attribute encoder data (such as the info about encoded
            // attributes).
            for (int i = 0; i < attributesEncoderIdsOrder.Length; i++)
            {
                int attEncoderId = attributesEncoderIdsOrder[i];
                attributesEncoders[attEncoderId].EncodeAttributesEncoderData(
                    buffer);
            }

            // Lastly encode all the attributes using the provided attribute encoders.
            EncodeAllAttributes();
        }

        /// <summary>
        /// Generate attribute encoders that are going to be used for encoding
        /// point attribute data. Calls GenerateAttributesEncoder() for every attribute
        /// of the encoded PointCloud.
        /// </summary>
        /// <exception cref="DrakoException">throws when failed to generate encoders</exception>
        protected virtual void GenerateAttributesEncoders()
        {

            for (int i = 0; i < pointCloud.NumAttributes; ++i)
            {
                GenerateAttributesEncoder(i);
            }
            attributeToEncoderMap = new int[pointCloud.NumAttributes];
            for (int i = 0; i < attributesEncoders.Count; ++i)
            {
                for (int j = 0; j < attributesEncoders[i].NumAttributes; ++j)
                {
                    attributeToEncoderMap[attributesEncoders[i].GetAttributeId(j)] = i;
                }
            }
        }

        /// <summary>
        /// Creates attribute encoder for a specific point attribute. This function
        /// needs to be implemented by the derived classes. The derived classes need
        /// to either 1. Create a new attribute encoder and add it using the
        /// AddAttributeEncoder method, or 2. add the attribute to an existing
        /// attribute encoder (using AttributesEncoder::AddAttributeId() method).
        /// </summary>
        protected abstract void GenerateAttributesEncoder(int attId);

        /// <summary>
        /// Encodes any data that is necessary to recreate a given attribute encoder.
        /// Note: this is called in order in which the attribute encoders are going to
        /// be encoded.
        /// </summary>
        protected virtual void EncodeAttributesEncoderIdentifier(int attEncoderId)
        {
        }

        /// <summary>
        /// Encodes all the attribute data using the created attribute encoders.
        /// </summary>
        protected virtual void EncodeAllAttributes()
        {

            for (int i = 0; i < attributesEncoderIdsOrder.Length; i++)
            {
                int attEncoderId = attributesEncoderIdsOrder[i];
                attributesEncoders[attEncoderId].EncodeAttributes(buffer);
            }
        }

        /// <summary>
        /// Rearranges attribute encoders and their attributes to reflect the
        /// underlying attribute dependencies. This ensures that the attributes are
        /// encoded in the correct order (parent attributes before their children).
        /// </summary>
        private void RearrangeAttributesEncoders()
        {

            // Find the encoding order of the attribute encoders that is determined by
            // the parent dependencies between individual encoders. Instead of traversing
            // a graph we encode the attributes in multiple iterations where encoding of
            // attributes that depend on other attributes may get posponed until the
            // parent attributes are processed.
            // This is simpler to implement than graph traversal and it automatically
            // detects any cycles in the dependency graph.
            // TODO(ostava): Current implementation needs to encode all attributes of a
            // single encoder to be encoded in a single "chunk", therefore we need to sort
            // attribute encoders before we sort individual attributes. This requirement
            // can be lifted for encoders that can encode individual attributes separately
            // but it will require changes in the current API.
            Array.Resize(ref attributesEncoderIdsOrder, attributesEncoders.Count);
            bool[] isEncoderProcessed = new bool[attributesEncoders.Count];
            int numProcessedEncoders = 0;
            while (numProcessedEncoders < attributesEncoders.Count)
            {
                // Flagged when any of the encoder get processed.
                bool encoderProcessed = false;
                for (int i = 0; i < attributesEncoders.Count; ++i)
                {
                    if (isEncoderProcessed[i])
                        continue; // Encoder already processed.
                    // Check if all parent encoders are already processed.
                    bool canBeProcessed = true;
                    for (int p = 0; p < attributesEncoders[i].NumAttributes; ++p)
                    {
                        int attId = attributesEncoders[i].GetAttributeId(p);
                        for (int ap = 0;
                            ap < attributesEncoders[i].NumParentAttributes(attId);
                            ++ap)
                        {
                            int parentAttId =
                                attributesEncoders[i].GetParentAttributeId(attId, ap);
                            int parentEncoderId = attributeToEncoderMap[parentAttId];
                            if (parentAttId != i && !isEncoderProcessed[parentEncoderId])
                            {
                                canBeProcessed = false;
                                break;
                            }
                        }
                    }
                    if (!canBeProcessed)
                        continue; // Try to process the encoder in the next iteration.
                    // Encoder can be processed. Update the encoding order.
                    attributesEncoderIdsOrder[numProcessedEncoders++] = i;
                    isEncoderProcessed[i] = true;
                    encoderProcessed = true;
                }
                if (!encoderProcessed &&
                    numProcessedEncoders < attributesEncoders.Count)
                {
                    // No encoder was processed but there are still some remaining unprocessed
                    // encoders.
                    throw DracoUtils.Failed();
                }
            }

            // Now for every encoder, reorder the attributes to satisfy their
            // dependencies (an attribute may still depend on other attributes within an
            // encoder).
            int[] attributeEncodingOrder = null;
            bool[] isAttributeProcessed = new bool[pointCloud.NumAttributes];
            int numProcessedAttributes;
            for (int aeOrder = 0; aeOrder < attributesEncoders.Count; ++aeOrder)
            {
                int ae = attributesEncoderIdsOrder[aeOrder];
                int numEncoderAttributes =
                    attributesEncoders[ae].NumAttributes;
                if (numEncoderAttributes < 2)
                    continue; // No need to resolve dependencies for a single attribute.
                numProcessedAttributes = 0;
                Array.Resize(ref attributeEncodingOrder, numEncoderAttributes);
                while (numProcessedAttributes < numEncoderAttributes)
                {
                    // Flagged when any of the attributes get processed.
                    bool attributeProcessed = false;
                    for (int i = 0; i < numEncoderAttributes; ++i)
                    {
                        int attId = attributesEncoders[ae].GetAttributeId(i);
                        if (isAttributeProcessed[i])
                            continue; // Attribute already processed.
                        // Check if all parent attributes are already processed.
                        bool canBeProcessed = true;
                        for (int p = 0;
                            p < attributesEncoders[ae].NumParentAttributes(attId);
                            ++p)
                        {
                            int parentAttId =
                                attributesEncoders[ae].GetParentAttributeId(attId, p);
                            if (!isAttributeProcessed[parentAttId])
                            {
                                canBeProcessed = false;
                                break;
                            }
                        }
                        if (!canBeProcessed)
                            continue; // Try to process the attribute in the next iteration.
                        // Attribute can be processed. Update the encoding order.
                        attributeEncodingOrder[numProcessedAttributes++] = i;
                        isAttributeProcessed[i] = true;
                        attributeProcessed = true;
                    }
                    if (!attributeProcessed &&
                        numProcessedAttributes < numEncoderAttributes)
                    {
                        // No attribute was processed but there are still some remaining
                        // unprocessed attributes.
                        throw DracoUtils.Failed();
                    }
                }
                // Update the order of the attributes within the encoder.
                attributesEncoders[ae].SetAttributeIds(attributeEncodingOrder);
            }
        }

        public PointAttribute GetPortableAttribute(int parent_att_id)
        {
            if (parent_att_id < 0 || parent_att_id >= pointCloud.NumAttributes)
                return null;
            int parent_att_encoder_id = attributeToEncoderMap[parent_att_id];
            return attributesEncoders[parent_att_encoder_id].GetPortableAttribute(parent_att_id);
        }
    }
}
