using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako
{

    /// <summary>
    /// DracoPointCloud is a collection of n-dimensional points that are described by a
    /// set of PointAttributes that can represent data such as positions or colors
    /// of individual points (see pointAttribute.h).
    /// </summary>
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    class DracoPointCloud
    {

        /// <summary>
        /// Attributes describing the point cloud.
        /// </summary>
        private List<PointAttribute> attributes = new List<PointAttribute>();

        private List<GeometryMetadata> metadatas = new List<GeometryMetadata>();

        /// <summary>
        /// Ids of named attributes of the given type.
        /// </summary>
        private IntList[] namedAttributeIndex = new IntList[(int) AttributeType.NamedAttributesCount];

        /// <summary>
        /// The number of n-dimensional points. All point attribute values are stored
        /// in corresponding PointAttribute instances in the |attributes| array.
        /// </summary>
        private int numPoints;
        /// <summary>
        /// Returns the number of named attributes of a given type.
        /// </summary>
        /// <param name="type">What type of attribute to count, defined by <see cref="AttributeType"/></param>
        /// <returns>The number of attributes in given type.</returns>
        public int NumNamedAttributes(AttributeType type)
        {
            int idx = (int) type;
            if (idx < 0 || idx >= namedAttributeIndex.Length)
                return 0;
            return namedAttributeIndex[idx].Count;
        }

        /// <summary>
        /// Returns attribute id of the first named attribute with a given type or -1
        /// when the attribute is not used by the point cloud.
        /// </summary>
        ///<param name="type">What type of attribute to find, defined by <see cref="AttributeType"/></param>
        /// <returns>The id of the first named attribute in given type</returns>
        public int GetNamedAttributeId(AttributeType type)
        {
            return GetNamedAttributeId(type, 0);
        }

        /// <summary>
        /// Returns the id of the i-th named attribute of a given type.
        /// </summary>
        ///<param name="type">What type of attribute to find, defined by <see cref="AttributeType"/></param>
        /// <param name="i">Index of the attribute.</param>
        /// <returns>The id of the first named attribute in given type</returns>
        public int GetNamedAttributeId(AttributeType type, int i)
        {
            int idx = (int) type;
            if (idx < 0 || idx >= namedAttributeIndex.Length)
                return -1;
            var attrs = namedAttributeIndex[idx];
            if (attrs == null)
                return -1;
            if (i < 0 || i >= attrs.Count)
                return -1;
            return attrs[i];
        }

        /// <summary>
        /// Returns the i-th named attribute of a given type.
        /// </summary>
        /// <param name="type">What type of attribute to find, defined by <see cref="AttributeType"/></param>
        /// <param name="i">Index of the attribute.</param>
        /// <returns>The instance of the i-th named attribute in given type</returns>
        public PointAttribute GetNamedAttribute(AttributeType type, int i = 0)
        {
            var id = GetNamedAttributeId(type, i);
            if (id == -1)
                return null;
            return attributes[id];
        }

        /// <summary>
        /// Returns the named attribute of a given custom id.
        /// </summary>
        /// <param name="type">What type of attribute to find, defined by <see cref="AttributeType"/></param>
        /// <param name="customId">Custom id of the attribute.</param>
        /// <returns>The instance of named attribute with given custom id</returns>
        public PointAttribute GetNamedAttributeByCustomId(
            AttributeType type, ushort customId)
        {
            int idx = (int) type;
            for (int attId = 0;
                attId < namedAttributeIndex[idx].Count;
                ++attId)
            {
                if (attributes[namedAttributeIndex[idx][attId]].UniqueId ==
                    customId)
                    return attributes[namedAttributeIndex[idx][attId]];
            }
            return null;
        }

        public int NumAttributes
        {
            get { return attributes.Count; }
        }

        public PointAttribute Attribute(int attId)
        {
            return attributes[attId];
        }

        /// <summary>
        /// Adds a new attribute to the point cloud.
        /// Returns the attribute id.
        /// </summary>
        /// <param name="pa">The new attribute to be added</param>
        /// <returns>The index of the new attribute.</returns>
        public virtual int AddAttribute(PointAttribute pa)
        {
            attributes.Add(pa);
            var attrs = namedAttributeIndex[(int) pa.AttributeType];
            if (attrs == null)
                attrs = namedAttributeIndex[(int) pa.AttributeType] = new IntList();
            int ret = attributes.Count - 1;
            attrs.Add(ret);
            pa.UniqueId = (ushort)ret;
            return ret;
        }


        /// <summary>
        /// Creates and adds a new attribute to the point cloud. The attribute has
        /// properties derived from the provided GeometryAttribute |att|.
        /// If |identityMapping| is set to true, the attribute will use identity
        /// mapping between point indices and attribute value indices (i.e., each point
        /// has a unique attribute value).
        /// If |identityMapping| is false, the mapping between point indices and
        /// attribute value indices is set to explicit, and it needs to be initialized
        /// manually using the PointAttribute::SetPointMapEntry() method.
        /// |numAttributeValues| can be used to specify the number of attribute
        /// values that are going to be stored in the newly created attribute.
        /// Returns attribute id of the newly created attribute.
        /// </summary>
        /// <param name="att">The instance of attribte to be added</param>
        /// <param name="identityMapping">Whether to use identity mapping between point indices and attribute value indices.</param>
        /// <param name="numAttributeValues">Specify the number of attribute values that will be stored in the new attribute.</param>
        /// <returns>The index of the new attribute.</returns>
        public int AddAttribute(GeometryAttribute att, bool identityMapping, int numAttributeValues)
        {

            AttributeType type = att.AttributeType;
            if (type == AttributeType.Invalid)
                return -1;
            int attId = AddAttribute((PointAttribute)att);
            PointAttribute pa = Attribute(attId);
            // Initialize point cloud specific attribute data.
            pa.IdentityMapping = identityMapping;
            if (!identityMapping)
            {
                // First create mapping between indices.
                pa.SetExplicitMapping(numPoints);
            }
            else
            {
                pa.Resize(numPoints);
            }
            if (numAttributeValues > 0)
            {
                pa.Reset(numAttributeValues);
            }
            return attId;
        }

        /// <summary>
        /// Deduplicates all attribute values (all attribute entries with the same
        /// value are merged into a single entry).
        /// </summary>
        /// <returns>true if deduplication successed.</returns>
        public virtual void DeduplicateAttributeValues()
        {

            if (numPoints == 0)
                return; // Unexcpected attribute size.
            // Deduplicate all attributes.
            for (int i = 0; i < attributes.Count; i++)
            {
                var attr = attributes[i];
                attr.DeduplicateValues();
            }
        }

        struct VertexIndex : IComparable<VertexIndex>
        {
            private readonly DracoPointCloud cloud;
            private readonly int p;

            public VertexIndex(DracoPointCloud cloud, int p)
            {
                this.cloud = cloud;
                this.p = p;
            }

            public override int GetHashCode()
            {
                int hash = 0;
                for (int i = 0; i < cloud.attributes.Count; ++i)
                {
                    int attId = cloud.attributes[i].MappedIndex(p);
                    hash = (attId.GetHashCode() << 2) ^ (hash.GetHashCode() << 1);
                }
                return hash;
            }

            public override bool Equals(object obj)
            {
                if (!(obj is VertexIndex))
                    return false;
                int ret = CompareTo((VertexIndex) obj);
                return ret == 0;
            }

            public int CompareTo(VertexIndex other)
            {
                for (int i = 0; i < cloud.attributes.Count; ++i)
                {
                    var attr = cloud.attributes[i];
                    var id0 = attr.MappedIndex(p);
                    var id1 = attr.MappedIndex(other.p);
                    if (id0 < id1)
                        return -1;
                    if (id0 > id1)
                        return 1;
                }
                return 0;
            }

            public override string ToString()
            {
                StringBuilder sb = new StringBuilder();
                sb.AppendFormat("p: ");
                for (int i = 0; i < cloud.attributes.Count; ++i)
                {
                    var attr = cloud.attributes[i];
                    var id0 = attr.MappedIndex(p);
                    if (i > 0)
                        sb.Append("/");
                    sb.Append(id0);
                }
                return sb.ToString();
            }
        }


        /// <summary>
        /// Removes duplicate point ids (two point ids are duplicate when all of their
        /// attributes are mapped to the same entry ids).
        /// </summary>
        public virtual void DeduplicatePointIds()
        {
            // Comparison function between two vertices.
            //TODO: benchmark class/struct here, since struct reduces GC but class reduce memory exchange.
            Dictionary<VertexIndex, int> uniquePointMap = new Dictionary<VertexIndex, int>();
            int numUniquePoints = 0;
            int[] indexMap = new int[numPoints];
            IntList uniquePoints = new IntList();
            // Go through all vertices and find their duplicates.
            for (int i = 0; i < numPoints; ++i)
            {
                int tmp;
                var p = new VertexIndex(this, i);
                if (uniquePointMap.TryGetValue(p, out tmp))
                {
                    indexMap[i] = tmp;
                }
                else
                {
                    uniquePointMap.Add(p, numUniquePoints);
                    indexMap[i] = numUniquePoints++;
                    uniquePoints.Add(i);
                }
            }
            if (numUniquePoints == numPoints)
                return; // All vertices are already unique.

            ApplyPointIdDeduplication(indexMap, uniquePoints);
            NumPoints = numUniquePoints;
        }


        /// <summary>
        /// Gets or sets the number of n-dimensional points stored within the point cloud.
        /// </summary>
        public int NumPoints
        {
            get { return numPoints; }
            set { numPoints = value; }
        }

        internal List<GeometryMetadata> Metadatas => metadatas;

        /// <summary>
        /// Applies id mapping of deduplicated points (called by DeduplicatePointIds).
        /// </summary>
        internal virtual void ApplyPointIdDeduplication(
            int[] idMap,
            IntList uniquePointIds)
        {
            int numUniquePoints = 0;
            for (int i = 0; i < uniquePointIds.Count; i++)
            {
                int newPointId = idMap[uniquePointIds[i]];
                if (newPointId >= numUniquePoints)
                {
                    // New unique vertex reached. Copy attribute indices to the proper
                    // position.
                    for (int a = 0; a < attributes.Count; ++a)
                    {
                        var attr = attributes[a];
                        attr.SetPointMapEntry(newPointId,
                            attr.MappedIndex(uniquePointIds[i]));
                    }
                    numUniquePoints = newPointId + 1;
                }
            }
            for (int a = 0; a < attributes.Count; ++a)
            {
                attributes[a].SetExplicitMapping(numUniquePoints);
            }
        }

        internal PointAttribute GetAttributeByUniqueId(int uniqueId)
        {
            for(int i = 0; i < attributes.Count; ++i)
            {
                var attr = attributes[i];
                if (attr.UniqueId == uniqueId)
                    return attr;
            }
            return null;
        }
    }
}
