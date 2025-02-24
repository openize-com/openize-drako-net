using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace FileFormat.Drako
{

    /// <summary>
    /// Supported attribute types.
    /// </summary>
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    enum AttributeType
    {
        Invalid = -1,
        /// <summary>
        /// Named attributes start here. The difference between named and generic
        /// attributes is that for named attributes we know their purpose and we
        /// can apply some special methods when dealing with them (e.g. during
        /// encoding).
        /// </summary>
        Position = 0,
        Normal,
        Color,
        TexCoord,
        /// <summary>
        /// A special id used to mark attributes that are not assigned to any known
        /// predefined use case. Such attributes are often used for a shader specific
        /// data.
        /// </summary>
        Generic,
        /// <summary>
        /// Total number of different attribute types.
        /// Always keep behind all named attributes.
        /// </summary>
        NamedAttributesCount
    }
    /// <summary>
    /// The class provides access to a specific attribute which is stored in a
    /// DataBuffer, such as normals or coordinates. However, the GeometryAttribute
    /// class does not own the buffer and the buffer itself may store other data
    /// unrelated to this attribute (such as data for other attributes in which case
    /// we can have multiple GeometryAttributes accessing one buffer). Typically,
    /// all attributes for a point (or corner, face) are stored in one block, which
    /// is advantageous in terms of memory access. The length of the entire block is
    /// given by the byteStride, the position where the attribute starts is given by
    /// the byteOffset, the actual number of bytes that the attribute occupies is
    /// given by the dataType and the number of components.
    /// </summary>
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    class GeometryAttribute
    {
        public int ComponentsCount { get; set; }
        public DataType DataType { get; set; }
        public bool Normalized { get; set; }
        public int ByteStride { get; set; }
        public int ByteOffset { get; set; }

        public AttributeType AttributeType { get; set; }

        /// <summary>
        /// User defined 16-bit id that can be for example used to identify generic
        /// attributes. By default |customId| == 0.
        /// </summary>
        public ushort UniqueId { get; set; }


        public override string ToString()
        {
            return string.Format("#{3} {0} : {1}[{2}]", AttributeType, DataType, ComponentsCount, UniqueId);
        }

        public virtual void CopyFrom(GeometryAttribute attr)
        {
            ComponentsCount = attr.ComponentsCount;
            DataType = attr.DataType;
            Normalized = attr.Normalized;
            ByteStride = attr.ByteStride;
            ByteOffset = attr.ByteOffset;
            AttributeType = attr.AttributeType;
            UniqueId = attr.UniqueId;
        }
    }
}
