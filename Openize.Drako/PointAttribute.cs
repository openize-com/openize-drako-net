using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako
{
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    class PointAttribute : GeometryAttribute
    {
        private const int kInvalidAttributeValueIndex = -1;
        private int[] indicesMap;
        private int numUniqueEntries;
        private DataBuffer buffer;

        public DataBuffer Buffer { get { return buffer; } }


        public int NumUniqueEntries
        {
            get { return numUniqueEntries; }
            set { numUniqueEntries = value; }
        }

        public PointAttribute()
        {
            
        }

        public PointAttribute(AttributeType type, DataType dataType, int components, bool normalized, int byteStride = -1, int byteOffset = 0, DataBuffer buffer = null)
        {
            this.AttributeType = type;
            this.ComponentsCount = components;
            this.DataType = dataType;
            this.Normalized = normalized;
            this.ByteOffset = byteOffset;
            this.buffer = buffer;
            if (byteStride == -1)
                this.ByteStride = DracoUtils.DataTypeLength(dataType) * components;
            else
                this.ByteStride = byteStride;
            if(buffer != null)
                numUniqueEntries = buffer.Length / this.ByteStride;
        }

#if !CSPORTER
        /// <summary>
        /// Wrap Vector2 to PointAttribute
        /// </summary>
        /// <param name="type">Attribute's type</param>
        /// <param name="vectors">Attribute data</param>
        /// <returns>Instance of <see cref="PointAttribute"/> wrapped from the vectors.</returns>
        public static PointAttribute Wrap(AttributeType type, Span<Vector2> vectors)
        {
            var bytes = MemoryMarshal.AsBytes(vectors);
            return new PointAttribute(
                type: type,
                dataType: DataType.FLOAT32, components: 2, // data type for position is  float[3]
                normalized: false,
                byteStride: -1, // -1 means auto infer the stride from data type
                byteOffset: 0, //offset in the buffer to the first position
                buffer: new DataBuffer(bytes) //construct a data buffer from Span<byte>, the number of the positions is calculated inside the constructor of PointAttribute
             );
        }
        /// <summary>
        /// Wrap Vector3 to PointAttribute
        /// </summary>
        /// <param name="type">Attribute's type</param>
        /// <param name="vectors">Attribute data</param>
        /// <returns>Instance of <see cref="PointAttribute"/> wrapped from the vectors.</returns>
        public static PointAttribute Wrap(AttributeType type, Span<Vector3> vectors)
        {
            var bytes = MemoryMarshal.AsBytes(vectors);
            return new PointAttribute(
                type: type,
                dataType: DataType.FLOAT32, components: 3, // data type for position is  float[3]
                normalized: false,
                byteStride: -1, // -1 means auto infer the stride from data type
                byteOffset: 0, //offset in the buffer to the first position
                buffer: new DataBuffer(bytes) //construct a data buffer from Span<byte>, the number of the positions is calculated inside the constructor of PointAttribute
             );
        }
        /// <summary>
        /// Wrap Vector3 to PointAttribute
        /// </summary>
        /// <param name="type">Attribute's type</param>
        /// <param name="components">Number of components of Attribute's type</param>
        /// <param name="values">Attribute data</param>
        /// <returns>Instance of <see cref="PointAttribute"/> wrapped from the vectors.</returns>
        public static PointAttribute Wrap(AttributeType type, int components, float[] values)
        {

#if CSPORTER
            var bytes = new byte[4 * values.Length];
            Unsafe.ToByteArray(values, 0, values.Length, bytes, 0);
#else
            var bytes = MemoryMarshal.AsBytes(values.AsSpan());
#endif

            return new PointAttribute(
                type: type,
                dataType: DataType.FLOAT32, components: components, 
                normalized: false,
                byteStride: -1, // -1 means auto infer the stride from data type
                byteOffset: 0, //offset in the buffer to the first position
                buffer: new DataBuffer(bytes) //construct a data buffer from Span<byte>, the number of the positions is calculated inside the constructor of PointAttribute
             );
        }
#else
        /// <summary>
        /// Wrap Vector2 to PointAttribute
        /// </summary>
        /// <param name="type">Attribute's type</param>
        /// <param name="vectors">Attribute data</param>
        /// <returns>Instance of <see cref="PointAttribute"/> wrapped from the vectors.</returns>
        public static PointAttribute Wrap(AttributeType type, Vector2[] vectors)
        {
            var bytes = new byte[4 * 2 * vectors.Length];
            for(int i = 0, p = 0; i < vectors.Length; i++)
            {
                Unsafe.PutLE32(bytes, p, Unsafe.FloatToUInt32((float)vectors[i].X));
                Unsafe.PutLE32(bytes, p + 4, Unsafe.FloatToUInt32((float)vectors[i].Y));
                p += 8;
            }
            return new PointAttribute(
                type,
                DataType.FLOAT32, 2, // data type for position is  float[3]
                false,
                -1, // -1 means auto infer the stride from data type
                0, //offset in the buffer to the first position
                new DataBuffer(bytes) //construct a data buffer from Span<byte>, the number of the positions is calculated inside the constructor of PointAttribute
             );
        }
        /// <summary>
        /// Wrap Vector3 to PointAttribute
        /// </summary>
        /// <param name="type">Attribute's type</param>
        /// <param name="vectors">Attribute data</param>
        /// <returns>Instance of <see cref="PointAttribute"/> wrapped from the vectors.</returns>
        public static PointAttribute Wrap(AttributeType type, Vector3[] vectors)
        {
            var bytes = new byte[4 * 3 * vectors.Length];
            for(int i = 0, p = 0; i < vectors.Length; i++)
            {
                Unsafe.PutLE32(bytes, p, Unsafe.FloatToUInt32((float)vectors[i].X));
                Unsafe.PutLE32(bytes, p + 4, Unsafe.FloatToUInt32((float)vectors[i].Y));
                Unsafe.PutLE32(bytes, p + 8, Unsafe.FloatToUInt32((float)vectors[i].Z));
                p += 12;
            }
            return new PointAttribute(
                type,
                DataType.FLOAT32, 3, // data type for position is  float[3]
                false,
                -1, // -1 means auto infer the stride from data type
                0, //offset in the buffer to the first position
                new DataBuffer(bytes) //construct a data buffer from Span<byte>, the number of the positions is calculated inside the constructor of PointAttribute
             );
        }

#endif

        /// <summary>
        /// Fills outData with the raw value of the requested attribute entry.
        /// outData must be at least byteStride long.
        /// </summary>
        /// <param name="attIndex">Index to the attribute entry</param>
        /// <param name="outData">Byte array to receive the attribute entry value</param>
        public void GetValue(int attIndex, byte[] outData)
        {
            int bytePos = ByteOffset + ByteStride * attIndex;
            buffer.Read(bytePos, outData, ByteStride);
        }

        public int GetBytePos(int attIndex)
        {
            return ByteOffset + ByteStride * attIndex;
        }

        public void GetValue(int attIndex, ushort[] v)
        {
            int bytePos = ByteOffset + ByteStride * attIndex;
            var data = buffer.GetBuffer();
            for (int i = 0; i < v.Length; i++)
            {
                v[i] = Unsafe.GetLE16(data, bytePos);
                bytePos += 2;
            }
        }
        public void GetValue(int attIndex, uint[] v)
        {
            int bytePos = ByteOffset + ByteStride * attIndex;
            var data = buffer.GetBuffer();
            for (int i = 0; i < v.Length; i++)
            {
                v[i] = Unsafe.GetLE32(data, bytePos);
                bytePos += 4;
            }
        }
        public void GetValue(int attIndex, float[] v)
        {
            int bytePos = ByteOffset + ByteStride * attIndex;
            for (int i = 0; i < v.Length; i++)
            {
                v[i] = buffer.ReadFloat(bytePos);
                bytePos += 4;
            }
        }
        public void GetValue(int attIndex, Span<float> v)
        {
            int bytePos = ByteOffset + ByteStride * attIndex;
            for (int i = 0; i < v.Length; i++)
            {
                v[i] = buffer.ReadFloat(bytePos);
                bytePos += 4;
            }
        }

        public Vector3 GetValueAsVector3(int attIndex)
        {
            int bytePos = ByteOffset + ByteStride * attIndex;
            var v = new Vector3();
            v.X = buffer.ReadFloat(bytePos);
            bytePos += 4;
            v.Y = buffer.ReadFloat(bytePos);
            bytePos += 4;
            v.Z = buffer.ReadFloat(bytePos);
            return v;
        }

        /// <summary>
        /// Prepares the attribute storage for the specified number of entries.
        /// </summary>
        /// <param name="numAttributeValues">Number of the attribute entries to preallocate</param>
        /// <returns>true means successed.</returns>
        public void Reset(int numAttributeValues)
        {
            if (buffer == null)
            {
                buffer = new DataBuffer();
            }
            int entrySize = DracoUtils.DataTypeLength(DataType) * ComponentsCount;
            buffer.Capacity = numAttributeValues*entrySize;
            buffer.Length = numAttributeValues*entrySize;
            this.ByteStride = entrySize;
            this.ByteOffset = 0;
            // Assign the new buffer to the parent attribute.
            numUniqueEntries = numAttributeValues;
        }

        public int MappedIndex(int pointIndex)
        {
            if (indicesMap == null)
                return pointIndex;
            return indicesMap[pointIndex];
        }

        /// <summary>
        /// Sets the new number of unique attribute entries for the attribute.
        /// </summary>
        /// <param name="newNumUniqueEntries">Number of unique entries</param>
        public void Resize(int newNumUniqueEntries)
        {
            numUniqueEntries = newNumUniqueEntries;
        }

        /// <summary>
        /// Functions for setting the type of mapping between point indices and
        /// attribute entry ids.
        /// This function sets the mapping to implicit, where point indices are equal
        /// to attribute entry indices.
        /// </summary>
        public bool IdentityMapping
        {
            get { return indicesMap == null; }
            set
            {
                if (value)
                    indicesMap = null;
                else if (indicesMap == null)
                    indicesMap = new int[0];
            }
        }

        public int[] IndicesMap { get { return indicesMap; } }
        internal AttributeTransformData AttributeTransformData { get; set; }

        /// <summary>
        /// This function sets the mapping to be explicitly using the indicesMap
        /// array that needs to be initialized by the caller.
        /// </summary>
        /// <param name="numPoints">Number of actual points</param>
        public void SetExplicitMapping(int numPoints)
        {
            var fillStart = 0;
            if (indicesMap == null)
                indicesMap = new int[numPoints];
            else
            {
                fillStart = indicesMap.Length;
                Array.Resize(ref indicesMap, numPoints);
            }
            for (int i = fillStart; i < numPoints; i++)
                indicesMap[i] = kInvalidAttributeValueIndex;

        }

        /// <summary>
        /// Set an explicit map entry for a specific point index.
        /// </summary>
        /// <param name="pointIndex">Index to the point</param>
        /// <param name="entryIndex">Index to the attribute entry</param>
        public void SetPointMapEntry(int pointIndex,
            int entryIndex)
        {
            indicesMap[pointIndex] = entryIndex;
        }

        public void ConvertValue(int attId, out int i)
        {
            int pos = ByteOffset + ByteStride * attId;
            i = buffer.ReadInt(pos);
        }
        internal LongVector3 ConvertValue(int attId)
        {
            int pos = ByteOffset + ByteStride * attId;
            if (DataType == DataType.INT32 || DataType == DataType.UINT32)
            {
                long x = buffer.ReadInt(pos);
                long y = buffer.ReadInt(pos + 4);
                long z = buffer.ReadInt(pos + 8);
                return new LongVector3(x, y, z);
            }
            else
                throw new NotImplementedException("Unsupported type cast");
        }

        struct ValueKey : IComparable<ValueKey>
        {
            private byte[] data;
            private int offset;
            private int size;
            private readonly long hashCode;

            public int Offset { get { return offset; } }

            public ValueKey(PointAttribute attribute, byte[] data, int index)
            {
                this.size = attribute.ByteStride;
                this.offset = attribute.ByteOffset + attribute.ByteStride * index;
                this.data = data;
                hashCode = DracoUtils.GetHashCode(data, offset, size);
            }

            public override string ToString()
            {
                StringBuilder sb = new StringBuilder();
                sb.AppendFormat("{0}, ", BitConverter.ToSingle(data, offset));
                sb.AppendFormat("{0}, ", BitConverter.ToSingle(data, offset + 4));
                sb.AppendFormat("{0}", BitConverter.ToSingle(data, offset + 8));

                return sb.ToString();
            }

            public override int GetHashCode()
            {
                return (int)hashCode;
            }

            public override bool Equals(object obj)
            {
                if (!(obj is ValueKey))
                    return false;
                return CompareTo((ValueKey) obj) == 0;
            }

            public int CompareTo(ValueKey rhs)
            {
                if (rhs.hashCode != hashCode)
                    return -1;
                if (rhs.offset == offset)
                    return -1;
                int ret = DracoUtils.Compare(data, offset, data, rhs.offset, size);
                return ret;
            }
        }

        public void DeduplicateValues()
        {
            int offset = ByteOffset;
            int stride = ByteStride;
            Dictionary<ValueKey, int> indiceMap = new Dictionary<ValueKey, int>();
            Dictionary<int, int> valueMap = new Dictionary<int, int>();
            int uniqueValues = 0;
            byte[] tmp = buffer.ToArray();
            for (int i = 0; i < this.NumUniqueEntries; i++, offset += stride)
            {
                ValueKey k = new ValueKey(this, tmp, i);
                int idx;
                if (!indiceMap.TryGetValue(k, out idx))
                {
                    SetAttributeValue(uniqueValues, Buffer.GetBuffer(), k.Offset);
                    idx = uniqueValues++;
                    indiceMap.Add(k, idx);
                }
                valueMap[i] = idx;
            }
            if (uniqueValues == numUniqueEntries)
                return ;//cannot deduplicate values
            if (IdentityMapping)
            {
                SetExplicitMapping(numUniqueEntries);
                for (int i = 0; i < numUniqueEntries; i++)
                {
                    SetPointMapEntry(i, valueMap[i]);
                }
            }
            else
            {
                for (int i = 0; i < this.indicesMap.Length; i++)
                {
                    SetPointMapEntry(i, valueMap[this.indicesMap[i]]);
                }
            }
            numUniqueEntries = uniqueValues;


        }

        /// <summary>
        /// Copy raw bytes from buffer with given offset to the attribute's internal buffer at specified element index
        /// </summary>
        /// <param name="index">Index of the attribute entry</param>
        /// <param name="buffer">Source raw bytes to copy from</param>
        /// <param name="offset">Offset to the buffer</param>
        public void SetAttributeValue(int index, byte[] buffer, int offset)
        {
            int dstOffset = ByteOffset + ByteStride * index;
            byte[] dst = this.buffer.GetBuffer();
            Array.Copy(buffer, offset, dst, dstOffset, ByteStride);
        }

        internal void SetAttributeValue(int index, uint[] vals)
        {
            int offset = ByteOffset + ByteStride * index;
            byte[] dst = this.buffer.GetBuffer();
            for (int i = 0; i < vals.Length; i++)
            {
                Unsafe.PutLE32(dst, offset, vals[i]);
                offset += 4;
            }
        }
        internal void SetAttributeValue(int index, ushort[] vals)
        {
            int offset = ByteOffset + ByteStride * index;
            byte[] dst = this.buffer.GetBuffer();
            for (int i = 0; i < vals.Length; i++)
            {
                Unsafe.PutLE16(dst, offset, vals[i]);
                offset += 2;
            }
        }


        public override void CopyFrom(GeometryAttribute attr)
        {

            if (buffer == null)
            {
                // If the destination attribute doesn't have a valid buffer, create it.
                buffer = new DataBuffer();
                ByteStride = 0;
                ByteOffset = 0;
            }

            base.CopyFrom(attr);
            var pa = (PointAttribute)attr;
            numUniqueEntries = pa.numUniqueEntries;
            if(pa.indicesMap != null) 
                indicesMap = (int[]) pa.indicesMap.Clone();
            /*
            if (pa.attributeTransformData)
            {
                attributeTransformData = pa.attributeTransformData.Clone();
            }
            else
            {
                attributeTransformData = null;
            }
            */
        }

        internal Span<byte> GetAddress(int attIndex)
        {
            var byte_pos = GetBytePos(attIndex);
            return buffer.AsSpan().Slice(byte_pos);
        }
    }
}
