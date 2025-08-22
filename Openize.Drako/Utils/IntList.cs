using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako.Utils
{
    /// <summary>
    /// Used to simulate std::vector&lt;int&gt;
    /// here we don't use List&lt;int&gt; because we need to simulate some behavior from C++ version.
    /// </summary>
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
    class IntList
    {
        private int count;
        internal int[] data;

        public int Count => count;

        public IntList()
        {
            data = new int[10];
        }
        public IntList(int size)
        {
            data = new int[size];
            count = size;
        }
        public int Capacity
        {
            get => data.Length;
            set
            {
                if (value <= data.Length)
                    return;
                var newData = new int[value];
                data.CopyTo(newData, 0);
                data = newData;
            }
        }
        public void Add(int v)
        {
            EnsureCapacity(count + 1);
            data[count++] = v;
        }
        public void AddRange(int[] v)
        {
            AddRange(v, v.Length);
        }
        public void AddRange(IntList v)
        {
            AddRange(v.data, v.count);
        }
        public void AddRange(int[] v, int length)
        {
            AddRange(v.AsSpan(0, length));
        }
        public void AddRange(ReadOnlySpan<int> v)
        {
            if (v.Length > 0)
            {
                EnsureCapacity(count + v.Length);
                v.CopyTo(data.AsSpan(count, v.Length));
                count += v.Length;
            }
        }
        public void Clear()
        {
            count = 0;
        }
        public void RemoveAt(int idx)
        {
            if (idx == count - 1)
                count--;
            else
            {
                Array.Copy(data, idx + 1, data, idx, count - idx - 1);
                count--;
            }
        }

        public int this[int idx]
        {
            get { return data[idx]; }
            set { data[idx] = value; }
        }
        public int Back => data[count - 1];
        public void PopBack()
        {
            if(count > 0)
                count--;
        }

        public void Reverse()
        {
            Array.Reverse(data, 0, count);
        }

        public void Resize(int newSize)
        {
            EnsureCapacity(newSize);
            count = newSize;
        }
        public void Resize(int newSize, int newValue)
        {
            EnsureCapacity(newSize);
            for(int i = count;  i < newSize; i++)
                data[i] = newValue;
            count = newSize;
        }
        private void EnsureCapacity(int newSize)
        {
            if (newSize < data.Length)
                return;
            var capacity = data.Length;
            while (capacity < newSize)
                capacity += capacity >> 1;
            Capacity = capacity;
        }
        public int[] ToArray()
        {
            var ret = new int[count];
            Array.Copy(data, 0, ret, 0, count);
            return ret;
        }


    }
}
