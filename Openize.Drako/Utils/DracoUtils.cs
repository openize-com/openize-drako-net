using System;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;

namespace Openize.Drako.Utils
{
    class DracoUtils
    {
        [Porter("intrinsic")]
        [Conditional("CSPORTER")]
        public static void __EmitJavaCode(String str)
        {
        }

        public static void Fill(int[] array, int value)
        {
            for (int i = 0; i < array.Length; i++)
                array[i] = value;
        }
        /// <summary>
        /// Returns the point id of |c| without using a corner table.
        /// </summary>
        public static int CornerToPointId(int c, DracoMesh mesh)
        {
            return mesh.ReadCorner(c);// (c/3)[c%3];
        }

        public static int IncrementMod(int I, int M)
        {
            return (((I) == ((M) - 1)) ? 0 : ((I) + 1));
        }
        public static int DataTypeLength(DataType dt)
        {
            switch (dt)
            {
                case DataType.INT8:
                case DataType.UINT8:
                    return 1;
                case DataType.INT16:
                case DataType.UINT16:
                    return 2;
                case DataType.INT32:
                case DataType.UINT32:
                    return 4;
                case DataType.INT64:
                case DataType.UINT64:
                    return 8;
                case DataType.FLOAT32:
                    return 4;
                case DataType.FLOAT64:
                    return 8;
                case DataType.BOOL:
                    return 1;
                default:
                    return -1;
            }
        }




        /// <summary>
        /// Copies the |bnbits| from the src integer into the |dst| integer using the
        /// provided bit offsets |dst_offset| and |src_offset|.
        /// </summary>
        public static void CopyBits32(ref uint dst, int dst_offset, uint src, int src_offset, int nbits)
        {
            uint mask = (~0u) >> (32 - nbits) << dst_offset;
            dst = (dst & (~mask)) | (((src >> src_offset) << dst_offset) & mask);
        }

        /// <summary>
        /// Returns the number of '1' bits within the input 32 bit integer.
        /// </summary>
        public static int CountOnes32(uint n)
        {
            n -= ((n >> 1) & 0x55555555);
            n = ((n >> 2) & 0x33333333) + (n & 0x33333333);
            return (int)((((n + (n >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24);
        }

        public static uint ReverseBits32(uint n)
        {
            n = ((n >> 1) & 0x55555555) | ((n & 0x55555555) << 1);
            n = ((n >> 2) & 0x33333333) | ((n & 0x33333333) << 2);
            n = ((n >> 4) & 0x0F0F0F0F) | ((n & 0x0F0F0F0F) << 4);
            n = ((n >> 8) & 0x00FF00FF) | ((n & 0x00FF00FF) << 8);
            return (n >> 16) | (n << 16);
        }

        public static bool IsIntegerType(DataType type)
        {
            switch (type)
            {
                case DataType.INT8:
                case DataType.UINT8:
                case DataType.INT16:
                case DataType.UINT16:
                case DataType.INT32:
                case DataType.UINT32:
                    return true;
                default:
                    return false;
            }
        }

        /// <summary>
        /// Returns the most location of the most significant bit in the input integer
        /// |n|.
        /// The funcionality is not defined for |n == 0|.
        /// </summary>
        public static int MostSignificantBit(uint n)
        {
            int msb = -1;
            while (n != 0)
            {
                msb++;
                n >>= 1;
            }
            return msb;
        }

        public static LongVector3 Add(LongVector3 a, LongVector3 b)
        {
            return new LongVector3(a.x + b.x, a.y + b.y, a.z + b.z);
        }
        public static LongVector3 Sub(LongVector3 a, LongVector3 b)
        {
            return new LongVector3(a.x - b.x, a.y - b.y, a.z - b.z);
        }
        public static LongVector3 Div(LongVector3 a, long b)
        {
            return new LongVector3(a.x / b, a.y / b, a.z / b);
        }
        public static LongVector3 Mul(LongVector3 a, long b)
        {
            return new LongVector3(a.x * b, a.y * b, a.z * b);
        }
        public static uint SquaredNorm(LongVector3 a)
        {
            return (uint)Dot(a, a);
        }
        public static long Dot(LongVector3 a, LongVector3 b)
        {
            long ret = a.x * b.x + a.y * b.y + a.z * b.z;
            return ret;

        }

        public static long AbsSum(Span<int> v)
        {
            long ret = 0;
            for (int i = 0; i < v.Length; i++)
                ret += Math.Abs(v[i]);
            return ret;
        }
        public static long AbsSum(LongVector3 v)
        {
            long ret = Math.Abs(v.x) + Math.Abs(v.y) + Math.Abs(v.z);
            return ret;
        }

        public static LongVector3 CrossProduct(LongVector3 u, LongVector3 v)
        {
            LongVector3 r = new LongVector3();
            r.x = (u.y * v.z) - (u.z * v.y);
            r.y = (u.z * v.x) - (u.x * v.z);
            r.z = (u.x * v.y) - (u.y * v.x);
            return r;
        }

        public static DrakoException Failed()
        {
            return new DrakoException();
        }

        internal static ulong IntSqrt(ulong number)
        {
            if (number == 0)
                return 0;
            ulong actNumber = number;
            ulong squareRoot = 1;
            while (actNumber >= 2)
            {
                squareRoot *= 2;
                actNumber /= 4;
            }
            // Perform Newton's (or Babylonian) method to find the true floor(sqrt()).
            do
            {
                squareRoot = (squareRoot + number / squareRoot) / 2;
            } while (squareRoot * squareRoot > number);
            return squareRoot;
        }

        public static bool VecEquals(LongVector3 u, LongVector3 v)
        {
            return u.x == v.x && u.y == v.y && u.z == v.z;
        }

        internal static long GetHashCode(byte[] bytes, int offset, int size)
        {
            //FNV hash algorithm
            unchecked
            {
                const int p = 16777619;
                long hash =  2166136261;

                int end = offset + size;
                for (int i = offset; i < end; i++)
                {
                    hash = (hash ^ bytes[i]) * p;
                }

                hash += hash << 13;
                hash ^= hash >> 7;
                hash += hash << 3;
                hash ^= hash >> 17;
                hash += hash << 5;

                return hash;
            }
        }

#if CSPORTER
        internal static unsafe int Compare(byte[] data1, byte[] data2, int size)
        {
            return Compare(data1, 0, data2, 0, size);
        }
        internal static unsafe int Compare(byte[] data1, int offset1, byte[] data2, int offset2, int size)
        {
            //TODO: range may out of bounds
            for(int i1 = offset1, i2 = offset2, i = 0; i < size; i++, i1++,i2++)
            {
                var v = data1[i1].CompareTo(data2[i2]);
                if (v != 0)
                    return v;
            }
            return 0;
        }

#else
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static unsafe int Compare(byte[] a, byte[] b, int length)
        {
            fixed(byte* pa = a, pb = b)
            {
                return UlongCompareInternal(pa, pb, length);
            }
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static unsafe int Compare(byte[] data1, int offset1, byte[] data2, int offset2, int size)
        {
            fixed(byte* pa = data1, pb = data2)
            {
                return UlongCompareInternal(pa + offset1, pb + offset2, size);
            }
        }

        private static unsafe int UlongCompareInternal(byte* xPtr, byte* yPtr, int length)
        {
            byte* lastAddr = xPtr + length;
            byte* lastAddrMinus32 = lastAddr - 32;
            while (xPtr < lastAddrMinus32) // Compare 32 bytes by loop
            {
                // Compare leading 64 bits
                var v = *(ulong*)xPtr - *(ulong*)yPtr;
                if (v != 0UL) return (int)v;
                // Second ulong
                v = *(ulong*)(xPtr + 8) - *(ulong*)(yPtr + 8);
                if (v != 0) return (int)v;
                // Third ulong
                v = *(ulong*)(xPtr + 16) - *(ulong*)(yPtr + 16);
                if (v != 0) return (int)v;
                // Fourth ulong
                v = *(ulong*)(xPtr + 24) - *(ulong*)(yPtr + 24);
                if (v != 0) return (int)v;
                // 256 bits compared, move to next 256 bits
                xPtr += 32;
                yPtr += 32;
            }
            // Compare rest bytes one by one
            while (xPtr < lastAddr)
            {
                if (*xPtr != *yPtr) return *xPtr - *yPtr;
                xPtr++;
                yPtr++;
            }
            return 0;
        }
#endif

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static bool IsZero(float v)
        {
            return v == 0;
        }

    }
}
