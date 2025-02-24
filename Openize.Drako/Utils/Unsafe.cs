using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace FileFormat.Drako.Utils
{
    unsafe class Unsafe
    {

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static float GetFloat(byte[] arr, int off)
        {
#if CSPORTER
            uint v = GetLE32(arr, off);
            DracoUtils.__EmitJavaCode("return java.lang.Float.intBitsToFloat(v);");
#else
            fixed(byte*p = arr)
            {
                return *(float*)(p + off);
            }
#endif
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ushort GetLE16(byte[] arr, int off = 0)
        {
            if (BitConverter.IsLittleEndian)
            {
#if CSPORTER
                ushort val = (ushort)(arr[off + 1] << 8);
                val |= arr[off];
                return val;
#else
                fixed (byte* p = arr)
                {
                    return *(ushort*)(p + off);
                }
#endif
            }
            else
            {
                ushort val = (ushort)(arr[off] << 8);
                val |= arr[off + 1];
                return val;
            }
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static uint GetLE24(byte[] arr, int off)
        {
            uint val = (uint) (arr[off + 2] << 16);
            val |= (uint) (arr[off + 1] << 8);
            val |= arr[off + 0];
            return val;
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static unsafe uint GetLE32(byte[] arr, int off)
        {
            if(BitConverter.IsLittleEndian)
            {
#if CSPORTER
                    uint val = (uint)(arr[off + 0] << 0);
                    val |= (uint)(arr[off + 1] << 8);
                    val |= (uint)(arr[off + 2] << 16);
                    val |= arr[off + 3] << 24;
                    return val;
#else
                fixed (byte* p = arr)
                {
                    return *(uint*)(p + off);
                }
#endif
            }
            else if (off + 3 < arr.Length)
            {
                unchecked
                {
                    uint val = (uint)(arr[off + 0] << 24);
                    val |= (uint)(arr[off + 1] << 16);
                    val |= (uint)(arr[off + 2] << 8);
                    val |= arr[off + 3];
                    return val;
                }
            }
            else
                return 0;
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ulong GetLE64(byte[] arr, int off)
        {
            if (BitConverter.IsLittleEndian)
            {
#if CSPORTER
                long val = 0;
                val |= ((long)arr[off + 7] << 56);
                val |= ((long)arr[off + 6] << 48);
                val |= ((long)arr[off + 5] << 40);
                val |= ((long)arr[off + 4] << 32);
                val |= ((long)arr[off + 3] << 24);
                val |= ((long)arr[off + 2] << 16);
                val |= ((long)arr[off + 1] << 8);
                val |= ((long)arr[off + 0] << 0);
                return val;
#else
                fixed (byte* p = arr)
                {
                    return *(ulong*)(p + off);
                }
#endif
            }
            else
            {
#pragma warning disable 0675
                long val = 0;
                val |= ((long)arr[off + 0] << 56);
                val |= ((long)arr[off + 1] << 48);
                val |= ((long)arr[off + 2] << 40);
                val |= ((long)arr[off + 3] << 32);
                val |= ((long)arr[off + 4] << 24);
                val |= ((long)arr[off + 5] << 16);
                val |= ((long)arr[off + 6] << 8);
#pragma warning restore 0675
                val |= arr[off + 7];
                return (ulong)val;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void PutLE16(byte[] arr, int off, ushort val)
        {
            if (BitConverter.IsLittleEndian)
            {
#if CSPORTER
                arr[off + 1] = (byte)((val >> 8) & 0xff);
                arr[off + 0] = (byte)((val >> 0) & 0xff);
#else
                fixed (byte* p = arr)
                {
                    *((ushort*)(p + off)) = val;
                }
#endif
            }
            else
            {
                arr[off] = (byte)((val >> 8) & 0xff);
                arr[off + 1] = (byte)((val >> 0) & 0xff);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void PutLE24(byte[] arr, int off, uint val)
        {
            arr[off] = (byte)((val >> 0) & 0xff);
            arr[off + 1] = (byte)((val >> 8) & 0xff);
            arr[off + 2] = (byte)((val >> 16) & 0xff);
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void PutLE32(byte[] arr, int off, uint val)
        {
            if (BitConverter.IsLittleEndian)
            {
#if CSPORTER
                arr[off + 0] = (byte)((val >> 0) & 0xff);
                arr[off + 1] = (byte)((val >> 8) & 0xff);
                arr[off + 2] = (byte)((val >> 16) & 0xff);
                arr[off + 3] = (byte)((val >> 24) & 0xff);
#else
                fixed (byte* p = arr)
                {
                    *((uint*)(p + off)) = val;
                }
#endif
            }
            else
            {

                arr[off + 3] = (byte)((val >> 0) & 0xff);
                arr[off + 2] = (byte)((val >> 8) & 0xff);
                arr[off + 1] = (byte)((val >> 16) & 0xff);
                arr[off + 0] = (byte)((val >> 24) & 0xff);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static uint FloatToUInt32(float value)
        {
#if CSPORTER
            DracoUtils.__EmitJavaCode("return java.lang.Float.floatToIntBits(value);");
#else
            return *((uint*) &value);
#endif
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void ToByteArray(float[] arr, int start, int len, byte[] dst, int offset)
        {
#if CSPORTER
            for(int i = 0, psrc = start, pdst = offset; i < len; psrc++, pdst += 4, i++)
            {
                var v = Unsafe.FloatToUInt32(arr[psrc]);
                PutLE32(dst, pdst, v);
            }
#else
            fixed (float* p = arr)
            {
                IntPtr pstart = new IntPtr((long) p + (start << 2));
                Marshal.Copy(pstart, dst, offset, len << 2);
            }
#endif
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static float[] ToFloatArray(byte[] array, float[] ret)
        {
            if (ret.Length == 0)
                return ret;
#if CSPORTER
            for(int i = 0, d = 0; i < ret.Length; i += 4, d++)
            {
                ret[d] = Unsafe.GetFloat(array, i);
            }
#else
            fixed (float* p = ret)
            {
                Marshal.Copy(array, 0, new IntPtr(p), ret.Length * sizeof(float));
            }
#endif
            return ret;
        }
    }
}
