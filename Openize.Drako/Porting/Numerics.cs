using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FileFormat.Drako
{
#if CSPORTER && !DRACO_EMBED_MODE
    public struct Vector2
    {
        public float X;
        public float Y;
        public Vector2(float x, float y)
        {
            this.X = x;
            this.Y = y;
        }

        public static bool operator ==(Vector2 a, Vector2 b)
        {
            return a.X == b.X && a.Y == b.Y;
        }
        public static bool operator !=(Vector2 a, Vector2 b)
        {
            return a.X != b.X || a.Y != b.Y;
        }
        public static float Dot(Vector2 v1, Vector2 v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y;
        }
        public float LengthSquared()
        {
            return X * X + Y * Y;
        }
        public static Vector2 operator*(Vector2 v1, float v2)
        {
            return new Vector2(
                v1.X * v2,
                v1.Y * v2
                );
        }
        public static Vector2 operator-(Vector2 v1, Vector2 v2)
        {
            return new Vector2(
                v1.X - v2.X,
                v1.Y - v2.Y
                );
        }
    }
    public struct Vector3
    {
        public float X;
        public float Y;
        public float Z;
        public Vector3(float x, float y, float z)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
        }

        public static float Dot(Vector3 v1, Vector3 v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
        }
        public float LengthSquared()
        {
            return X * X + Y * Y + Z * Z;
        }
        public static Vector3 operator*(Vector3 v1, float v2)
        {
            return new Vector3(
                v1.X * v2,
                v1.Y * v2,
                v1.Z * v2);
        }
        public static Vector3 operator-(Vector3 v1, Vector3 v2)
        {
            return new Vector3(
                v1.X - v2.X,
                v1.Y - v2.Y,
                v1.Z - v2.Z);
        }
    }
#endif
}
