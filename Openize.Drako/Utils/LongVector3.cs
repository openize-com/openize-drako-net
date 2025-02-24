using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Openize.Drako.Utils
{
    struct IntVector
    {
        public int x;
        public int y;
        public IntVector(int x, int y)
        {
            this.x = x;
            this.y = y;
        }
            public static IntVector operator -(IntVector a, IntVector b)
            {
                return new IntVector(a.x - b.x, a.y - b.y);
            }

            public static IntVector operator +(IntVector a, IntVector b)
            {
                return new IntVector(a.x + b.x, a.y + b.y);
            }

    }
    struct LongVector3
    {
        public long x;
        public long y;
        public long z;
        public LongVector3(long x, long y)
        {
            this.x = x;
            this.y = y;
            this.z = 0;
        }
        public LongVector3(long x, long y, long z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

    }
}
