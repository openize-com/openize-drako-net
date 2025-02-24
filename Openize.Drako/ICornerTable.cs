using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako
{
    abstract class ICornerTable
    {
        public abstract int NumFaces { get; }
        public abstract int NumVertices { get; }

        public static int LocalIndex(int corner)
        {
            return corner % 3;
        }

        public int Next(int corner)
        {
            if (corner < 0)
                return corner;
            return LocalIndex(++corner) != 0 ? corner : corner - 3;
        }

        public int Previous(int corner)
        {
            if (corner < 0)
                return corner;
            return LocalIndex(corner) != 0 ? corner - 1 : corner + 2;
        }



        public abstract int Vertex(int corner);
        public abstract bool IsOnBoundary(int vert);

        public abstract int Opposite(int corner);
        public abstract int GetRightCorner(int cornerId);
        public abstract int GetLeftCorner(int cornerId);

        public abstract int LeftMostCorner(int v);
        public abstract int SwingRight(int corner);
        public abstract int SwingLeft(int corner);

        //int ConfidentVertex(int corner);
        //int Valence(int v);
    }
}
