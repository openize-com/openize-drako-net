using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako
{
    struct VertexCornersIterator
    {
        // Create the iterator from the provided corner table and the central vertex.
        public static VertexCornersIterator FromVertex(ICornerTable table, int vert_id)
        {
            var ret = new VertexCornersIterator();
            ret.corner_table_ = table;
            ret.start_corner_ = table.LeftMostCorner(vert_id);
            ret.corner_ = ret.start_corner_;
            ret.left_traversal_ = true;
            return ret;

        } // Gets the last visited corner.
        public static VertexCornersIterator FromCorner(ICornerTable table, int corner_id)
        {
            var ret = new VertexCornersIterator();
            ret.corner_table_ = table;
            ret.start_corner_ = corner_id;
            ret.corner_ = ret.start_corner_;
            ret.left_traversal_ = true;
            return ret;
        } 

        public int Corner
        {
            get { return corner_; }
        }

        // Returns true when all ring vertices have been visited.
        public bool End
        {
            get { return corner_ == CornerTable.kInvalidCornerIndex; }
        }

        // Proceeds to the next corner if possible.
        public void Next()
        {
            if (left_traversal_)
            {
                corner_ = corner_table_.SwingLeft(corner_);
                if (corner_ == CornerTable.kInvalidCornerIndex)
                {
                    // Open boundary reached.
                    corner_ = corner_table_.SwingRight(start_corner_);
                    left_traversal_ = false;
                }
                else if (corner_ == start_corner_)
                {
                    // End reached.
                    corner_ = CornerTable.kInvalidCornerIndex;
                }
            }
            else
            {
                // Go to the right until we reach a boundary there (no explicit check
                // is needed in this case).
                corner_ = corner_table_.SwingRight(corner_);
            }
        }


        private ICornerTable corner_table_;

        // The first processed corner.
        private int start_corner_;

        // The last processed corner.
        private int corner_;

        // Traversal direction.
        private bool left_traversal_;
    }
}
