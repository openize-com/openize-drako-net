using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Compression
{
    /// <summary>
    /// Struct used for storing data about a source face that connects to an
    /// already traversed face that was either the initial face or a face encoded
    /// with either topology S (split) symbol. Such connection can be only caused by
    /// topology changes on the traversed surface (if its genus != 0, i.e. when the
    /// surface has topological handles or holes).
    /// For each occurence of such event we always encode the split symbol id, source
    /// symbol id and source edge id (left, or right). There will be always exectly
    /// two occurences of this event for every topological handle on the traversed
    /// mesh and one occurence for a hole.
    /// </summary>
    class TopologySplitEventData
    {
        internal int splitSymbolId;
        internal int sourceSymbolId;
        // We need to use uint instead of EdgeFaceName because the most recent
        // version of gcc does not allow that when optimizations are turned on.
        internal byte sourceEdge;
        internal byte splitEdge;
    }

    /// <summary>
    /// Hole event is used to store info about the first symbol that reached a
    /// vertex of so far unvisited hole. This can happen only on either the initial
    /// face or during a regular traversal when TOPOLOGYS is encountered.
    /// </summary>
    struct HoleEventData
    {
        internal int symbolId;

        public HoleEventData(int symId)
        {
            this.symbolId = symId;
        }
    }
}
