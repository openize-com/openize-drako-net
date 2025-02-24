using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Compression
{

    // ValenceCache provides support for the caching of valences off of some kind of
    // CornerTable 'type' of class.
    // No valences should be queried before Caching is
    // performed and values should be removed/recached when changes to the
    // underlying mesh are taking place.
    class ValenceCache
    {
        private CornerTable table_;
        private sbyte[] vertex_valence_cache_8_bit_;
        private int[] vertex_valence_cache_32_bit_;
        public ValenceCache(CornerTable table)
        {
            this.table_ = table;
        }

        // Do not call before CacheValences() / CacheValencesInaccurate().
        public sbyte ValenceFromCacheInaccurateC(int c)
        {
            if (c == -1)
                return -1;
            return ValenceFromCacheInaccurateV(table_.Vertex(c));
        }
        public int ValenceFromCacheC(int c)
        {
            if (c == -1)
                return -1;
            return ValenceFromCacheV(table_.Vertex(c));
        }
        public int ConfidentValenceFromCacheV(int v)
        {
            //DRACO_DCHECK_LT(v, table_.NumVertices);
            //DRACO_DCHECK_EQ(vertex_valence_cache_32_bit_.Length, table_.NumVertices);
            return vertex_valence_cache_32_bit_[v];
        }

        // Collect the valence for all vertices so they can be reused later.  The
        // 'inaccurate' versions of this family of functions clips the true valence
        // of the vertices to 8 signed bits as a space optimization.  This clipping
        // will lead to occasionally wrong results.  If accurate results are required
        // under all circumstances, do not use the 'inaccurate' version or else
        // use it and fetch the correct result in the event the value appears clipped.
        // The topology of the mesh should be a constant when Valence Cache functions
        // are being used.  Modification of the mesh while cache(s) are filled will
        // not guarantee proper results on subsequent calls unless they are rebuilt.
        public void CacheValencesInaccurate()
        {
            if (vertex_valence_cache_8_bit_ == null)
            {
                int vertex_count = table_.NumVertices;
                vertex_valence_cache_8_bit_ = new sbyte[vertex_count];
                for (int v = 0; v < vertex_count; v += 1)
                    vertex_valence_cache_8_bit_[v] = (sbyte)(
                        Math.Min(sbyte.MaxValue, table_.Valence(v)));
            }
        }
        public void CacheValences()
        {
            if (vertex_valence_cache_32_bit_ == null)
            {
                int vertex_count = table_.NumVertices;
                vertex_valence_cache_32_bit_ = new int[vertex_count];
                for (int v = 0; v < vertex_count; v += 1)
                    vertex_valence_cache_32_bit_[v] = table_.Valence(v);
            }
        }

        public sbyte ConfidentValenceFromCacheInaccurateC(int c)
        {
            //DRACO_DCHECK_GE(c, 0);
            return ConfidentValenceFromCacheInaccurateV(table_.ConfidentVertex(c));
        }
        public int ConfidentValenceFromCacheC(int c)
        {
            //DRACO_DCHECK_GE(c, 0);
            return ConfidentValenceFromCacheV(table_.ConfidentVertex(c));
        }
        public sbyte ValenceFromCacheInaccurateV(int v)
        {
            //DRACO_DCHECK_EQ(vertex_valence_cache_8_bit_.Length, table_.NumVertices);
            if (v == -1 || v >= table_.NumVertices)
                return -1;
            return ConfidentValenceFromCacheInaccurateV(v);
        }
        public sbyte ConfidentValenceFromCacheInaccurateV(int v)
        {
            return vertex_valence_cache_8_bit_[v];
        }

        // TODO(draco-eng) Add unit tests for ValenceCache functions.
        public int ValenceFromCacheV(int v)
        {
            //DRACO_DCHECK_EQ(vertex_valence_cache_32_bit_.Length, table_.NumVertices);
            if (v == -1 || v >= table_.NumVertices)
                return -1;
            return ConfidentValenceFromCacheC(v);
        }

        // Clear the cache of valences and deallocate the memory.
        public void ClearValenceCacheInaccurate()
        {
            vertex_valence_cache_8_bit_ = null;
        }
        public void ClearValenceCache()
        {
            vertex_valence_cache_32_bit_ = null;
        }

        public bool IsCacheEmpty()
        {
            return vertex_valence_cache_8_bit_ == null &&
                 vertex_valence_cache_32_bit_ == null;
        }


    }
}
