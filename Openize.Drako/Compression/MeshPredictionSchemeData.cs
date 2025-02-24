using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Compression
{
    class MeshPredictionSchemeData
    {
        private DracoMesh mesh;
        private ICornerTable cornerTable;

        /// <summary>
        /// Mapping between vertices and their encoding order. I.e. when an attribute
        /// entry on a given vertex was encoded.
        /// </summary>
        internal int[] vertexToDataMap;

        /// <summary>
        /// Array that stores which corner was processed when a given attribute entry
        /// was encoded or decoded.
        /// </summary>
        internal IntList dataToCornerMap;

        public MeshPredictionSchemeData(DracoMesh mesh, ICornerTable table, IntList dataToCornerMap,
            int[] vertexToDataMap)
        {
            this.mesh = mesh;
            cornerTable = table;
            this.dataToCornerMap = dataToCornerMap;
            this.vertexToDataMap = vertexToDataMap;
        }

        public ICornerTable CornerTable
        {
            get { return cornerTable;}
        }
    }
}
