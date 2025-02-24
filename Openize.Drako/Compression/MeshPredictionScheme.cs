using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Compression
{

    /// <summary>
    /// Base class for all mesh prediction schemes that use the mesh connectivity
    /// data. |MeshDataT| can be any class that provides the same interface as the
    /// PredictionSchemeMeshData class.
    /// </summary>
    abstract class MeshPredictionScheme : PredictionScheme
    {
        protected MeshPredictionSchemeData meshData;
        protected MeshPredictionScheme(PointAttribute attribute,
            PredictionSchemeTransform transform, MeshPredictionSchemeData meshData)
            : base(attribute, transform)
        {
            this.meshData = meshData;
        }

        public override bool Initialized { get { return true; } }
        protected static void GetParallelogramEntries(
            int ci, ICornerTable table,
            int[] vertexToDataMap, ref int oppEntry,
            ref int nextEntry, ref int prevEntry)
        {
            // One vertex of the input |table| correspond to exactly one attribute value
            // entry. The |table| can be either CornerTable for per-vertex attributes,
            // or MeshAttributeCornerTable for attributes with interior seams.
            oppEntry = vertexToDataMap[table.Vertex(ci)];
            nextEntry = vertexToDataMap[table.Vertex(table.Next(ci))];
            prevEntry = vertexToDataMap[table.Vertex(table.Previous(ci))];
        }
    }
}
