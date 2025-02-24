using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako
{
    class MeshAttributeCornerTable : ICornerTable
    {
        private const int kInvalidVertexIndex = -1;
        private const int kInvalidCornerIndex = -1;
        private bool[] isEdgeOnSeam;
        private bool[] isVertexOnSeam;

        /// <summary>
        /// If this is set to true, it means that there are no attribute seams between
        /// two faces. This can be used to speed up some algorithms.
        /// </summary>
        private bool noInteriorSeams = true;

        private int[] cornerToVertexMap;

        /// <summary>
        /// Map between vertices and their associated left most corners. A left most
        /// corner is a corner that is adjecent to a boundary or an attribute seam from
        /// right (i.e., SwingLeft from that corner will return an invalid corner). If
        /// no such corner exists for a given vertex, then any corner attached to the
        /// vertex can be used.
        /// </summary>
        private IntList vertexToLeftMostCornerMap = new IntList();

        /// <summary>
        /// Map between vertex ids and attribute entry ids (i.e. the values stored in
        /// the attribute buffer). The attribute entry id can be retrieved using the
        /// VertexParent() method.
        /// </summary>
        private IntList vertexToAttributeEntryIdMap = new IntList();
        private CornerTable cornerTable;

        public MeshAttributeCornerTable(CornerTable table)
        {

            isEdgeOnSeam = new bool[table.NumCorners];
            isVertexOnSeam = new bool[table.NumVertices];
            cornerToVertexMap = new int[table.NumCorners];
            for (int i = 0; i < cornerToVertexMap.Length; i++)
                cornerToVertexMap[i] = kInvalidCornerIndex;
            vertexToAttributeEntryIdMap.Capacity = table.NumVertices;
            vertexToLeftMostCornerMap.Capacity = table.NumVertices;
            cornerTable = table;
            noInteriorSeams = true;
        }

        public MeshAttributeCornerTable(DracoMesh mesh, CornerTable table, PointAttribute att)
            : this(table)
        {

            // Find all necessary data for encoding attributes. For now we check which of
            // the mesh vertices is part of an attribute seam, because seams require
            // special handling.
            for (int c = 0; c < cornerTable.NumCorners; ++c)
            {
                int f = cornerTable.Face(c);
                if (cornerTable.IsDegenerated(f))
                    continue; // Ignore corners on degenerated faces.
                int oppCorner = cornerTable.Opposite(c);
                if (oppCorner < 0)
                {
                    // Boundary. Mark it as seam edge.
                    isEdgeOnSeam[c] = true;
                    // Mark seam vertices.
                    int v = cornerTable.Vertex(cornerTable.Next(c));
                    isVertexOnSeam[v] = true;
                    v = cornerTable.Vertex(cornerTable.Previous(c));
                    isVertexOnSeam[v] = true;
                    continue;
                }
                if (oppCorner < c)
                    continue; // Opposite corner was already processed.

                int actC = c, actSiblingC = oppCorner;
                for (int i = 0; i < 2; ++i)
                {
                    // Get the sibling corners. I.e., the two corners attached to the same
                    // vertex but divided by the seam edge.
                    actC = cornerTable.Next(actC);
                    actSiblingC = cornerTable.Previous(actSiblingC);
                    int pointId = DracoUtils.CornerToPointId(actC, mesh);
                    int siblingPointId = DracoUtils.CornerToPointId(actSiblingC, mesh);
                    if (att.MappedIndex(pointId) != att.MappedIndex(siblingPointId))
                    {
                        noInteriorSeams = false;
                        isEdgeOnSeam[c] = true;
                        isEdgeOnSeam[oppCorner] = true;
                        // Mark seam vertices.
                        isVertexOnSeam[cornerTable.Vertex(cornerTable.Next(c))] = true;
                        isVertexOnSeam[cornerTable.Vertex(cornerTable.Previous(c))] = true;
                        isVertexOnSeam[cornerTable.Vertex(cornerTable.Next(oppCorner))] = true;
                        isVertexOnSeam[cornerTable.Vertex(cornerTable.Previous(oppCorner))] = true;
                        break;
                    }
                }
            }
            RecomputeVertices(mesh, att);
        }

        public void AddSeamEdge(int c)
        {
            isEdgeOnSeam[c] = true;
            // Mark seam vertices.
            isVertexOnSeam[cornerTable.Vertex(cornerTable.Next(c))] = true;
            isVertexOnSeam[cornerTable.Vertex(cornerTable.Previous(c))] = true;
            int oppCorner = cornerTable.Opposite(c);
            if (oppCorner >= 0)
            {
                noInteriorSeams = false;
                isEdgeOnSeam[oppCorner] = true;
                isVertexOnSeam[cornerTable.Vertex(cornerTable.Next(oppCorner))] = true;
                isVertexOnSeam[cornerTable.Vertex(cornerTable.Previous(oppCorner))] = true;
            }
        }

        /// <summary>
        /// Recomputes vertices using the newly added seam edges (needs to be called
        /// whenever the seam edges are updated).
        /// |mesh| and |att| can be null, in which case mapping between vertices and
        /// attribute value ids is set to identity.
        /// </summary>
        public void RecomputeVertices(DracoMesh mesh, PointAttribute att)
        {

            if (mesh != null && att != null)
            {
                RecomputeVerticesInternal(true, mesh, att);
            }
            else
            {
                RecomputeVerticesInternal(false, null, null);
            }
        }

        public void RecomputeVerticesInternal(bool initVertexToAttributeEntryMap, DracoMesh mesh, PointAttribute att)
        {
            int numNewVertices = 0;
            for (int v = 0; v < cornerTable.NumVertices; ++v)
            {
                int c = cornerTable.LeftMostCorner(v);
                if (c < 0)
                    continue; // Isolated vertex?
                int firstVertId = numNewVertices++;
                if (initVertexToAttributeEntryMap)
                {
                    int pointId = DracoUtils.CornerToPointId(c, mesh);
                    vertexToAttributeEntryIdMap.Add(att.MappedIndex(pointId));
                }
                else
                {
                    // Identity mapping
                    vertexToAttributeEntryIdMap.Add(firstVertId);
                }
                int firstC = c;
                int actC;
                // Check if the vertex is on a seam edge, if it is we need to find the first
                // attribute entry on the seam edge when traversing in the ccw direction.
                if (isVertexOnSeam[v])
                {
                    // Try to swing left on the modified corner table. We need to get the
                    // first corner that defines an attribute seam.
                    actC = SwingLeft(firstC);
                    while (actC >= 0)
                    {
                        firstC = actC;
                        actC = SwingLeft(actC);
                    }
                }
                cornerToVertexMap[firstC] = firstVertId;
                vertexToLeftMostCornerMap.Add(firstC);
                actC = cornerTable.SwingRight(firstC);
                while (actC >= 0 && actC != firstC)
                {
                    if (IsCornerOppositeToSeamEdge(cornerTable.Next(actC)))
                    {
                        firstVertId = numNewVertices++;
                        if (initVertexToAttributeEntryMap)
                        {
                            int pointId = DracoUtils.CornerToPointId(actC, mesh);
                            vertexToAttributeEntryIdMap.Add(att.MappedIndex(pointId));
                        }
                        else
                        {
                            // Identity mapping.
                            vertexToAttributeEntryIdMap.Add(firstVertId);
                        }
                        vertexToLeftMostCornerMap.Add(actC);
                    }
                    cornerToVertexMap[actC] = firstVertId;
                    actC = cornerTable.SwingRight(actC);
                }
            }
        }

        public bool IsCornerOppositeToSeamEdge(int corner)
        {
            return isEdgeOnSeam[corner];
        }

        public override int Opposite(int corner)
        {
            if (IsCornerOppositeToSeamEdge(corner))
                return kInvalidCornerIndex;
            return cornerTable.Opposite(corner);
        }

        /// <summary>
        /// Returns true when a corner is attached to any attribute seam.
        /// </summary>
        public bool IsCornerOnSeam(int corner)
        {
            return isVertexOnSeam[cornerTable.Vertex(corner)];
        }

        /// <summary>
        /// Similar to CornerTable::GetLeftCorner and CornerTable::GetRightCorner, but
        /// does not go over seam edges.
        /// </summary>
        public override int GetLeftCorner(int corner)
        {
            return Opposite(Previous(corner));
        }

        public override int GetRightCorner(int corner)
        {
            return Opposite(Next(corner));
        }

        /// <summary>
        /// Similar to CornerTable::SwingRight, but it does not go over seam edges.
        /// </summary>
        /// <param name="corner"></param>
        /// <returns></returns>
        public override int SwingRight(int corner)
        {
            return Previous(Opposite(Previous(corner)));
        }

        /// <summary>
        /// Similar to CornerTable.SwingLeft, but it does not go over seam edges.
        /// </summary>
        /// <param name="corner"></param>
        /// <returns></returns>
        public override int SwingLeft(int corner)
        {
            return Next(Opposite(Next(corner)));
        }

        public override int NumVertices
        {
            get { return vertexToAttributeEntryIdMap.Count; }
        }

        public override int NumFaces
        {
            get { return cornerTable.NumFaces; }
        }

        public override int Vertex(int corner)
        {
            return cornerToVertexMap[corner];
        }

        // Returns the attribute entry id associated to the given vertex.
        public int VertexParent(int vert)
        {
            return vertexToAttributeEntryIdMap[vert];
        }

        public override int LeftMostCorner(int v)
        {
            return vertexToLeftMostCornerMap[v];
        }

        public override bool IsOnBoundary(int vert)
        {
            int corner = LeftMostCorner(vert);
            if (corner < 0)
                return true;
            return IsCornerOnSeam(corner);
        }

        public bool NoInteriorSeams
        {
            get { return noInteriorSeams; }
        }

        CornerTable CornerTable
        {
            get { return cornerTable; }
        }

    }


}
