using FileFormat.Drako.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Compression
{
    enum NormalPredictionMode : byte
    {
        OneTriangle,
        TriangleArea

    }

    abstract class MeshPredictionSchemeGeometricNormalPredictorBase
    {
        public PointAttribute pos_attribute_;
        public int[] entry_to_point_id_map_;
        protected MeshPredictionSchemeData mesh_data_;
        protected NormalPredictionMode normal_prediction_mode_;

        public MeshPredictionSchemeGeometricNormalPredictorBase(MeshPredictionSchemeData meshData)
        {
            this.mesh_data_ = meshData;

        }

        public bool IsInitialized
        {
            get
            {
                if (pos_attribute_ == null)
                    return false;
                if (entry_to_point_id_map_ == null)
                    return false;
                return true;
            }
        }

        public abstract bool SetNormalPredictionMode(NormalPredictionMode mode);

        public virtual NormalPredictionMode GetNormalPredictionMode()
        {
            return normal_prediction_mode_;
        }

        protected LongVector3 GetPositionForDataId(int data_id)
        {
            var point_id = entry_to_point_id_map_[data_id];
            var pos_val_id = pos_attribute_.MappedIndex(point_id);
            //long[] pos = new long[3];
            return pos_attribute_.ConvertValue(pos_val_id);
        }

        protected LongVector3 GetPositionForCorner(int ci)
        {
            var corner_table = mesh_data_.CornerTable;
            var vert_id = corner_table.Vertex(ci);
            var data_id = mesh_data_.vertexToDataMap[vert_id];
            return GetPositionForDataId(data_id);
        }

        protected int[] GetOctahedralCoordForDataId(int data_id, int[] data)
        {
            int data_offset = data_id * 2;
            return new int[] {data[data_offset], data[data_offset + 1]};
        }

        // Computes predicted octahedral coordinates on a given corner.
        public abstract void ComputePredictedValue(int corner_id, int[] prediction);

    }

    class MeshPredictionSchemeGeometricNormalPredictorArea : MeshPredictionSchemeGeometricNormalPredictorBase
    {
        public MeshPredictionSchemeGeometricNormalPredictorArea(MeshPredictionSchemeData meshData)
            : base(meshData)
        {
            SetNormalPredictionMode(NormalPredictionMode.TriangleArea);
        }

        // Computes predicted octahedral coordinates on a given corner.
        public override void ComputePredictedValue(int corner_id, int[] prediction)
        {
            //typedef typename MeshDataT::CornerTable CornerTable;
            ICornerTable corner_table = mesh_data_.CornerTable;
            // Going to compute the predicted normal from the surrounding triangles
            // according to the connectivity of the given corner table.
            VertexCornersIterator cit = VertexCornersIterator.FromCorner(corner_table, corner_id);
            // Position of central vertex does not change in loop.
            var pos_cent = this.GetPositionForCorner(corner_id);
            // Computing normals for triangles and adding them up.

            var normal = new LongVector3();
            int c_next, c_prev;
            while (!cit.End)
            {
                // Getting corners.
                if (this.normal_prediction_mode_ == NormalPredictionMode.OneTriangle)
                {
                    c_next = corner_table.Next(corner_id);
                    c_prev = corner_table.Previous(corner_id);
                }
                else
                {
                    c_next = corner_table.Next(cit.Corner);
                    c_prev = corner_table.Previous(cit.Corner);
                }

                var pos_next = this.GetPositionForCorner(c_next);
                var pos_prev = this.GetPositionForCorner(c_prev);

                // Computing delta vectors to next and prev.
                var delta_next= DracoUtils.Sub(pos_next, pos_cent);
                var delta_prev = DracoUtils.Sub(pos_prev, pos_cent);

                // Computing cross product.
                var cross = DracoUtils.CrossProduct(delta_next, delta_prev);
                normal = DracoUtils.Add(normal, cross);
                cit.Next();
            }

            // Convert to int32_t, make sure entries are not too large.
            long upper_bound = 1 << 29;
            if (this.normal_prediction_mode_ == NormalPredictionMode.OneTriangle)
            {
                int abs_sum = (int)DracoUtils.AbsSum(normal);
                if (abs_sum > upper_bound)
                {
                    long quotient = abs_sum / upper_bound;
                    normal = DracoUtils.Div(normal, quotient);
                }
            }
            else
            {
                long abs_sum = DracoUtils.AbsSum(normal);
                if (abs_sum > upper_bound)
                {
                    long quotient = abs_sum / upper_bound;
                    normal = DracoUtils.Div(normal, quotient);
                }
            }

            prediction[0] = (int) normal.x;
            prediction[1] = (int) normal.y;
            prediction[2] = (int) normal.z;
        }

        public override bool SetNormalPredictionMode(NormalPredictionMode mode)
        {
            if (mode == NormalPredictionMode.OneTriangle)
            {
                this.normal_prediction_mode_ = mode;
                return true;
            }
            else if (mode == NormalPredictionMode.TriangleArea)
            {
                this.normal_prediction_mode_ = mode;
                return true;
            }

            return false;
        }
    }
}
