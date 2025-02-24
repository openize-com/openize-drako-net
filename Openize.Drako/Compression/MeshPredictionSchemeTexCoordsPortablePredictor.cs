using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Compression
{
    class MeshPredictionSchemeTexCoordsPortablePredictor
    {
        private const int kNumComponents = 2;
        public PointAttribute pos_attribute_;
        public int[] entry_to_point_id_map_;
        public int[] predicted_value_ = new int[2];
        // Encoded / decoded array of UV flips.
        // TODO(ostava): We should remove this and replace this with in-place encoding
        // and decoding to avoid unnecessary copy.
        private bool[] orientations_ = new bool[10];
        private int orientationCount = 0;
        MeshPredictionSchemeData mesh_data_;


        public MeshPredictionSchemeTexCoordsPortablePredictor(MeshPredictionSchemeData md)
        {
            pos_attribute_ = null;
            entry_to_point_id_map_ = null;
            mesh_data_ = md;
        }

        bool IsInitialized() { return pos_attribute_ != null; }

        LongVector3 GetPositionForEntryId(int entry_id)
        {
            int point_id = entry_to_point_id_map_[entry_id];
            var pos = pos_attribute_.ConvertValue(pos_attribute_.MappedIndex(point_id));
            return pos;
        }

        LongVector3 GetTexCoordForEntryId(int entry_id, Span<int> data)
        {
            int data_offset = entry_id * 2;
            return new LongVector3 (data[data_offset], data[data_offset + 1] );
        }

        // Computes predicted UV coordinates on a given corner. The coordinates are
        // stored in |predicted_value_| member.
        public bool ComputePredictedValue(bool is_encoder_t, int corner_id, Span<int> data, int data_id)
        {
            // Compute the predicted UV coordinate from the positions on all corners
            // of the processed triangle. For the best prediction, the UV coordinates
            // on the next/previous corners need to be already encoded/decoded.
            int next_corner_id = mesh_data_.CornerTable.Next(corner_id);
            int prev_corner_id =
                mesh_data_.CornerTable.Previous(corner_id);
            // Get the encoded data ids from the next and previous corners.
            // The data id is the encoding order of the UV coordinates.
            int next_data_id, prev_data_id;

            int next_vert_id, prev_vert_id;
            next_vert_id = mesh_data_.CornerTable.Vertex(next_corner_id);
            prev_vert_id = mesh_data_.CornerTable.Vertex(prev_corner_id);

            next_data_id = mesh_data_.vertexToDataMap[next_vert_id];
            prev_data_id = mesh_data_.vertexToDataMap[prev_vert_id];

            if (prev_data_id < data_id && next_data_id < data_id)
            {
                // Both other corners have available UV coordinates for prediction.
                var n_uv = GetTexCoordForEntryId(next_data_id, data);
                var p_uv = GetTexCoordForEntryId(prev_data_id, data);
                if (DracoUtils.VecEquals(p_uv, n_uv))
                {
                    // We cannot do a reliable prediction on degenerated UV triangles.
                    predicted_value_[0] = (int) p_uv.x;
                    predicted_value_[1] = (int) p_uv.y;
                    return true;
                }

                // Get positions at all corners.
                var tip_pos = GetPositionForEntryId(data_id);
                var next_pos = GetPositionForEntryId(next_data_id);
                var prev_pos = GetPositionForEntryId(prev_data_id);
                // We use the positions of the above triangle to predict the texture
                // coordinate on the tip corner C.
                // To convert the triangle into the UV coordinate system we first compute
                // position X on the vector |prev_pos - next_pos| that is the projection of
                // point C onto vector |prev_pos - next_pos|:
                //
                //              C
                //             /.  \
                //            / .     \
                //           /  .        \
                //          N---X----------P
                //
                // Where next_pos is point (N), prev_pos is point (P) and tip_pos is the
                // position of predicted coordinate (C).
                //
                var pn = DracoUtils.Sub(prev_pos, next_pos);
                uint pn_norm2_squared = DracoUtils.SquaredNorm(pn);
                if (pn_norm2_squared != 0)
                {
                    // Compute the projection of C onto PN by computing dot product of CN with
                    // PN and normalizing it by length of PN. This gives us a factor |s| where
                    // |s = PN.Dot(CN) / PN.SquaredNorm2()|. This factor can be used to
                    // compute X in UV space |X_UV| as |X_UV = N_UV + s * PN_UV|.
                    var cn = DracoUtils.Sub(tip_pos, next_pos);
                    long cn_dot_pn = DracoUtils.Dot(pn, cn);

                    var pn_uv = DracoUtils.Sub(p_uv, n_uv);
                    // Because we perform all computations with integers, we don't explicitly
                    // compute the normalized factor |s|, but rather we perform all operations
                    // over UV vectors in a non-normalized coordinate system scaled with a
                    // scaling factor |pn_norm2_squared|:
                    //
                    //      x_uv = X_UV * PN.Norm2Squared()
                    //
                    var x_uv = DracoUtils.Add(DracoUtils.Mul(n_uv, pn_norm2_squared), DracoUtils.Mul(pn_uv, cn_dot_pn));

                    // Compute squared length of vector CX in position coordinate system:
                    var x_pos =
                        DracoUtils.Add(next_pos, DracoUtils.Div(DracoUtils.Mul(pn, cn_dot_pn), pn_norm2_squared));
                    long cx_norm2_squared = DracoUtils.SquaredNorm(DracoUtils.Sub(tip_pos, x_pos));

                    // Compute vector CX_UV in the uv space by rotating vector PN_UV by 90
                    // degrees and scaling it with factor CX.Norm2() / PN.Norm2():
                    //
                    //     CX_UV = (CX.Norm2() / PN.Norm2()) * Rot(PN_UV)
                    //
                    // To preserve precision, we perform all operations in scaled space as
                    // explained above, so we want the final vector to be:
                    //
                    //     cx_uv = CX_UV * PN.Norm2Squared()
                    //
                    // We can then rewrite the formula as:
                    //
                    //     cx_uv = CX.Norm2() * PN.Norm2() * Rot(PN_UV)
                    //
                    var cx_uv = new LongVector3(pn_uv.y, -pn_uv.x); // Rotated PN_UV.
                    // Compute CX.Norm2() * PN.Norm2()
                    uint norm_squared =
                        (uint)DracoUtils.IntSqrt((ulong)(cx_norm2_squared * pn_norm2_squared));
                    // Final cx_uv in the scaled coordinate space.
                    cx_uv = DracoUtils.Mul(cx_uv, norm_squared);

                    // Predicted uv coordinate is then computed by either adding or
                    // subtracting CX_UV to/from X_UV.
                    LongVector3 predicted_uv;
                    if (is_encoder_t)
                    {
                        // When encoding, compute both possible vectors and determine which one
                        // results in a better prediction.
                        // Both vectors need to be transformed back from the scaled space to
                        // the real UV coordinate space.
                        var predicted_uv_0 = DracoUtils.Div(DracoUtils.Add(x_uv, cx_uv), pn_norm2_squared);
                        var predicted_uv_1 = DracoUtils.Div(DracoUtils.Sub(x_uv, cx_uv), pn_norm2_squared);
                        var c_uv = GetTexCoordForEntryId(data_id, data);
                        if(orientationCount == orientations_.Length)
                        {
                            Array.Resize(ref orientations_, orientations_.Length + (orientations_.Length >> 1));
                        }
                        if (DracoUtils.SquaredNorm(DracoUtils.Sub(c_uv, predicted_uv_0)) <
                            DracoUtils.SquaredNorm(DracoUtils.Sub(c_uv, predicted_uv_1)))
                        {
                            predicted_uv = predicted_uv_0;
                            orientations_[orientationCount++] = true;
                        }
                        else
                        {
                            predicted_uv = predicted_uv_1;
                            orientations_[orientationCount++] = false;
                        }
                    }
                    else
                    {
                        // When decoding the data, we already know which orientation to use.
                        if (orientationCount == 0)
                            return false;
                        //remove last orientation
                        bool orientation = orientations_[orientationCount - 1];
                        orientationCount--;
                        if (orientation)
                            predicted_uv = DracoUtils.Div(DracoUtils.Add(x_uv, cx_uv), pn_norm2_squared);
                        else
                            predicted_uv = DracoUtils.Div(DracoUtils.Sub(x_uv, cx_uv), pn_norm2_squared);
                    }

                    predicted_value_[0] = (int) (predicted_uv.x);
                    predicted_value_[1] = (int) (predicted_uv.y);
                    return true;
                }
            }

            // Else we don't have available textures on both corners or the position data
            // is invalid. For such cases we can't use positions for predicting the uv
            // value and we resort to delta coding.
            int data_offset = 0;
            if (prev_data_id < data_id)
            {
                // Use the value on the previous corner as the prediction.
                data_offset = prev_data_id * kNumComponents;
            }

            if (next_data_id < data_id)
            {
                // Use the value on the next corner as the prediction.
                data_offset = next_data_id * kNumComponents;
            }
            else
            {
                // None of the other corners have a valid value. Use the last encoded value
                // as the prediction if possible.
                if (data_id > 0)
                {
                    data_offset = (data_id - 1) * kNumComponents;
                }
                else
                {
                    // We are encoding the first value. Predict 0.
                    for (int i = 0; i < kNumComponents; ++i)
                    {
                        predicted_value_[i] = 0;
                    }

                    return true;
                }
            }

            for (int i = 0; i < kNumComponents; ++i)
            {
                predicted_value_[i] = data[data_offset + i];
            }

            return true;
        }

        public bool orientation(int i) { return orientations_[i]; }
        public void set_orientation(int i, bool v) { orientations_[i] = v; }
        public int num_orientations() { return orientationCount; }
        public void ResizeOrientations(int num_orientations)
        {
            Array.Resize(ref orientations_, num_orientations);
            orientationCount = num_orientations;
        }

    }

}
