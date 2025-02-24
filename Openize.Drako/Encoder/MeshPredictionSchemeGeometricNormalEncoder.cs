using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Compression;
using FileFormat.Drako.Encoder;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Compression
{
    partial class MeshPredictionSchemeGeometricNormal
    {
        private RAnsBitEncoder flip_normal_bit_encoder_ = new RAnsBitEncoder();
        public override void ComputeCorrectionValues(Span<int> in_data, Span<int> out_corr, int size, int num_components,
            int[] entry_to_point_id_map)
        {

            octahedron_tool_box_.SetQuantizationBits(((PredictionSchemeNormalOctahedronTransformBase) transform_).QuantizationBits);
            predictor_.entry_to_point_id_map_ = entry_to_point_id_map;

            flip_normal_bit_encoder_.StartEncoding();

            int corner_map_size = this.meshData.dataToCornerMap.Count;

            var pred_normal_3d = new int[3];
            Span<int> pos_pred_normal_oct = stackalloc int[2];
            Span<int> neg_pred_normal_oct = stackalloc int[2];
            Span<int> pos_correction = stackalloc int[2];
            Span<int> neg_correction = stackalloc int[2];
            for (int data_id = 0; data_id < corner_map_size; ++data_id)
            {
                int corner_id =
                    this.meshData.dataToCornerMap[data_id];
                predictor_.ComputePredictedValue(corner_id, pred_normal_3d);

                // Compute predicted octahedral coordinates.
                octahedron_tool_box_.CanonicalizeIntegerVector(pred_normal_3d);

                // Compute octahedral coordinates for both possible directions.
                int s, t;
                octahedron_tool_box_.IntegerVectorToQuantizedOctahedralCoords(pred_normal_3d, out s, out t);
                pos_pred_normal_oct[0] = s;
                pos_pred_normal_oct[1] = t;
                pred_normal_3d[0] = -pred_normal_3d[0];
                pred_normal_3d[1] = -pred_normal_3d[1];
                pred_normal_3d[2] = -pred_normal_3d[2];

                octahedron_tool_box_.IntegerVectorToQuantizedOctahedralCoords(pred_normal_3d, out s, out t);
                neg_pred_normal_oct[0] = s;
                neg_pred_normal_oct[1] = t;

                // Choose the one with the best correction value.
                 int data_offset = data_id * 2;

                this.transform_.ComputeCorrection(in_data.Slice(data_offset),
                    pos_pred_normal_oct,
                    pos_correction, 0);
                this.transform_.ComputeCorrection(in_data.Slice(data_offset),
                    neg_pred_normal_oct,
                    neg_correction, 0);
                pos_correction[0] = octahedron_tool_box_.ModMax(pos_correction[0]);
                pos_correction[1] = octahedron_tool_box_.ModMax(pos_correction[1]);
                neg_correction[0] = octahedron_tool_box_.ModMax(neg_correction[0]);
                neg_correction[1] = octahedron_tool_box_.ModMax(neg_correction[1]);
                if (DracoUtils.AbsSum(pos_correction) < DracoUtils.AbsSum(neg_correction))
                {
                    flip_normal_bit_encoder_.EncodeBit(false);
                    out_corr[data_offset] = octahedron_tool_box_.MakePositive(pos_correction[0]);
                    out_corr[data_offset + 1] = octahedron_tool_box_.MakePositive(pos_correction[1]);
                }
                else
                {
                    flip_normal_bit_encoder_.EncodeBit(true);
                    out_corr[data_offset] = octahedron_tool_box_.MakePositive(neg_correction[0]);
                    out_corr[data_offset + 1] = octahedron_tool_box_.MakePositive(neg_correction[1]);
                }
            }

        }

        public override void EncodePredictionData(EncoderBuffer buffer)
        {

            this.transform_.EncodeTransformData(buffer);

            // Encode normal flips.
            flip_normal_bit_encoder_.EndEncoding(buffer);
        }
    }

}
