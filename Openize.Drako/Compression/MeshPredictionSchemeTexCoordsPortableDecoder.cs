using Openize.Drako.Decoder;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;
using Openize.Drako.Encoder;

namespace Openize.Drako.Compression
{
    class MeshPredictionSchemeTexCoordsPortableDecoder : MeshPredictionScheme
    {

        MeshPredictionSchemeTexCoordsPortablePredictor predictor_;
        public MeshPredictionSchemeTexCoordsPortableDecoder(PointAttribute attribute, PredictionSchemeTransform transform, MeshPredictionSchemeData meshData) : base(attribute, transform, meshData)
        {
            predictor_ = new MeshPredictionSchemeTexCoordsPortablePredictor(meshData);
        }

        public override void ComputeCorrectionValues(Span<int> in_data, Span<int> out_corr, int size, int num_components,
            int[] entry_to_point_id_map)
        {
            predictor_.entry_to_point_id_map_ = entry_to_point_id_map;
            this.transform_.InitializeEncoding(in_data, num_components);
            // We start processing from the end because this prediction uses data from
            // previous entries that could be overwritten when an entry is processed.
            for (int p = this.meshData.dataToCornerMap.Count - 1;
                p >= 0;
                --p)
            {
                int corner_id = this.meshData.dataToCornerMap[p];
                predictor_.ComputePredictedValue(true, corner_id, in_data, p);

                int dst_offset = p * num_components;
                this.transform_.ComputeCorrection(in_data, dst_offset,
                    predictor_.predicted_value_.AsSpan(), 0,
                    out_corr, dst_offset, 0);
            }
        }

        public override void EncodePredictionData(EncoderBuffer buffer)
        {

            // Encode the delta-coded orientations using arithmetic coding.
            int num_orientations = predictor_.num_orientations();
            buffer.Encode(num_orientations);
            bool last_orientation = true;
            RAnsBitEncoder encoder = new RAnsBitEncoder();
            encoder.StartEncoding();
            for (int i = 0; i < num_orientations; ++i)
            {
                bool orientation = predictor_.orientation(i);
                encoder.EncodeBit(orientation == last_orientation);
                last_orientation = orientation;
            }

            encoder.EndEncoding(buffer);
            base.EncodePredictionData(buffer);
        }


        public override void SetParentAttribute(PointAttribute att)
        {
            if (att == null || att.AttributeType != AttributeType.Position)
                throw DracoUtils.Failed();  // Invalid attribute type.
            if (att.ComponentsCount != 3)
                throw DracoUtils.Failed();  // Currently works only for 3 component positions.
            predictor_.pos_attribute_ = att;
        }

        public override void ComputeOriginalValues(Span<int> inCorr, Span<int> outData, int size, int numComponents, int[] entryToPointIdMap)
        {
            predictor_.entry_to_point_id_map_ = entryToPointIdMap;
            this.transform_.InitializeDecoding(numComponents);

            int corner_map_size = this.meshData.dataToCornerMap.Count;
            for (int p = 0; p < corner_map_size; ++p) {
                int corner_id = this.meshData.dataToCornerMap[p];
                if (!predictor_.ComputePredictedValue(false, corner_id, outData, p))
                    throw DracoUtils.Failed();

                int dst_offset = p * numComponents;
                this.transform_.ComputeOriginalValue(predictor_.predicted_value_.AsSpan(), 0,
                                                       inCorr, dst_offset,
                                                       outData, dst_offset);
            }
        }
        public override void DecodePredictionData(DecoderBuffer buffer) {
            // Decode the delta coded orientations.
            int num_orientations = buffer.DecodeI32();
            if (num_orientations < 0)
                throw DracoUtils.Failed();
            predictor_.ResizeOrientations(num_orientations);
            bool last_orientation = true;
            RAnsBitDecoder decoder = new RAnsBitDecoder();
            decoder.StartDecoding(buffer);
            for (int i = 0; i < num_orientations; ++i) {
                if (!decoder.DecodeNextBit())
                    last_orientation = !last_orientation;
                predictor_.set_orientation(i, last_orientation);
            }
            decoder.EndDecoding();
            base.DecodePredictionData(buffer);
        }
        public override AttributeType GetParentAttributeType(int i)
        {
            return AttributeType.Position;
        }
        public override int NumParentAttributes => 1;

        public override PredictionSchemeMethod PredictionMethod
        {
            get { return PredictionSchemeMethod.TexCoordsPortable; }
        }
    }
}
