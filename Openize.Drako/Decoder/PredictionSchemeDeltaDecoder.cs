using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Compression;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Decoder
{
    class PredictionSchemeDeltaDecoder : PredictionScheme
    {
        public override bool Initialized => true;
        public override PredictionSchemeMethod PredictionMethod => PredictionSchemeMethod.Difference;

        public PredictionSchemeDeltaDecoder(PointAttribute attribute, PredictionSchemeTransform transform)
            : base(attribute, transform)
        {

        }
        public override void ComputeCorrectionValues(Span<int> in_data, Span<int> out_corr, int size, int num_components, int[] entry_to_point_id_map)
        {
            throw new NotImplementedException();
        }

        public override void ComputeOriginalValues(Span<int> in_corr, Span<int> out_data, int size, int num_components,
            int[] entry_to_point_id_map)
        {

            this.transform_.InitializeDecoding(num_components);
            // Decode the original value for the first element.
            Span<int> zero_vals = stackalloc int[num_components];
            this.transform_.ComputeOriginalValue(zero_vals, in_corr, out_data);

            // Decode data from the front using D(i) = D(i) + D(i - 1).
            for (int i = num_components; i < size; i += num_components)
            {
                this.transform_.ComputeOriginalValue(out_data, i - num_components,
                    in_corr, i, out_data, i);
            }

        }
    }
}
