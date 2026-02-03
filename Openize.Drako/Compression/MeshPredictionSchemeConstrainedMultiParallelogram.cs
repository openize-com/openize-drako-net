using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Decoder;
using Openize.Drako.Utils;

namespace Openize.Drako.Compression
{
    internal class MeshPredictionSchemeConstrainedMultiParallelogram : MeshPredictionScheme
    {

        const int OPTIMAL_MULTI_PARALLELOGRAM = 0;
        const int kMaxNumParallelograms = 4;
        // Crease edges are used to store whether any given edge should be used for
        // parallelogram prediction or not. New values are added in the order in which
        // the edges are processed. For better compression, the flags are stored in
        // in separate contexts based on the number of available parallelograms at a
        // given vertex.
        List<bool>[] is_crease_edge_ = new List<bool>[kMaxNumParallelograms];
        
        public MeshPredictionSchemeConstrainedMultiParallelogram(PointAttribute attribute,
            PredictionSchemeTransform transform,
            MeshPredictionSchemeData meshData)
            :base(attribute, transform, meshData)
        {
            for (int i = 0; i < kMaxNumParallelograms; ++i)
            {
                is_crease_edge_[i] = new List<bool>();
            }
        }

        public override PredictionSchemeMethod PredictionMethod => PredictionSchemeMethod.ConstrainedMultiParallelogram;

        struct PredictionConfiguration {
            public Error error;
            public byte configuration;  // Bitfield, 1 use parallelogram, 0 don't use it.
            public int num_used_parallelograms;
            public int[] predicted_value;
            public int[] residuals;
        }

        uint[] entropy_symbols_;
        ShannonEntropyTracker entropy_tracker_ = new ShannonEntropyTracker();
        
        public override void ComputeCorrectionValues(Span<int> in_data, Span<int> out_corr, int size, int num_components, int[] entry_to_point_id_map)
        {
            this.transform_.InitializeEncoding(in_data,  num_components);
            var table = this.meshData.CornerTable;
            var vertex_to_data_map =
                this.meshData.vertexToDataMap;

            // Predicted values for all simple parallelograms encountered at any given
            // vertex.
            int[][] pred_vals = new int[kMaxNumParallelograms] [];
            for (int i = 0; i < kMaxNumParallelograms; ++i) {
                pred_vals[i] = new int[num_components];
            }
            // Used to store predicted value for various multi-parallelogram predictions
            // (combinations of simple parallelogram predictions).
            Span<int> multi_pred_vals = stackalloc int[num_components];
            entropy_symbols_ = new uint[num_components];

            // Bit-field used for computing permutations of excluded edges
            // (parallelograms).
            Span<bool> exluded_parallelograms = stackalloc bool[kMaxNumParallelograms];

            // Data about the number of used parallelogram and total number of available
            // parallelogram for each context. Used to compute overhead needed for storing
            // the parallelogram choices made by the encoder.
            long[] total_used_parallelograms = new long[kMaxNumParallelograms];
            long[] total_parallelograms = new long[kMaxNumParallelograms];

            Span<int> current_residuals = stackalloc int[num_components];

            // We start processing the vertices from the end because this prediction uses
            // data from previous entries that could be overwritten when an entry is
            // processed.
            for (int p =
                     this.meshData.dataToCornerMap.Count - 1;
                 p > 0; --p) {
                var start_corner_id =
                    this.meshData.dataToCornerMap[p];

                // Go over all corners attached to the vertex and compute the predicted
                // value from the parallelograms defined by their opposite faces.
                var corner_id = start_corner_id;
                int num_parallelograms = 0;
                bool first_pass = true;
                while (corner_id != CornerTable.kInvalidCornerIndex) {
                    if (MeshPredictionSchemeParallelogram.ComputeParallelogramPrediction(
                            p, corner_id, table, vertex_to_data_map, in_data, num_components,
                            pred_vals[num_parallelograms])) {
                        // Parallelogram prediction applied and stored in
                        // |pred_vals[num_parallelograms]|
                        ++num_parallelograms;
                        // Stop processing when we reach the maximum number of allowed
                        // parallelograms.
                        if (num_parallelograms == kMaxNumParallelograms) {
                            break;
                        }
                    }

                    // Proceed to the next corner attached to the vertex. First swing left
                    // and if we reach a boundary, swing right from the start corner.
                    if (first_pass) {
                        corner_id = table.SwingLeft(corner_id);
                    } else {
                        corner_id = table.SwingRight(corner_id);
                    }
                    if (corner_id == start_corner_id) {
                        break;
                    }
                    if (corner_id == CornerTable.kInvalidCornerIndex && first_pass) {
                        first_pass = false;
                        corner_id = table.SwingRight(start_corner_id);
                    }
                }

                // Offset to the target (destination) vertex.
                int dst_offset = p * num_components;
                Error error;

                // Compute all prediction errors for all possible configurations of
                // available parallelograms.

                // Variable for holding the best configuration that has been found so far.
                PredictionConfiguration best_prediction = new PredictionConfiguration();

                // Compute delta coding error (configuration when no parallelogram is
                // selected).
                int src_offset = (p - 1) * num_components;
                error = ComputeError(in_data.Slice(src_offset), in_data.Slice(dst_offset),
                                     current_residuals, num_components);

                if (num_parallelograms > 0) {
                    total_parallelograms[num_parallelograms - 1] += num_parallelograms;
                    var new_overhead_bits =
                        ComputeOverheadBits(total_used_parallelograms[num_parallelograms - 1],
                                            total_parallelograms[num_parallelograms - 1]);
                    error.num_bits += new_overhead_bits;
                }

                best_prediction.error = error;
                best_prediction.configuration = 0;
                best_prediction.num_used_parallelograms = 0;
                best_prediction.predicted_value = new int[num_components];
                Array.Copy(in_data.ToArray(), src_offset, best_prediction.predicted_value, 0, num_components);
                best_prediction.residuals = new int[num_components];
                Array.Copy(current_residuals.ToArray(), 0, best_prediction.residuals, 0, num_components);

                // Compute prediction error for different cases of used parallelograms.
                for (int num_used_parallelograms = 1;
                     num_used_parallelograms <= num_parallelograms;
                     ++num_used_parallelograms) {
                    // Mark all parallelograms as excluded.
                    for (int j = 0; j < num_parallelograms; ++j)
                    {
                        exluded_parallelograms[j] = true;
                    }
                    // Mark the first |num_used_parallelograms| as not excluded.
                    for (int j = 0; j < num_used_parallelograms; ++j) {
                        exluded_parallelograms[j] = false;
                    }
                    // Permute over the excluded edges and compute error for each
                    // configuration (permutation of excluded parallelograms).
                    do {
                        // Reset the multi-parallelogram predicted values.
                        for (int j = 0; j < num_components; ++j) {
                            multi_pred_vals[j] = 0;
                        }
                        int configuration = 0;
                        for (int j = 0; j < num_parallelograms; ++j) {
                            if (exluded_parallelograms[j]) {
                                continue;
                            }
                            for (int c = 0; c < num_components; ++c) {
                                multi_pred_vals[c] += pred_vals[j][c];
                            }
                            // Set jth bit of the configuration.
                            configuration |= (1 << j);
                        }

                        for (int j = 0; j < num_components; ++j) {
                            multi_pred_vals[j] /= num_used_parallelograms;
                        }
                        error = ComputeError(multi_pred_vals, in_data.Slice(dst_offset),
                                             current_residuals, num_components);
                        if (num_parallelograms > 0) {
                            var new_overhead_bits = ComputeOverheadBits(
                                total_used_parallelograms[num_parallelograms - 1] +
                                    num_used_parallelograms,
                                total_parallelograms[num_parallelograms - 1]);

                            // Add overhead bits to the total error.
                            error.num_bits += new_overhead_bits;
                        }
                        if (error < best_prediction.error) {
                            best_prediction.error = error;
                            best_prediction.configuration = (byte)configuration;
                            best_prediction.num_used_parallelograms = num_used_parallelograms;
                            best_prediction.predicted_value = new int[num_components];
                            for (int j = 0; j < num_components; ++j)
                            {
                                best_prediction.predicted_value[j] = multi_pred_vals[j];
                            }
                            best_prediction.residuals = new int[num_components];
                            Array.Copy(current_residuals.ToArray(), 0, best_prediction.residuals, 0, num_components);
                        }
                    } while (NextPermutation(exluded_parallelograms, num_parallelograms));
                }
                if (num_parallelograms > 0) {
                    total_used_parallelograms[num_parallelograms - 1] +=
                        best_prediction.num_used_parallelograms;
                }

                // Update the entropy stream by adding selected residuals as symbols to the
                // stream.
                for (int i = 0; i < num_components; ++i) {
                    entropy_symbols_[i] =
                        ConvertSignedIntToSymbol(best_prediction.residuals[i]);
                }
                entropy_tracker_.Push(entropy_symbols_, num_components);

                for (int i = 0; i < num_parallelograms; ++i) {
                    if ((best_prediction.configuration & (1 << i)) == 0) {
                        // Parallelogram not used, mark the edge as crease.
                        is_crease_edge_[num_parallelograms - 1].Add(true);
                    } else {
                        // Parallelogram used. Add it to the predicted value and mark the
                        // edge as not a crease.
                        is_crease_edge_[num_parallelograms - 1].Add(false);
                    }
                }
                this.transform_.ComputeCorrection(in_data.Slice(dst_offset),
                                                    best_prediction.predicted_value,
                                                    out_corr.Slice(dst_offset), 0);
            }
            // First element is always fixed because it cannot be predicted.
            for (int i = 0; i < num_components; ++i) {
                pred_vals[0][i] = 0;
            }
            this.transform_.ComputeCorrection(in_data, 0, pred_vals[0], 0, out_corr, 0, 0);
        }

        /// <summary>
        /// Function used to compute number of bits needed to store overhead of the
        /// predictor. In this case, we consider overhead to be all bits that mark
        /// whether a parallelogram should be used for prediction or not. The input
        /// to this method is the total number of parallelograms that were evaluated so
        /// far(total_parallelogram), and the number of parallelograms we decided to
        /// use for prediction (total_used_parallelograms).
        /// Returns number of bits required to store the overhead.
        /// </summary>
        long ComputeOverheadBits(long total_used_parallelograms,
                                    long total_parallelogram)
        {
            // For now we assume RAns coding for the bits where the total required size
            // is directly correlated to the binary entropy of the input stream.
            // TODO(ostava): This should be generalized in case we use other binary
            // coding scheme.
            double entropy = ComputeBinaryShannonEntropy(
                (uint)(total_parallelogram),
                (uint)(total_used_parallelograms));

            // Round up to the nearest full bit.
            return (long)(
                Math.Ceiling((double)(total_parallelogram) * entropy));
        }

        double ComputeBinaryShannonEntropy(uint num_values,
                                           uint num_true_values)
        {
            if (num_values == 0)
            {
                return 0;
            }

            // We can exit early if the data set has 0 entropy.
            if (num_true_values == 0 || num_values == num_true_values)
            {
                return 0;
            }
            double true_freq =
                (double)(num_true_values) / (double)(num_values);
            double false_freq = 1.0 - true_freq;
            return -(true_freq * Math.Log(true_freq, 2) +
                     false_freq * Math.Log(false_freq, 2));
        }
        
        // Computes error for predicting |predicted_val| instead of |actual_val|.
        // Error is computed as the number of bits needed to encode the difference
        // between the values.
        Error ComputeError(Span<int> predicted_val,
                           Span<int> actual_val, Span<int> out_residuals,
                           int num_components)
        {
            Error error = new Error();

            for (int i = 0; i < num_components; ++i)
            {
                int dif = (predicted_val[i] - actual_val[i]);
                error.residual_error += Math.Abs(dif);
                out_residuals[i] = dif;
                // Entropy needs unsigned symbols, so convert the signed difference to an
                // unsigned symbol.
                entropy_symbols_[i] = ConvertSignedIntToSymbol(dif);
            }

            // Generate entropy data for case that this configuration was used.
            // Note that the entropy stream is NOT updated in this case.
            var entropy_data =
                entropy_tracker_.Peek(entropy_symbols_, num_components);

            error.num_bits = entropy_tracker_.GetNumberOfDataBits(entropy_data) +
                             entropy_tracker_.GetNumberOfRAnsTableBits(entropy_data);
            return error;
        }

        /// <summary>
        /// Helper function that converts a single signed integer value into an unsigned
        /// integer symbol that can be encoded using an entropy encoder.
        /// </summary>
        uint ConvertSignedIntToSymbol(int val)
        {
            // Early exit if val is positive.
            if (val >= 0)
            {
                return (uint)(val) << 1;
            }
            val = -(val + 1);  // Map -1 to 0, -2 to -1, etc..
            uint ret = (uint)(val);
            ret <<= 1;
            ret |= 1;
            return ret;
        }
        
        /// <summary>
        /// Generates the next permutation of the boolean array in lexicographic order.
        /// Returns true if a next permutation exists, false otherwise.
        /// </summary>
        private bool NextPermutation(Span<bool> array, int length)
        {
            // Find the largest index i such that array[i] < array[i + 1]
            // Convert bool to int for comparison (false=0, true=1)
            int i = length - 2;
            while (i >= 0 && (!array[i] || array[i + 1]))
            {
                i--;
            }

            if (i < 0)
            {
                return false; // No next permutation
            }

            // Find the largest index j such that array[i] < array[j]
            int j = length - 1;
            while (!array[j] || array[i])
            {
                j--;
            }

            // Swap array[i] and array[j]
            bool temp = array[i];
            array[i] = array[j];
            array[j] = temp;

            // Reverse the suffix starting at array[i + 1]
            int left = i + 1;
            int right = length - 1;
            while (left < right)
            {
                temp = array[left];
                array[left] = array[right];
                array[right] = temp;
                left++;
                right--;
            }

            return true;
        }
        
        public override void DecodePredictionData(DecoderBuffer buffer)
        {

            if (buffer.BitstreamVersion < 22)
            {
                // Decode prediction mode.
                var mode = buffer.DecodeU8();

                if (mode != OPTIMAL_MULTI_PARALLELOGRAM)
                {
                    // Unsupported mode.
                    throw DracoUtils.Failed();
                }
            }

            // Encode selected edges using separate rans bit coder for each context.
            for (int i = 0; i < kMaxNumParallelograms; ++i)
            {
                var num_flags = buffer.DecodeVarintU32();
                if (num_flags > this.meshData.CornerTable.NumCorners)
                {
                    throw DracoUtils.Failed();
                }
                if (num_flags > 0)
                {
                    is_crease_edge_[i] = new List<bool>((int)num_flags);
                    RAnsBitDecoder decoder = new RAnsBitDecoder();
                    decoder.StartDecoding(buffer);
                    for (var j = 0; j < num_flags; ++j)
                    {
                        is_crease_edge_[i].Add(decoder.DecodeNextBit());
                    }
                    decoder.EndDecoding();
                }
            }
            base.DecodePredictionData(buffer);
        }
        
        public override void ComputeOriginalValues(Span<int> in_corr, Span<int> out_data, int size, int num_components, int[] entry_to_point_id_map)
        {
            this.transform_.InitializeDecoding(num_components);

            // Predicted values for all simple parallelograms encountered at any given
            // vertex.
            int[][] pred_vals = new int[kMaxNumParallelograms][];
            for (int i = 0; i < kMaxNumParallelograms; ++i)
            {
                pred_vals[i] = new int[num_components];
            }
            this.transform_.ComputeOriginalValue(pred_vals[0], in_corr,
                                                   out_data);

            var table = this.meshData.CornerTable;
            var vertex_to_data_map = this.meshData.vertexToDataMap;

            // Current position in the |is_crease_edge_| array for each context.
            var is_crease_edge_pos = new int[kMaxNumParallelograms];

            // Used to store predicted value for multi-parallelogram prediction.
            var multi_pred_vals = new int[num_components];

            int corner_map_size =
                this.meshData.dataToCornerMap.Count;
            for (int p = 1; p < corner_map_size; ++p)
            {
                var start_corner_id =
                    this.meshData.dataToCornerMap[p];

                var corner_id = start_corner_id;
                int num_parallelograms = 0;
                bool first_pass = true;
                while (corner_id != CornerTable.kInvalidCornerIndex)
                {
                    if (MeshPredictionSchemeParallelogram.ComputeParallelogramPrediction(
                            p, corner_id, table, vertex_to_data_map, out_data,
                            num_components, pred_vals[num_parallelograms]))
                    {
                        // Parallelogram prediction applied and stored in
                        // |pred_vals[num_parallelograms]|
                        ++num_parallelograms;
                        // Stop processing when we reach the maximum number of allowed
                        // parallelograms.
                        if (num_parallelograms == kMaxNumParallelograms)
                        {
                            break;
                        }
                    }

                    // Proceed to the next corner attached to the vertex. First swing left
                    // and if we reach a boundary, swing right from the start corner.
                    if (first_pass)
                    {
                        corner_id = table.SwingLeft(corner_id);
                    }
                    else
                    {
                        corner_id = table.SwingRight(corner_id);
                    }
                    if (corner_id == start_corner_id)
                    {
                        break;
                    }
                    if (corner_id == CornerTable.kInvalidCornerIndex && first_pass)
                    {
                        first_pass = false;
                        corner_id = table.SwingRight(start_corner_id);
                    }
                }

                // Check which of the available parallelograms are actually used and compute
                // the final predicted value.
                int num_used_parallelograms = 0;
                if (num_parallelograms > 0)
                {
                    for (int i = 0; i < num_components; ++i)
                    {
                        multi_pred_vals[i] = 0;
                    }
                    // Check which parallelograms are actually used.
                    for (int i = 0; i < num_parallelograms; ++i)
                    {
                        int context = num_parallelograms - 1;
                        int pos = is_crease_edge_pos[context]++;
                        if (is_crease_edge_[context].Count <= pos)
                        {
                            throw DracoUtils.Failed();
                        }
                        bool is_crease = is_crease_edge_[context][pos];
                        if (!is_crease)
                        {
                            ++num_used_parallelograms;
                            for (int j = 0; j < num_components; ++j)
                            {
                                multi_pred_vals[j] = (int)((ulong)multi_pred_vals[j] + (ulong)pred_vals[i][j]);
                            }
                        }
                    }
                }
                int dst_offset = p * num_components;
                if (num_used_parallelograms == 0)
                {
                    // No parallelogram was valid.
                    // We use the last decoded point as a reference.
                    int src_offset = (p - 1) * num_components;
                    this.transform_.ComputeOriginalValue(
                        out_data.Slice(src_offset), in_corr.Slice(dst_offset), out_data.Slice(dst_offset));
                }
                else
                {
                    // Compute the correction from the predicted value.
                    for (int c = 0; c < num_components; ++c)
                    {
                        multi_pred_vals[c] /= num_used_parallelograms;
                    }
                    this.transform_.ComputeOriginalValue(
                        multi_pred_vals, in_corr.Slice(dst_offset), out_data.Slice(dst_offset));
                }
            }
        }
    }

    /// <summary>
    /// Struct for holding error information about prediction.
    /// </summary>
    struct Error
    {
        public long residual_error;
        public long num_bits;

        public static bool operator <(Error a, Error b)
        {
            if (a.num_bits != b.num_bits)
                return a.num_bits < b.num_bits;
            return a.residual_error < b.residual_error;
        }

        public static bool operator >(Error a, Error b)
        {
            if (a.num_bits != b.num_bits)
                return a.num_bits > b.num_bits;
            return a.residual_error > b.residual_error;
        }
    }
}
