using System;
using System.Collections.Generic;

namespace Openize.Drako.Utils
{
    /// <summary>
    /// Class that can be used to keep track of the Shannon entropy on streamed data.
    /// As new symbols are pushed to the tracker, the entropy is automatically
    /// recomputed. The class also support recomputing the entropy without actually
    /// pushing the symbols to the tracker through the Peek() method.
    /// </summary>
    public class ShannonEntropyTracker
    {
        /// <summary>
        /// Struct for holding entropy data about the symbols added to the tracker.
        /// It can be used to compute the number of bits needed to store the data using
        /// the method:
        ///   ShannonEntropyTracker.GetNumberOfDataBits(entropy_data);
        /// or to compute the approximate size of the frequency table needed by the
        /// rans coding using method:
        ///   ShannonEntropyTracker.GetNumberOfRAnsTableBits(entropy_data);
        /// </summary>
        public struct EntropyData
        {
            public double entropy_norm;
            public int num_values;
            public int max_symbol;
            public int num_unique_symbols;

            public EntropyData(double entropyNorm, int numValues, int maxSymbol, int numUniqueSymbols)
            {
                entropy_norm = entropyNorm;
                num_values = numValues;
                max_symbol = maxSymbol;
                num_unique_symbols = numUniqueSymbols;
            }

            public static EntropyData Default => new EntropyData(0.0, 0, 0, 0);
        }

        private List<int> frequencies_;
        private EntropyData entropy_data_;

        public ShannonEntropyTracker()
        {
            frequencies_ = new List<int>();
            entropy_data_ = EntropyData.Default;
        }

        /// <summary>
        /// Adds new symbols to the tracker and recomputes the entropy accordingly.
        /// </summary>
        public EntropyData Push(uint[] symbols, int num_symbols)
        {
            return UpdateSymbols(symbols, num_symbols, true);
        }

        /// <summary>
        /// Returns new entropy data for the tracker as if |symbols| were added to the
        /// tracker without actually changing the status of the tracker.
        /// </summary>
        public EntropyData Peek(uint[] symbols, int num_symbols)
        {
            return UpdateSymbols(symbols, num_symbols, false);
        }

        /// <summary>
        /// Gets the number of bits needed for encoding symbols added to the tracker.
        /// </summary>
        public long GetNumberOfDataBits()
        {
            return GetNumberOfDataBits(entropy_data_);
        }

        /// <summary>
        /// Gets the number of bits needed for encoding frequency table using the rans
        /// encoder.
        /// </summary>
        public long GetNumberOfRAnsTableBits()
        {
            return GetNumberOfRAnsTableBits(entropy_data_);
        }

        /// <summary>
        /// Gets the number of bits needed for encoding given |entropy_data|.
        /// </summary>
        public long GetNumberOfDataBits(EntropyData entropy_data)
        {
            if (entropy_data.num_values < 2)
            {
                return 0;
            }
            // We need to compute the number of bits required to represent the stream
            // using the entropy norm. Note that:
            //
            //   entropy = log2(num_values) - entropy_norm / num_values
            //
            // and number of bits required for the entropy is: num_values * entropy
            //
            return (long)Math.Ceiling(
                entropy_data.num_values * Math.Log(entropy_data.num_values, 2) -
                entropy_data.entropy_norm);
        }

        /// <summary>
        /// Gets the number of bits needed for encoding frequency table using the rans
        /// encoder for the given |entropy_data|.
        /// </summary>
        public long GetNumberOfRAnsTableBits(EntropyData entropy_data)
        {
            return ApproximateRAnsFrequencyTableBits(entropy_data.max_symbol + 1,
                                                     entropy_data.num_unique_symbols);
        }

        private EntropyData UpdateSymbols(uint[] symbols, int num_symbols, bool push_changes)
        {
            EntropyData ret_data = entropy_data_;
            ret_data.num_values += num_symbols;

            for (int i = 0; i < num_symbols; ++i)
            {
                uint symbol = symbols[i];
                
                // Ensure the frequencies list is large enough
                while (frequencies_.Count <= symbol)
                {
                    frequencies_.Add(0);
                }

                // Update the entropy of the stream. Note that entropy of |N| values
                // represented by |S| unique symbols is defined as:
                //
                //  entropy = -sum_over_S(symbol_frequency / N * log2(symbol_frequency / N))
                //
                // To avoid the need to recompute the entire sum when new values are added,
                // we can instead update a so called entropy norm that is defined as:
                //
                //  entropy_norm = sum_over_S(symbol_frequency * log2(symbol_frequency))
                //
                // In this case, all we need to do is update entries on the symbols where
                // the frequency actually changed.
                //
                // Note that entropy_norm and entropy can be easily transformed to the
                // actual entropy as:
                //
                //  entropy = log2(N) - entropy_norm / N
                //
                double old_symbol_entropy_norm = 0;
                int frequency = frequencies_[(int)symbol];
                
                if (frequency > 1)
                {
                    old_symbol_entropy_norm = frequency * Math.Log(frequency, 2);
                }
                else if (frequency == 0)
                {
                    ret_data.num_unique_symbols++;
                    if (symbol > ret_data.max_symbol)
                    {
                        ret_data.max_symbol = (int)symbol;
                    }
                }
                
                frequency++;
                frequencies_[(int)symbol] = frequency;
                double new_symbol_entropy_norm = frequency * Math.Log(frequency, 2);

                // Update the final entropy.
                ret_data.entropy_norm += new_symbol_entropy_norm - old_symbol_entropy_norm;
            }

            if (push_changes)
            {
                // Update entropy data of the stream.
                entropy_data_ = ret_data;
            }
            else
            {
                // We are only peeking so do not update the stream.
                // Revert changes in the frequency table.
                for (int i = 0; i < num_symbols; ++i)
                {
                    uint symbol = symbols[i];
                    frequencies_[(int)symbol]--;
                }
            }

            return ret_data;
        }

        /// <summary>
        /// Compute approximate frequency table size needed for storing the provided
        /// symbols.
        /// </summary>
        private static long ApproximateRAnsFrequencyTableBits(int max_value, int num_unique_symbols)
        {
            // Approximate number of bits for storing zero frequency entries using the
            // run length encoding (with max length of 64).
            long table_zero_frequency_bits =
                8 * (num_unique_symbols + (max_value - num_unique_symbols) / 64);
            return 8 * num_unique_symbols + table_zero_frequency_bits;
        }
    }
}
