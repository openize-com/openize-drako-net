using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileFormat.Drako.Utils;

namespace FileFormat.Drako.Decoder
{
    /// <summary>
    /// Class for performing rANS decoding using a desired number of precision bits.
    /// The number of precision bits needs to be the same as with the RAnsEncoder
    /// that was used to encode the input data.
    /// </summary>
    class RAnsDecoder
    {

        private const int ioBase = 256;
        struct ransDecSym
        {
            internal uint val;
            internal uint prob;
            internal uint cumProb; // not-inclusive
        }

        struct ransSym
        {
            internal uint prob;
            internal uint cumProb; // not-inclusive
        }

        private int ransPrecisionBits;
        private int ransPrecision;
        private int lRansBase;
        private uint[] lutTable;
        private ransSym[] probabilityTable;

        private BytePointer buf;
        private int bufOffset;
        private uint state;

        public RAnsDecoder(int ransPrecisionBits)
        {
            this.ransPrecisionBits = ransPrecisionBits;
            this.ransPrecision = 1 << ransPrecisionBits;
            this.lRansBase = ransPrecision * 4;
        }

        /// <summary>
        /// Initializes the decoder from the input buffer. The |offset| specifies the
        /// number of bytes encoded by the encoder. A non zero return value is an
        /// error.
        /// </summary>
        public int readInit(BytePointer buf, int offset)
        {
            uint x;
            if (offset < 1)
                return 1;
            this.buf = buf;
            x = (uint)(buf[offset - 1] >> 6);
            if (x == 0)
            {
                this.bufOffset = offset - 1;
                this.state = (uint)(buf[offset - 1] & 0x3F);
            }
            else if (x == 1)
            {
                if (offset < 2)
                    return 1;
                this.bufOffset = offset - 2;
                this.state = (uint)buf.ToUInt16LE(offset - 2) & 0x3FFF;
            }
            else if (x == 2)
            {
                if (offset < 3)
                    return 1;
                this.bufOffset = offset - 3;
                this.state = (uint)buf.ToUInt24LE(offset - 3) & 0x3FFFFF;
            }
            else if (x == 3)
            {
                this.bufOffset = offset - 4;
                this.state = buf.ToUInt32LE(offset - 4) & 0x3FFFFFFF;
            }
            else
            {
                return 1;
            }
            this.state += (uint)lRansBase;
            if (this.state >= lRansBase * ioBase)
                return 1;
            return 0;
        }

        public bool ReadEnd()
        {
            return this.state == lRansBase;
        }

        public bool readerHasError()
        {
            return this.state < lRansBase && this.bufOffset == 0;
        }

        public int Read()
        {
            uint rem;
            uint quo;
            ransDecSym sym = new ransDecSym();
            while (this.state < lRansBase && this.bufOffset > 0)
            {
                this.state = (uint)(this.state * ioBase + this.buf[--this.bufOffset]);
            }
            // |ransPrecision| is a power of two compile time constant, and the below
            // division and modulo are going to be optimized by the compiler.
            quo = (uint)(this.state / ransPrecision);
            rem = (uint)(this.state % ransPrecision);
            fetchSym(ref sym, (int)rem);
            this.state = (uint)(quo * sym.prob + rem - sym.cumProb);
            return (int)sym.val;
        }

        /// <summary>
        /// Construct a look up table with |ransPrecision| number of entries.
        /// Returns false if the table couldn't be built (because of wrong input data).
        /// </summary>
        public bool BuildLookupTable(uint[] tokenProbs, int numSymbols)
        {
            lutTable = new uint[ransPrecision];
            probabilityTable = new ransSym[numSymbols];
            uint cumProb = 0;
            uint actProb = 0;
            for (int i = 0; i < numSymbols; ++i)
            {
                ransSym sym = new ransSym();
                sym.prob = tokenProbs[i];
                sym.cumProb = cumProb;
                probabilityTable[i] = sym;
                cumProb += tokenProbs[i];
                if (cumProb > ransPrecision)
                {
                    return false;
                }
                for (int j = (int)actProb; j < cumProb; ++j)
                {
                    lutTable[j] = (uint)i;
                }
                actProb = cumProb;
            }
            if (cumProb != ransPrecision)
            {
                return false;
            }
            return true;
        }

        private void fetchSym(ref ransDecSym res, int rem)
        {
            uint symbol = lutTable[rem];
            res.val = symbol;
            res.prob = probabilityTable[(int)symbol].prob;
            res.cumProb = probabilityTable[(int)symbol].cumProb;
        }

    }
}
