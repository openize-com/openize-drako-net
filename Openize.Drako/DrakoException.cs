using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Openize.Drako
{
    /// <summary>
    /// Exception when failed to encode or decode draco files.
    /// </summary>
    
#if DRACO_EMBED_MODE
    internal
#else
    public
#endif
        class DrakoException : Exception
    {
    }
}
