using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Encoder
{
    /// <summary>
    /// Abstract base class for all mesh encoders. It provides some basic
    /// funcionality that's shared between different encoders.
    /// </summary>
    abstract class MeshEncoder : PointCloudEncoder
    {
        private DracoMesh mesh;


        public override EncodedGeometryType GeometryType
        {
            get { return EncodedGeometryType.TriangularMesh; }
        }

        /// <summary>
        /// Returns the base connectivity of the encoded mesh (or nullptr if it is not
        /// initialized).
        /// </summary>
        public virtual CornerTable CornerTable
        {
            get { return null; }
        }

        /// <summary>
        /// Returns the attribute connectivity data or nullptr if it does not exist.
        /// </summary>
        public virtual MeshAttributeCornerTable GetAttributeCornerTable(int attId)
        {
            return null;
        }

        /// <summary>
        /// Returns the encoding data for a given attribute or nullptr when the data
        /// does not exist.
        /// </summary>
        public virtual MeshAttributeIndicesEncodingData GetAttributeEncodingData(int attId)
        {
            return null;
        }

        public DracoMesh Mesh
        {
            get { return mesh; }
            set { mesh = value;
                PointCloud = value;
            }
        }

        protected override void EncodeGeometryData()
        {
            EncodeConnectivity();
        }

        /// <summary>
        /// Needs to be implemented by the derived classes.
        /// </summary>
        protected abstract void EncodeConnectivity();

        // TODO(ostava): Prediction schemes need refactoring.
        /*
        // This method should be overriden by derived class to perform custom
        // initialization of various prediction schemes.
        virtual bool InitPredictionSchemeInternal(
            const MeshAttributeEncoder *attEncoder,
            PredictionSchemeInterface *scheme) {
          return true;
        }
        */
    }
}
