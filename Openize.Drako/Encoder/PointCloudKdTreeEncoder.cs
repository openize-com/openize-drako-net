using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FileFormat.Drako.Encoder
{
    class PointCloudKdTreeEncoder : PointCloudEncoder
    {
        public override DracoEncodingMethod GetEncodingMethod()
        {
            return DracoEncodingMethod.KdTree;
        }

        protected override void GenerateAttributesEncoder(int attId)
        {

            if (NumAttributesEncoders == 0)
            {
                // Create a new attribute encoder only for the first attribute.
                AddAttributesEncoder(new KdTreeAttributesEncoder(attId));
                return ;
            }

            // Add a new attribute to the attribute encoder.
            AttributesEncoder(0).AddAttributeId(attId);
        }

        protected override void EncodeGeometryData()
        {
            int num_points = PointCloud.NumPoints;
            Buffer.Encode(num_points);
        }

    }
}
