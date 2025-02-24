using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Openize.Drako.Decoder
{
    class PointCloudSequentialDecoder : PointCloudDecoder
    {
        public PointCloudSequentialDecoder()
            : base(EncodedGeometryType.PointCloud)
        {

        }

        protected override void DecodeGeometryData()
        {
            int num_points = buffer.DecodeI32();
            PointCloud.NumPoints = num_points;
        }

        protected override void CreateAttributesDecoder(int attrDecoderId)
        {
            // Always create the basic attribute decoder.
            SetAttributesDecoder(
                attrDecoderId,
                new SequentialAttributeDecodersController(
                    new LinearSequencer(PointCloud.NumPoints)));
        }
    }
}
