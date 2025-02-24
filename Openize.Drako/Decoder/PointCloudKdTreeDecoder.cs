using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Decoder
{
    class PointCloudKdTreeDecoder : PointCloudDecoder
    {
        public PointCloudKdTreeDecoder() : base(EncodedGeometryType.PointCloud)
        {
        }

        protected override void DecodeGeometryData()
        {
            int num_points = buffer.DecodeI32();
            if (num_points < 0)
                throw DracoUtils.Failed();
            PointCloud.NumPoints = num_points;
        }

        protected override void CreateAttributesDecoder(int attrDecoderId)
        {
            // Always create the basic attribute decoder.
            SetAttributesDecoder(attrDecoderId, new KdTreeAttributesDecoder());
        }
    }
}
