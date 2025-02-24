using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using Openize.Drako.Utils;

namespace Openize.Drako.Decoder
{
    class SequentialQuantizationAttributeDecoder : SequentialIntegerAttributeDecoder
    {
        /// <summary>
        /// Max number of quantization bits used to encode each component of the
        /// attribute.
        /// </summary>
        private int quantizationBits;
        private float[] minValue;
        private float maxValueDif;

        public override void Initialize(PointCloudDecoder decoder, int attributeId)
        {
            base.Initialize(decoder, attributeId);
            PointAttribute attribute = decoder.PointCloud.Attribute(attributeId);
            // Currently we can quantize only floating point arguments.
            if (attribute.DataType != DataType.FLOAT32)
                throw DracoUtils.Failed();
        }

        public override void DecodeIntegerValues(int[] pointIds, DecoderBuffer inBuffer)
        {
            if (Decoder.BitstreamVersion < 20)
            {
                DecodeQuantizedDataInfo();
            }
            base.DecodeIntegerValues(pointIds, inBuffer);
        }

        protected override void StoreValues(int numValues)
        {
            DequantizeValues(numValues);
        }

        public override void DecodeDataNeededByPortableTransform(int[] pointIds, DecoderBuffer in_buffer)
        {
            if (Decoder.BitstreamVersion >= 20)
            {
                // Decode quantization data here only for files with bitstream version 2.0+
                DecodeQuantizedDataInfo();
            }
            // Store the decoded transform data in portable attribute;
            var transform = new AttributeQuantizationTransform();
            transform.SetParameters(quantizationBits, minValue, attribute.ComponentsCount, maxValueDif);
            transform.TransferToAttribute(PortableAttribute);
        }

        private void DecodeQuantizedDataInfo()
        {
            int numComponents = Attribute.ComponentsCount;
            minValue = new float[numComponents];
            if (!Decoder.Buffer.Decode(minValue))
                throw DracoUtils.Failed();
            maxValueDif = decoder.Buffer.DecodeF32();
            byte quantizationBits = Decoder.Buffer.DecodeU8();
            if (quantizationBits > 31)
                throw DracoUtils.Failed();
            this.quantizationBits = quantizationBits;
        }

        private void DequantizeValues(int numValues)
        {
            // Convert all quantized values back to floats.
            int maxQuantizedValue = (1 << (quantizationBits)) - 1;
            int numComponents = Attribute.ComponentsCount;
            int entrySize = sizeof(float) * numComponents;
            float[] attVal = new float[numComponents];
            int quantValId = 0;
            int outBytePos = 0;
            Dequantizer dequantizer = new Dequantizer(maxValueDif, maxQuantizedValue);
            //IntArray values = IntArray.Wrap(PortableAttribute.Buffer.GetBuffer(), 0, numValues * numComponents);
            var values = MemoryMarshal.Cast<byte, int>(PortableAttribute.Buffer.GetBuffer().AsSpan(0, numValues * numComponents * 4));
            for (uint i = 0; i < numValues; ++i)
            {
                for (int c = 0; c < numComponents; ++c)
                {
                    float value = dequantizer.DequantizeFloat(values[quantValId++]);
                    value = value + minValue[c];
                    attVal[c] = value;
                }
                // Store the floating point value into the attribute buffer.
                Attribute.Buffer.Write(outBytePos, attVal);
                outBytePos += entrySize;
            }
        }
    }
}
