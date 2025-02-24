# Openize.Drako


Openize.Drako was ported from [Google Draco](https://github.com/google/draco).


### 1. Installation

```
dotnet add package Openize.Drako
```


### 2. Example

```csharp
using Openize.Drako;
Vector3[] controlPoints = new Vector3[]
{
        new Vector3( -5, 0, 5.0f),
        new Vector3( 5, 0, 5.0f),
        new Vector3( 5, 10, 5.0f),
        new Vector3( -5, 10, 5.0f),
        new Vector3( -5, 0, -5.0f),
        new Vector3( 50, 0, -5.0f),
        new Vector3( 5, 10, -5.0f),
        new Vector3( -5, 10, -5.0f)
};

int[] indices = new int[]
{
        0,1,2, 0, 2, 3, // Front face (Z+)
        1,5,6, 1, 6, 2, // Right side (X+)
        5,4,7, 5, 7, 6, // Back face (Z-)
        4,0,3, 4, 3, 7, // Left side (X-)
        0,4,5, 0, 5, 1, // Bottom face (Y-)
        3,2,6, 3, 6, 7 // Top face (Y+)
};

var mesh = new DracoMesh();
//construct an attribute for position, with type float[3], 
var attrPos = PointAttribute.Wrap(AttributeType.Position, controlPoints);
mesh.AddAttribute(attrPos);
//add triangle indices
mesh.Indices.AddRange(indices);
//number of the control points, it's required for the encoder to produce correct result.
mesh.NumPoints = 8;
//You can also use following methods to deduplicate the attributes to reduce the file size
//mesh.DeduplicateAttributeValues();
//mesh.DeduplicatePointIds();

var opt = new DracoEncodeOptions();
var drcBytes = Openize.Draco.Draco.Encode(mesh, opt);
```

## License
Openize.Drako is available under [Openize License](LICENSE).
> [!CAUTION]
> Openize does not and cannot grant You a patent license for the utilization of [Google Draco](https://github.com/google/draco) compression/decompression technologies.

## OSS Notice
Sample files used for tests and located in the "./Openize.Drako.Tests/TestsData" folder belong to [Google Draco](https://github.com/google/draco) and are used according to [Apache License 2.0](https://github.com/google/draco/blob/main/LICENSE)


## Coming updates
Openize.Drako will receive new features and regular updates to stay in sync with the latest versions of [Google Draco](https://github.com/google/draco). We appreciate your patience as we work on these improvements. Stay tuned for more updates soon.

