using System.Runtime.InteropServices;
using System.Text;

namespace FileFormat.Drako.App
{
    class Args
    {
        public string InputFile { get; set; }
        public string OutputFile { get; set; }
    }

    internal class Program
    {
        static void Main(string[] args)
        {
            try
            {
                var cliArgs = ParseArgs(args);
                if (string.IsNullOrEmpty(cliArgs.InputFile))
                {
                    Help();
                    return;
                }
                ConvertToObj(cliArgs.InputFile, cliArgs.OutputFile);
            }
            catch (InvalidDataException ex)
            {
                Console.Error.WriteLine("ERROR:" + ex.Message);
            }
        }
        static Args? ParseArgs(string[] args)
        {
            var ret = new Args();
            for (int i = 0; i < args.Length; i++)
            {
                var arg = args[i];
                if (arg == "-i")
                {
                    if (i + 1 < args.Length)
                    {
                        ret.InputFile = args[++i];
                        if (!File.Exists(ret.InputFile))
                            throw new InvalidDataException($"Cannot find input file {ret.InputFile}");
                    }
                    else
                        throw new InvalidDataException("ERROR: input file is required after -i");
                }
                else if (arg == "-o")
                {
                    if (i + 1 < args.Length)
                        ret.OutputFile = args[++i];
                    else
                        throw new InvalidDataException("ERROR: output file is required after -o");
                }
            }
            return ret;
        }
        static void Help()
        {
            Console.WriteLine("  Convert a .drc file to .obj");
            Console.WriteLine("    FileFormat.Drako.App -i [input-file] -o [output-file]");
        }
        static void ConvertToObj(string inputFile, string outputFile)
        {
            if (string.IsNullOrEmpty(outputFile))
            {
                ConvertToObj(inputFile, Path.ChangeExtension(inputFile, ".obj"));
                return;
            }
            //load draco file
            var bytes = File.ReadAllBytes(inputFile);
            var mesh = Draco.Decode(bytes) as DracoMesh;
            if (mesh == null)
                throw new InvalidDataException("Input file is not a valid draco file.");
            var attrPos = mesh.GetNamedAttribute(AttributeType.Position);
            var points = MemoryMarshal.Cast<byte, float>(attrPos.Buffer.AsSpan());
            var sb = new StringBuilder();
            for (int i = 0; i < points.Length; i += 3)
            {
                sb.AppendLine($"v {points[i]} {points[i + 1]} {points[i + 2]}");
            }
            Span<int> face = stackalloc int[3];
            for (int i = 0; i < mesh.NumFaces; i++)
            {
                mesh.ReadFace(i, face);
                var a = attrPos.MappedIndex(face[0]) + 1;
                var b = attrPos.MappedIndex(face[1]) + 1;
                var c = attrPos.MappedIndex(face[2]) + 1;
                sb.AppendLine($"f {a} {b} {c}");
            }
            File.WriteAllText(outputFile, sb.ToString());
            Console.WriteLine($"File {inputFile} has been converted to {outputFile}.");
        }

    }
}
