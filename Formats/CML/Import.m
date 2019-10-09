Begin["CMLConverter`"]

If[(FileFormatDump`AddFormat["CML", False, False, False, False, False, {"*.cml"}, {"CHEMICAL/X-CML"}, None, {}]; $VersionNumber < 12.),
(* older versions *)
ImportExport`RegisterImport[
          "CML",
          CMLConverter`cmlImport,
          {
            "Graphics3D" :> CMLConverter`cmlImport3D,
            "StructureDiagram" :> ImportExport`StructureDiagram,
            Automatic :> CMLConverter`cmlToDefault
           },
           "AvailableElements" -> {"EdgeRules", "EdgeTypes", "FormalCharges", "Graphics3D",
                                               "MassNumbers", "StructureDiagram", "VertexCoordinates", "VertexTypes"},
           "DefaultElement" -> Automatic,
           "FunctionChannels" -> {"FileNames"},
           "Sources" -> FileNameJoin[{DirectoryName[$InputFileName], "Converter.m"}]
],
(* newer versions with Molecule *)
ImportExport`RegisterImport[
          "CML",
          {
            "Molecule" :> CMLConverter`cmlMolecule, 
            "StructureDiagram" :> CMLConverter`cmlImport2D,
            "Graphics3D" :> CMLConverter`cmlImport3D, 
            CMLConverter`cmlImport
            },
           "DefaultElement" -> "Molecule",
           "FunctionChannels" -> {"FileNames"},
           "Sources" -> FileNameJoin[{DirectoryName[$InputFileName], "Converter.m"}]]
]

End[]