# CMLConvert
### A package for importing and exporting Chemical Markup Language files.

CMLConvert extends the `Import[]`/`Export[]` functionality of *Mathematica* to handle [Chemical Markup Language (CML)](http://xml-cml.org/) files.

Download the paclet from the [releases](https://github.com/tpfto/CMLConvert/releases) page, and install it by evaluating the following in *Mathematica*:

    Needs["PacletManager`"]
    PacletInstall["/path/to/paclet/CMLConvert-1.0.paclet"];
    Get[FileNameJoin[{$UserBaseDirectory, "Paclets", "Repository", "CMLConvert-1.0", "CMLConvert.m"}]];

where `"/path/to/paclet/CMLConvert-1.0.paclet"` should be replaced with the actual path to the downloaded paclet file,
and `"CMLConvert-1.0"` should be replaced with the actual version to be installed.

After installation, try it out:

    (* for version 12 *)
    m = Import["http://chem-file.sourceforge.net/data/polycyclic_alkanes/abietic_acid.cml", "CML"];
    MoleculePlot3D[m]
    
    Export["abietic.cml", m, "CML"]
    
    (* for version 11 and earlier *)
    Import["http://chem-file.sourceforge.net/data/polycyclic_alkanes/abietic_acid.cml", "CML"]

    asp = Import["ExampleData/aspirin.mol",
                 {"MOL", {"EdgeRules", "EdgeTypes", "FormalCharges",
                          "MassNumbers", "VertexCoordinates", "VertexTypes"}}];
    Export["aspirin.cml", asp, {"CML", {"EdgeRules", "EdgeTypes", "FormalCharges",
                                        "MassNumbers", "VertexCoordinates", "VertexTypes"}}]
