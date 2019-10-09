Begin["CMLConverter`"]

ImportExport`RegisterExport[
          "CML",
          CMLConverter`cmlExport,
          "FunctionChannels" -> {"FileNames"},
          "Sources" -> FileNameJoin[{DirectoryName[$InputFileName], "Converter.m"}]
]

End[]