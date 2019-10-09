
(* :Title: CML Converter *)

(* :Author: J. M. *)

(* :Summary:

     This package extends Import and Export to handle CML (Chemical Markup Language) files.
     For more information on the CML format, please visit http://xml-cml.org/ or http://cml.sourceforge.net/.

 *)

(* :Copyright:

     © 2017-2019 by J. M. (pleasureoffiguring(AT)gmail(DOT)com)
     This work is free. It comes without any warranty, to the extent permitted by applicable law.
     You can redistribute and/or modify it under the terms of the MIT License.
 *)

(* :Mathematica Version: 8.0 *)

(* :History:

     1.0 - initial release
*)

(* :Keywords:
     CML, Chemical Markup Language, export, files, import *)

(* :Limitations:
     currently not able to handle CML files with crystal, spectra, or reaction information.
     currently not able to handle CML files in array mode.
     currently does not recognize wedge/hash bonds as stereochemistry indicators
*)

(* :Warnings:
     adds definitions to Import and Export to enable processing of CML files.
*)

FileFormatDump`AddFormat["CML", False, False, False, False, False, {"*.cml"}, {"CHEMICAL/X-CML"}, None, {}];

Begin["CMLConverter`"]

Unprotect[cmlImport3D, cmlImport, cmlToDefault, cmlExport, cmlIndent];

$coordTol = 1.*^-4;

If[$VersionNumber < 12.,

(* for older versions *)

cmlImport[fn_String, opts___] :=
     Module[{atoms, atomList, atomData, bonds, bondData, bondList, bondType,
                  coords2d, coords3d, idx, mList, nAtoms, qList, molObj, vc, xmlObj},

                 xmlObj = Import[fn, {"XML", "XMLObject"}];
                 If[xmlObj === $Failed, Return[$Failed, Module]];

                 molObj = Flatten[Cases[xmlObj, XMLElement["molecule", tag_, rest__] :> rest, Infinity]];
                 If[molObj === {}, Return[$Failed, Module]];

                 atomList = Flatten[Cases[molObj, XMLElement["atomArray", _, rest__] :> rest, Infinity]];
                 If[atomList === {}, Return[$Failed, Module]];

                 bondList = Flatten[Cases[molObj, XMLElement["bondArray", _, rest__] :> rest, Infinity]];

                 atomData = Cases[atomList, XMLElement["atom", data_, rest__] :> Append[data, "atomTag" -> rest]];
                 atomData = SortBy[atomData, StringCases["id" /. #, s : (NumberString ..) :> FromDigits[s]] &];
                 idx = MapIndexed[#1 -> #2[[1]] &, "id" /. atomData];

                 atoms = "elementType" /. atomData; nAtoms = Length[atoms];

                 qList = If[FreeQ[atomData, "formalCharge"], Table[0, {nAtoms}],
                                "formalCharge" /. atomData /. {s_String /; StringMatchQ[s, NumberString] :> ToExpression[s],
                                                                                "formalCharge" -> 0}];

                 mList = Which[FreeQ[atomData, "isotopeNumber" | "isotope"], Table[None, {nAtoms}],
                                       ! FreeQ[atomData, "isotope"], "isotope" /. atomData,
                                       ! FreeQ[atomData, "isotopeNumber"], 
                                       "isotopeNumber" /. atomData] /.
                             {s_String /; StringMatchQ[s, NumberString] :> ToExpression[s], "isotope" | "isotopeNumber" -> None};
  
                 coords2d = Map[If[StringQ[#], Internal`StringToDouble[#], Null] &,
                                          ({"x2", "y2"} /. atomData) /. Thread[{"x2", "y2"} -> Null], {2}];
                 If[! FreeQ[coords2d, Null], coords2d = {}];

                 coords3d = Map[If[StringQ[#], Internal`StringToDouble[#], Null] &,
                                          ({"x3", "y3", "z3"} /. atomData) /. Thread[{"x3", "y3", "z3"} -> Null], {2}];
                 If[FreeQ[coords3d, Null], coords3d *= 100, coords3d = {}];

                 vc = Which[Length[coords3d] > 0 && MatrixQ[coords3d, NumberQ], 
                                   coords3d, Length[coords2d] > 0 && MatrixQ[coords2d, NumberQ], 
                                   coords2d,
                                   True, Missing["NotAvailable"]];

                 bondData = Cases[bondList, XMLElement["bond", data_, rest__] :> data];
                 If[bondData =!= {},
                     bonds = Rule @@@ (Map[StringSplit, "atomRefs2" /. bondData] /. idx);
                     bondType = Replace["order" /. bondData,
                                                    Thread[{"S" | "1", "D" | "2", "T" | "3", "A"} -> {"Single", "Double", "Triple", "Aromatic"}], {1}],
                     bonds = bondType = {}];
  
                 Thread[{"EdgeRules", "EdgeTypes", "FormalCharges", "MassNumbers", "VertexCoordinates", "VertexTypes"} ->
                            {bonds, bondType, qList, mList, vc, atoms}]];

cmlImport3D[rules_, opts___] := Block[{res},
     res = DataPaclets`ChemicalDataDump`MoleculePlot3D[rules, "MultiBond",
                                                                                     Sequence @@ Flatten[{opts, "ImageSizeScaled" -> True}]];
     If[Head[res] =!= Graphics3D, res = ImportExport`MoleculePlot3D[rules, opts]];
     res];

cmlToDefault[elems_, opts___?OptionQ] := Module[{coords},
     coords = "VertexCoordinates" /. elems;
     If[MatchQ[coords, Missing["NotAvailable"] | $Failed], Return[coords, Module]];
     If[! MatrixQ[coords], Return[$Failed, Module]];
     Switch[Last[Dimensions[coords]],
                2, ImportExport`StructureDiagram[elems, opts],
                3, cmlImport3D[elems, Sequence @@ Flatten[{opts, "InferBonds" -> False}]],
                True, $Failed]];
,
(* for newer versions *)

cmlImportMolecule[fn_String, opts___] := 
     Module[{atoms, atomList, atomData, atomStereo, bonds, bondData, bondList, bondStereo,
                   coords2d, coords3d, dbnd, idx, mList, nAtoms, qList, molObj, tets, xmlObj},

                 xmlObj = Import[fn, { "XML", "XMLObject"}];
                 If[xmlObj === $Failed, Return[$Failed, Module]];

                 molObj = Flatten[Cases[xmlObj, XMLElement["molecule", tag_, rest__] :> rest, Infinity]];
                 If[molObj === {}, Return[$Failed, Module]];

                 atomList = Flatten[Cases[molObj, XMLElement["atomArray", _, rest__] :> rest, Infinity]];
                 If[atomList === {}, Return[$Failed, Module]];

                 bondList = Flatten[Cases[molObj, XMLElement["bondArray", _, rest__] :> rest, Infinity]];
  
                 atomData = Cases[atomList, XMLElement["atom", data_, rest__] :> Append[data, "atomTag" -> rest]];
                 atomData = SortBy[atomData, StringCases["id" /. #, s : (NumberString ..) :> FromDigits[s]] &];
                 nAtoms = Length[atomData]; idx = MapIndexed[#1 -> #2[[1]] &, "id" /. atomData];

                 qList = Table[Null, {nAtoms}];
                 If[! FreeQ[atomData, "formalCharge"],
                     qList = "formalCharge" /. atomData /. {s_String /; StringMatchQ[s, NumberString] :> ToExpression[s], "formalCharge" -> Null}];
                 qList = Thread["FormalCharge" -> qList];

                 mList = Which[FreeQ[atomData, "isotope" | "isotopeNumber"], Table[Null, {nAtoms}],
                                       ! FreeQ[atomData, "isotope"], "isotope" /. atomData,
                                       ! FreeQ[atomData, "isotopeNumber"], "isotopeNumber" /. atomData];
                 mList = Thread["MassNumber" -> (mList /. {s_String /; StringMatchQ[s, NumberString] :> ToExpression[s],
                                                                                    "isotope" | "isotopeNumber" -> Null})];

                 atoms = MapThread[Atom[#1, ##2] &, {"elementType" /. atomData, qList, mList}] /. HoldPattern[_ -> Null] -> Unevaluated[];

                 coords2d = Map[If[StringQ[#], Internal`StringToDouble[#], Null] &,
                                          ({"x2", "y2"} /. atomData) /.  Thread[{"x2", "y2"} -> Null], {2}];
                 If[! FreeQ[coords2d, Null], coords2d = {}];

                 coords3d = Map[If[StringQ[#], Internal`StringToDouble[#], Null] &,
                                          ({"x3", "y3", "z3"} /. atomData) /. Thread[{"x3", "y3", "z3"} -> Null], {2}];
                 coords3d = If[FreeQ[coords3d, Null], QuantityArray[coords3d, "Angstroms"], {}];
  
                 atomStereo = Select[atomData, With[{t = "atomTag" /. #}, Length[t] != 0 && MemberQ[t, XMLElement["atomParity", __]]] &];
                 If[atomStereo =!= {},
                     atomStereo = {"id" /. atomStereo, Cases["atomTag" /. atomStereo, XMLElement["atomParity", rest__] :> Join[rest], Infinity]};

                     tets = {atomStereo[[1]], StringSplit /@ ("atomRefs4" /. Take[atomStereo[[2]], All, 1]), ToExpression /@ atomStereo[[2, All, -1]]} /. idx;
                     tets = MapThread[If[FreeQ[#2, #1], {#1, First[#2], Rest[#2], #3},
                                                    With[{tp = Position[#2, #1][[1, 1]]},
                                                            Join[{#1}, Through[{First, Rest}[Delete[#2, tp]]], {(1 - 2 Mod[tp, 2]) #3}]]] &, tets];
                     tets = MapAt[(# /. {1 -> "Clockwise", -1 -> "Counterclockwise"}) &, tets, {All, 4}];

                     atomStereo = Association[Prepend[Thread[{"ChiralCenter", "FiducialAtom", "Ligands", "Direction"} -> #], "StereoType" -> "Tetrahedral"]] & /@ tets];

                 If[bondList =!= {},
                     bondData = Cases[bondList, XMLElement["bond", data_, rest__] :> data];
                     bonds = Apply[Bond[StringSplit[#1] /. idx, #2 /. Thread[{"S" | "1", "D" | "2", "T" | "3", "A"} -> {"Single", "Double", "Triple", "Aromatic"}]] &,
                                           {"atomRefs2", "order"} /. bondData, {1}];

                     bondStereo = Cases[bondList, XMLElement["bondStereo", rest__] :> Flatten[{rest}], Infinity];
                     If[bondStereo =!= {} && ! FreeQ[bondStereo, "atomRefs4"],
                         dbnd = Map[StringSplit, "atomRefs4" /. Take[bondStereo, All, 1]][[All, {2, 3, 1, 4}]] /. idx;
                         dbnd = Transpose[{Take[dbnd, All, 2], Take[dbnd, All, -2], bondStereo[[All, -1]] /. {"C" -> "Together", "T" -> "Opposite"}}];

                         bondStereo = Association[Prepend[Thread[{"StereoBond", "Ligands", "Value"} -> #], "StereoType" -> "DoubleBond"]] & /@ dbnd,
                         bondStereo = {}],
                     bonds = bondStereo = {}];
  
                 Molecule @@ ({atoms, bonds, AtomCoordinates -> coords3d, AtomDiagramCoordinates -> coords2d,
                                        StereochemistryElements -> Join[atomStereo, bondStereo]} /. HoldPattern[_ -> {}] -> Unevaluated[])];

cmlMolecule[fileName_String, opts___] := ("Molecule" -> cmlImportMolecule[fileName]);

cmlImport[fileName_String, opts___] := Module[{mol = cmlImportMolecule[fileName], mop2, mop3, vc},

     mop2 = Options[mol, AtomDiagramCoordinates];
     If[mop2 === {}, 
         mop2 = {AtomDiagramCoordinates -> Missing["NotAvailable"]}];

     mop3 = Options[mol, AtomCoordinates];
     If[mop3 === {}, mop3 = {AtomCoordinates -> Missing["NotAvailable"]}];

     vc = {AtomCoordinates, AtomDiagramCoordinates} /. Join[mop3, mop2];
     vc = Switch[vc,
                       {_?MatrixQ, _?MatrixQ | _Missing}, First[vc],
                       {_Missing, _?MatrixQ}, Last[vc],
                       _, Missing["NotAvailable"]];

     {"EdgeRules" -> MoleculeValue[mol, "EdgeRules"], "EdgeTypes" -> MoleculeValue[mol, "EdgeTypes"], "FormalCharges" -> MoleculeValue[mol, "FormalCharges"],
       "MassNumbers" -> MoleculeValue[mol, "MassNumbers"], "VertexCoordinates" -> vc, "VertexTypes" -> MoleculeValue[mol, "VertexTypes"]}]

cmlImport2D[fileName_String, opts___] :=
     ("StructureDiagram" -> MoleculePlot[cmlImportMolecule[fileName], Sequence @@ FilterRules[{opts}, Options[MoleculePlot]]]);

cmlImport3D[fileName_String, opts___] :=
     ("Graphics3D" -> MoleculePlot3D[cmlImportMolecule[fileName],
                                                       Sequence @@ FilterRules[{opts} /. HoldPattern["Rendering" -> val_] :> (PlotTheme -> val), Options[MoleculePlot3D]]]);
]

cmlIndent[s_] := If[s === {"http://www.xml-cml.org/schema", "cml"}, False, Automatic]

cmlExport[fn_String, data : {(_Rule | _RuleDelayed) ..}, opts___] := 
     Module[{atoms, atomArray, bonds, bondArray, bondTypes, coords,
                   mList, nAtoms, nBonds, qList, sch, xmlObj},

                 {atoms, qList, mList, coords} = {"VertexTypes", "FormalCharges", "MassNumbers", "VertexCoordinates"} /. data;

                 If[NonPositive[nAtoms = Length[atoms]] || ! MatchQ[atoms, {__String}], Return[$Failed, Module]];

                 If[qList === "FormalCharges", qList = Table[0, {nAtoms}]];
                 If[mList === "MassNumbers", mList = Table[None, {nAtoms}]];
                 If[coords === "VertexCoordinates", coords = Table[{None}, {nAtoms}]];

                 Which[$VersionNumber < 12. && MatrixQ[coords, NumberQ] && Last[Dimensions[coords]] == 3,
                           coords /= 100.,
                           MatrixQ[coords, QuantityQ] && Last[Dimensions[coords]] == 3, 
                           coords = Normal[coords/Quantity["Angstroms"]]];

                 If[! (nAtoms == Length[qList] == Length[mList] == Length[coords] && MatchQ[qList, {__Integer}] &&
                        MatchQ[mList, {(_Integer | None) ..}] && 1 <= Last[Dimensions[coords]] <= 3 &&
                        (MatchQ[coords, {{None} ..}] || MatrixQ[coords, NumberQ])),
                    Return[$Failed, Module]];
  
                 {bonds, bondTypes} = {"EdgeRules", "EdgeTypes"} /. data;

                 If[(nBonds = Length[bonds]) != Length[bondTypes], Return[$Failed, Module]];
                 If[nBonds > 0,
                     If[! (MatchQ[bonds, {(_Integer -> _Integer) ..}] && (1 <= Min[#] <= Max[#] <= nAtoms) &[List @@@ bonds] &&
                            MatchQ[bondTypes, {("Single" | "Double" | "Triple" | "Aromatic") ..}]),
                         Return[$Failed, Module]]];

                 coords = Map[If[NumberQ[#], ToString[Round[#, $coordTol]], #] &, coords, {2}];
                 qList = Replace[ToString /@ qList, "0" -> None, {1}];
                 mList = If[# =!= None, ToString[#], #] & /@ mList;
                 atomArray = MapIndexed[XMLElement["atom", DeleteCases[Join[{"id" -> StringJoin["v", IntegerString @@ #2]},
                                                                                                                 Thread[{"elementType", "formalCharge", "isotope"} -> Take[#1, 3]],
                                                                                                                 Thread[Switch[Length[Last[#1]], 1, {"x"}, 2, {"x2", "y2"}, 3, {"x3", "y3", "z3"}] -> Last[#1]]],
                                                                                                          HoldPattern[_ -> None]], {}] &,
                                                       Transpose[{atoms, qList, mList, coords}]];
                 atomArray = XMLElement["atomArray", {}, atomArray];

                 If[nBonds > 0,
                     bonds = StringJoin[Insert[#, " ", 2]] & /@ Map[StringJoin["v", IntegerString[#]] &, List @@@ bonds, {2}];
                     bondTypes = bondTypes /. Thread[{"Single", "Double", "Triple", "Aromatic"} -> {"1", "2", "3", "A"}];
                     bondArray = MapIndexed[XMLElement["bond", Join[{"id" -> StringJoin["e", IntegerString @@ #2]},
                                                                                                 Thread[{"atomRefs2", "order"} -> #1]], {}] &,
                                                           Transpose[{bonds, bondTypes}]];
                     bondArray = XMLElement["bondArray", {}, bondArray],
                     bondArray = Unevaluated[]];

                 sch = {{"http://www.w3.org/2000/xmlns/", "xmlns"} -> "http://www.xml-cml.org/schema",
                            {"http://www.w3.org/2000/xmlns/", "cmlDict"} -> "http://www.xml-cml.org/dictionary/cml/",
                            {"http://www.w3.org/2000/xmlns/", "xsd"} -> "http://www.w3c.org/2001/XMLSchema",
                            {"http://www.w3.org/2000/xmlns/", "convention"} -> "http://www.xml-cml.org/convention", 
                            "convention" -> "convention:molecular"};

                 xmlObj = XMLObject["Document"][{XMLObject["Declaration"]["Version" -> "1.0"]},
                                                                     XMLElement["cml", sch,
                                                                                        {XMLElement["molecule",
                                                                                                             {"title" -> "Made with Mathematica " <> ToString[$VersionNumber], "id" -> "mol1"},
                                                                                                             {atomArray, bondArray}]}],
                                                                     {}];

                 Export[fn, xmlObj, "XML", opts, "AttributeQuoting" -> "\"", "ElementFormatting" -> cmlIndent]]

If[$VersionNumber >= 12.,

cmlExport[fn_String, data_Molecule?MoleculeQ, opts___] :=  
     Module[{atoms, atomArray, bonds, bondArray, bondTypes, coords2d, coords3d,
                  dbls, mList, nAtoms, nBonds, qList, sch, stereos, tets, xmlObj},

                 {atoms, qList, mList, coords2d, coords3d} =
                 MoleculeValue[data, {"VertexTypes", "FormalCharges", "MassNumbers", "AtomDiagramCoordinates", "AtomCoordinates"}];
                 coords3d = Normal[coords3d/Quantity["Angstroms"]];

                 {bonds, bondTypes} = MoleculeValue[data, {"EdgeRules", "EdgeTypes"}];
  
                 coords2d = Map[If[NumberQ[#], ToString[Round[#, $coordTol]], #] &, coords2d, {2}];
                 coords3d = Map[If[NumberQ[#], ToString[Round[#, $coordTol]], #] &, coords3d, {2}];
                 qList = Replace[ToString /@ qList, "0" -> None, {1}];
                 mList = If[# =!= None, ToString[#], #] & /@ mList;

                 atomArray = MapIndexed[XMLElement["atom", DeleteCases[Join[{"id" -> StringJoin["v", IntegerString @@ #2]},
                                                                                                                 Thread[{"elementType", "formalCharge", "isotope"} -> Take[#1, 3]],
                                                                                                                 Thread[{"x2", "y2"} -> #1[[4]]],
                                                                                                                 Thread[{"x3", "y3", "z3"} -> #1[[5]]]],
                                                                                                          HoldPattern[_ -> None]], {}] &, 
                                                       Transpose[{atoms, qList, mList, coords2d, coords3d}]];
                 atomArray = XMLElement["atomArray", {}, atomArray];

                 If[(nBonds = Length[bonds]) > 0,
                     bonds = StringJoin[Insert[#, " ", 2]] & /@ Map[StringJoin["v", IntegerString[#]] &, List @@@ bonds, {2}];
                     bondTypes = bondTypes /. Thread[{"Single", "Double", "Triple", "Aromatic"} -> {"1", "2", "3", "A"}];
                     bondArray = MapIndexed[XMLElement["bond", Join[{"id" -> StringJoin["e", IntegerString @@ #2]},
                                                                                                 Thread[{"atomRefs2", "order"} -> #1]], {}] &,
                                                           Transpose[{bonds, bondTypes}]];
                     bondArray = XMLElement["bondArray", {}, bondArray],
                     bondArray = Nothing];

                 stereos = StereochemistryElements /. Options[data, StereochemistryElements] /. StereochemistryElements -> {};
                 If[Length[stereos] > 0,
                     tets = Select[stereos, StringMatchQ[Lookup[#, "StereoType"], "Tetrahedral"] &];
                     If[Length[tets] > 0,
                         tets = With[{t = Lookup[#, "ChiralCenter"], l = Lookup[#, "Ligands"]},
                                           {StringJoin["v", IntegerString[t]], 
                                             XMLElement["atomParity", {"atomRefs4" -> StringRiffle[StringJoin["v", IntegerString[#]] & /@ PadLeft[Prepend[l, Lookup[#, "FiducialAtom"]], 4, t], " "]},
                                                                {Lookup[#, "Direction"] /. Thread[{"Clockwise", "Counterclockwise"} -> If[OddQ[Length[l]], Identity, Reverse][{"1", "-1"}]]}]}] & /@ tets;
                         atomArray = atomArray /. ((XMLElement["atom", dat_, att_] /; StringMatchQ["id" /. dat, #1] :> XMLElement["atom", dat, Append[att, #2]]) & @@@ tets)];

                     dbls = Select[stereos, StringMatchQ[Lookup[#, "StereoType"], "DoubleBond"] &];
                     If[Length[dbls] > 0,
                         dbls = With[{d = Lookup[#, "StereoBond"]},
                                            {StringRiffle[StringJoin["v", IntegerString[#]] & /@ d, " "], 
                                              XMLElement["bondStereo", {"atomRefs4" -> StringRiffle[StringJoin["v", IntegerString[#]] & /@ (Join[d, Lookup[#, "Ligands"]][[{3, 1, 2, 4}]]), " "]},
                                                                 {Lookup[#, "Value"] /. Thread[{"Together", "Opposite"} -> {"C", "T"}]}]}] & /@ dbls;
                         bondArray = bondArray /. ((XMLElement["bond", dat_, att_] /; StringMatchQ["atomRefs2" /. dat, #1] :> XMLElement["bond", dat, Append[att, #2]]) & @@@ dbls)]];

                 sch = {{"http://www.w3.org/2000/xmlns/", "xmlns"} -> "http://www.xml-cml.org/schema",
                            {"http://www.w3.org/2000/xmlns/", "cmlDict"} -> "http://www.xml-cml.org/dictionary/cml/",
                            {"http://www.w3.org/2000/xmlns/", "xsd"} -> "http://www.w3c.org/2001/XMLSchema",
                            {"http://www.w3.org/2000/xmlns/", "convention"} -> "http://www.xml-cml.org/convention", 
                            "convention" -> "convention:molecular"};

                 xmlObj = XMLObject["Document"][{XMLObject["Declaration"]["Version" -> "1.0"]}, 
                                                                     XMLElement["cml", sch,
                                                                                        {XMLElement["molecule",
                                                                                                             {"title" -> "Made with Mathematica " <> ToString[$VersionNumber], "id" -> "mol1"},
                                                                                                             {atomArray, bondArray}]}],
                                                                     {}];

                 Export[fn, xmlObj, "XML", opts, "AttributeQuoting" -> "\"", "ElementFormatting" -> cmlIndent]]
]

cmlExport[fn_String, _, opts___] := $Failed

SetAttributes[{cmlImport3D, cmlImport, cmlToDefault, cmlExport, cmlIndent}, ReadProtected];

Protect[cmlImport3D, cmlImport, cmlToDefault, cmlExport, cmlIndent];

End[]