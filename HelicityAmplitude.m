(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: Helicity Amplitude Calculations              					*)

(* :Summary:	This is a Mathematica package for symbolic evaluation
				of Feynman diagrams and algebraic calculations in quantum
				field theory and elementary particle physics.
				Specifying in helicity amplitudes.           				*)

(* ------------------------------------------------------------------------ *)


BeginPackage[ "HelicityAmplitude`"];
Needs["FeynCalc`"]


(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)
(* :Section: Basic Four-Vector Operations              					*)
(* ------------------------------------------------------------------------ *)


Len::usage = "Len[a_,b_] computes the norm of the four-vectors given";
DotProduct::usage = "DotProduct[a_,b_] computes the product of the two four-vectors given";
PT::usage = "PT[a_,b_] computes the norm of the transverse momenta";


Begin[ "Private`"];


Len[a_,b_]:=Len[a,b]=Sqrt[a[[2]]*b[[2]]+a[[3]]*b[[3]]+a[[4]]*b[[4]]];
DotProduct[a_,b_]:=DotProduct[a,b]=a[[1]]*b[[1]]-(Len[a,b])^2;
PT[a_,b_]:=PT[a,b]=Sqrt[a[[2]]*b[[2]]+a[[3]]*b[[3]]];


End[];


ContractLV[a_,b_,c_,d_]:=ContractLV[a,b,c,d]=Contract[LeviCivita[\[Mu],\[Nu],\[Alpha],\[Beta]] . a . b . c . d];
LVexpr[a_,b_,c_,d_]:=LVexpr[a,b,c,d]=LeviCivitaTensor[4] . a . b . c . d;


(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)
(* :Section: FeynCalc Spin and Polarisation            					*)
(* ------------------------------------------------------------------------ *)


PV::usage = "PV[k,\[Mu]] is the abbriviation of PolarizationVector[k,\[Mu]]";
cPV::usage = "cPV[k,\[Mu]] is the abbriviation of ComplexConjugate[PolarizationVector[k,\[Mu]]]";
ContractLV::usage = "Symbolically Contract Levi-Civita with the four-vectors given";
LVexpr::usage = "Contract four-vectors (with given E,x,y,z values) with levi-civita";
PolVec::usage = "PolVecDef[k,\[Lambda]] computes the polarisation vector for spin-1 particle with four-momentum k and helicity \[Lambda]";
SpinVecML::usage = "SpinVecML[p,\[Lambda]]computes the spin vector for massless spin-1/2 particle with four-momentum k and helicity \[Lambda]";
SpinVecMS::usage = "SpinVecMS[p,\[Lambda]]computes the spin vector for massive spin-1/2 particle with four-momentum k and helicity \[Lambda]";


Begin[ "Private`"];


PV[k_,\[Mu]_]:=PV[k,\[Mu]]=PolarizationVector[k,\[Mu]];
cPV[k_,\[Mu]_]:=cPV[k,\[Mu]]=ComplexConjugate[PolarizationVector[k,\[Mu]]];
PolVec[k_,\[Lambda]_]:=PolVec[k,\[Lambda]]=1/Sqrt[2]*(-\[Lambda]*{0,k[[2]]*k[[4]],k[[3]]*k[[4]],-PT[k,k]*PT[k,k]}/(PT[k,k]*Len[k,k])-I*{0,-k[[3]],k[[2]],0}/PT[k,k]);
SpinVecML[p_,\[Lambda]_]:=SpinVecML[p,\[Lambda]]=\[Lambda]*{1,p[[2]]/Len[p,p],p[[3]]/Len[p,p],p[[4]]/Len[p,p]};
SpinVecMS[p_,\[Lambda]_]:=SpinVecMS[p,\[Lambda]]=(\[Lambda]*p[[1]])/(Sqrt[DotProduct[p,p]]*Len[p,p])*{Len[p,p]^2/p[[1]],p[[2]],p[[3]],p[[4]]};


End[];


EndPackage[];
