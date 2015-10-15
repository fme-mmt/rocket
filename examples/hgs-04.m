%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 04: H2O equilibrium dissociation for different values of
%             temperature (T)
%
% H20 <-> H2 + O2 + H + O + OH
 
clear; clc;
 
format compact
 
 
p=1;                        % bar
T=2700 % K
 
% loop to compute composition using hgseq
 
 
comp=hgseq({'H2','O2','H2O','H','O','OH'},[2;1;0;0;0;0],T,1);
 
% nombre de mol, no concentracio
 
nH2=comp(1)
nO2=comp(2)
nH2O=comp(3)
nH=comp(4)
nO=comp(5)
nOH=comp(6)
 
 
 
%{
 
T =
        2700
nH2 =
    0.1596
nO2 =
    0.0560
nH2O =
    1.7721
nH =
    0.0332
nO =
    0.0125
nOH =
    0.1033
 
%}
