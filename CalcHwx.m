%This function calculates the
% SLOPE of a third order polynomial
% at a point X

%Functions will work on Octave, FreeMat
% and Matlab
%Created by Mohammad Tawfik
%mohammad.tawfik@gmail.com 
%In assotiation with research paper
% published on ResearchGate.Net
%Author: Mohammad Tawfik
%Title: Panel Flutter
%DOI: 10.13140/RG.2.1.1537.6807
%Updated text link:
%https://www.researchgate.net/publication/275712979_Panel_Flutter
%More code abpout other topics in the text
% may be downloaded from:
% https://github.com/mohammadtawfik/PanelFlutter 

function Hw=CalcHwx(x)

Hw=zeros(1,4);
f1=0;
f2=1;
f3=2*x;
f4=3*x*x;
Hw(1)=f1;
Hw(2)=f2;
Hw(3)=f3;
Hw(4)=f4;
%NOTE:
% This may seem as over-killing the evaluation,
% IT IS ... But this is the way we are going to 
% develop the finite element models in all our
% future work, so just get familiar with it
% in order to understand the 2-D and higher
% order problems later