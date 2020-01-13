%This function calculates the
% [N2] Nonlinear stiffness matrix
% for the Beam-panel-futtler problem

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

function [N2]=CalcNonLinearBeam(A,TbInv,LengthX,LengthY,Wb)

%Select the number of Gauss-integration
% points you want
%5 points are more than enough for 
% a 2-node beam
NGc=8; %Number of Gauss integration points
%Get the Gauss points and weights from the 
% function below:
GaussConstants=GetGC(NGc); 

%Initialize the matrix
N2  =zeros(4,4);
%Start the numerical integration procedure
for Xi=1:NGc
  %Evaluating the physical X on the element
  X=LengthX*(GaussConstants(2,Xi)+1)/2;
  %Evaluating the needed vectors for the 
  % calculation of the lement matrices
  ct=CalcHwx(X);
  Bt=ct*TbInv; %1*4
  Theta=Bt*Wb; %1*1
  N2s = Bt'*Theta'*A*Theta*Bt; %4*4
  %performing the weighted summation
  N2  =N2  +GaussConstants(1,Xi)*N2s;
  %End of Calculation loop body
end
%Multiplying by Jacobian
N2  =N2*1.5*LengthX/2;