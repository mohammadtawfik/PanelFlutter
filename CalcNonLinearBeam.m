function [N2]=CalcNonLinearBeam(A,TbInv,LengthX,LengthY,Wb)

NGc=8;
GaussConstants=GetGC(NGc);
N2  =zeros(4,4);
%Start the numerical integration procedure
for Xi=1:8
  X=LengthX*(GaussConstants(2,Xi)+1)/2;
  %**************************
  ct=CalcHwx(X);
  %**************************
  Bt=ct*TbInv; %1*4
  Theta=Bt*Wb; %1*1
  N2s = Bt'*Theta'*A*Theta*Bt; %4*4
  %performing the weighted summation
  N2  =N2  +GaussConstants(1,Xi)*N2s;
  %End of Calculation loop body
end
%Multiplying by Jacobian
N2  =N2*1.5*LengthX/2;
