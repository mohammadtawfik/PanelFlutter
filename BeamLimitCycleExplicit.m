%This code evaluates and plots the stability
% boundaries for a panel using a beam model
%The equations use the linear tems only

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

%Clearing the memory and display
clear all
clc
close all
%***************************
%Problem Data
NE=10; %number of elements
Length=.50; %beam length
Width=0.02; %beam width
Thickness=0.001; %beam thickness
Modulus=71e9; %Modulus of Elasticity Aluminum (GPa)
Rho=2700; %Density (Kg/m^3)
Alpha=22.5e-6; %Thermal Expansion coefficientt
%***************************
%Calculating basic quantities
%Cross-section area
Area=Width*Thickness;
%Second moment of area
Imoment=Width*Thickness*Thickness*Thickness/12;
Le=Length/NE; %Element Length
%Fundamental frequency
Wo=sqrt(Modulus*Imoment/Rho/Area/(Length^4));
%***************************
%Calculating transformation Matrix
TbInv=inv([1,0 ,0    ,0; ...
           0,1 ,0    ,0; ...
           1,Le,Le*Le,Le*Le*Le; ...
           0,1 ,2*Le ,3*Le*Le]);%***************************
%Element matrices
%Stiffness Matrix
Ke=Modulus*Imoment*[12  ,6*Le   ,-12  ,6*Le; ...
6*Le,4*Le*Le,-6*Le,2*Le*Le; ...
-12 ,-6*Le  ,12   ,-6*Le; ...
6*Le,2*Le*Le,-6*Le,4*Le*Le]/Le/Le/Le;

M=Le*[156  ,22*Le    ,54    ,-13*Le; ...
22*Le, 4*Le*Le ,13*Le ,-3*Le*Le; ...
54   ,13*Le    ,156   ,-22*Le; ...
-13*Le,-3*Le*Le ,-22*Le,4*Le*Le]/420;
%Mass Matrix
Me =Rho*Area*M; 
%Aerodynamic Damping Matrix
Mg =Modulus*Imoment/Wo/(Length^4)*M; 
%Geometric Matrix
Kg =Alpha*Modulus*Area*[36  ,3*Le   ,-36  ,3*Le   ; ...
3*Le,4*Le*Le,-3*Le,-Le*Le ; ...
-36 ,-3*Le  ,36   ,-3*Le  ; ...
3*Le,-Le*Le ,-3*Le,4*Le*Le]/30/Le; 
%Aerodynamic Stiffness matrix
Ka =Modulus*Imoment/(Length^3)*[-30  ,6*Le ,30   ,-6*Le; ...
-6*Le,0    ,6*Le ,-Le*Le; ...
-30  ,-6*Le,30   ,6*Le; ...
6*Le ,Le*Le,-6*Le,0]/60;
%***************************
%Global stiffness and mass matrix assembly
%***************************
%Initializing an empty matrix
KeGlobal=zeros(2*(NE+1),2*(NE+1));
KgGlobal=zeros(2*(NE+1),2*(NE+1));
KaGlobal=zeros(2*(NE+1),2*(NE+1));
MeGlobal=zeros(2*(NE+1),2*(NE+1));
MgGlobal=zeros(2*(NE+1),2*(NE+1));
%Assembling the global matrix
for ii=1:NE
  KeGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
  KeGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Ke;
  KgGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
  KgGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Kg;
  KaGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
  KaGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Ka;
  MeGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
  MeGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Me;
  MgGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
  MgGlobal(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Mg;
end
%***************************
%The boundary conditions
%***************************
%For a Fixed-Fixed beam the first two and last two
% degree of freedom are fixed
BCs=[1,2,2*NE+1,2*NE+2];
%For simply supported beam, the first and before-last
% degees of freedom are fixed
%BCs=[1,2*NE+1];
BCsC=1:1:2*(NE+1);
BCsC(BCs)=[];
Displacements=zeros(2*(NE+1),1);

%Applying the boundary conditions
KeGlobal(BCs,:)=[];
KeGlobal(:,BCs)=[];
KgGlobal(BCs,:)=[];
KgGlobal(:,BCs)=[];
KaGlobal(BCs,:)=[];
KaGlobal(:,BCs)=[];
MeGlobal(BCs,:)=[];
MeGlobal(:,BCs)=[];
MgGlobal(BCs,:)=[];
MgGlobal(:,BCs)=[];

%***************************
%Calculating the buckling temperature
% for the static case - without air
TBucklingStatic=min(eig(inv(KgGlobal)*KeGlobal))
%***************************

%Looping for different values of the dynamic pressure
% to evaluate the buckling temperature at each speed
for ii=1:101
  LambdaBuc(ii)=(ii-1)*2.0;
  TTBuc(ii)=min(eig(inv(KgGlobal)*(KeGlobal+LambdaBuc(ii)*KaGlobal)));
endfor
TTBuc=TTBuc/TBucklingStatic;

%Looping for different temperatures to evaluate
% the critical dynamic pressure
for ii=1:11
  TempFlut=(ii-1)*0.15*TBucklingStatic; 
  %The number 0.017 in the above line is
  % totally arbitrary for scaling the graph
  %**********************
  %Searching for the flutter pressure
  for jj=1:1001
    LL=1*(jj-1); %Value of Lamda to test for flutter
    KK=KeGlobal-TempFlut*KgGlobal+LL*KaGlobal;
    MM=MeGlobal+MgGlobal;
    [vv,Kappa]=eig(inv(MM)*KK);
    Kappa=diag(Kappa);
    IK=imag(Kappa);
    %If any of the imaginary eigenvalues
    % is NOT ZERO then record the critical
    % pressure
    if IK'*IK!=0
      TTFlut(ii)=TempFlut/TBucklingStatic;
      LambdaFlut(ii)=LL;
      [Kappa,SortIndex]=sort(Kappa);
      Phi=real(vv(:,SortIndex(1)));
      Phi=Phi/max(abs(Phi));
      break
    endif
  endfor
  LamdaCycle(1,ii)=LL;
  Displacements(BCsC)=Phi;
  Cvect(1)=0;
  for kk=1:5
    c0=0.5*kk;
    Cvect(kk+1)=c0;
    %Calculating the Nonlinear Stiffness of each element
    N2G=zeros(2*(NE+1),2*(NE+1));
    %looping for the elements
    for jj=1:NE
      %Extracting Node Displacements
      Wbe=c0*Thickness*Displacements(2*jj-1:2*jj+2);
      N2e=CalcNonLinearBeam(Modulus*Area,TbInv,Le,Width,Wbe);
      %Assembling Non-Linear Stiffness Matrices
      N2G(2*jj-1:2*(jj+1),2*jj-1:2*(jj+1))= ...
      N2G(2*jj-1:2*(jj+1),2*jj-1:2*(jj+1))+N2e;
      %END of Non-Linear Assembly
    end %for elements
    N2G(BCs,:)=[]; N2G(:,BCs)=[];
    for jj=LL:2*LL
      KK=KeGlobal-TempFlut*KgGlobal+jj*KaGlobal+N2G;
      MM=MeGlobal+MgGlobal;
      [vv,Kappa]=eig(inv(MM)*KK);
      Kappa=diag(Kappa);
      IK=imag(Kappa);
      %If any of the imaginary eigenvalues
      % is NOT ZERO then record the critical
      % pressure
      if IK'*IK!=0
        [kk ii jj]
        LamdaCycle(kk+1,ii)=jj;
        break
      endif
    endfor
  endfor
endfor
%Plorring the resulting graphs
figure(1)  
plot(real(TTBuc),LambdaBuc,TTFlut,LambdaFlut)
grid
ylabel("Lambda")
xlabel("Temp/TempCritical")

figure(2)
plot(LamdaCycle,Cvect)
grid
ylabel("Limit Cycle amplitude (co) / Thickness")
xlabel("Dynamic Pressure")

figure(3)
plot(TTFlut, LamdaCycle')
grid
ylabel("Dynamic Pressure")
xlabel("Temperature/TCritical")

save XYZ.txt LamdaCycle -ascii