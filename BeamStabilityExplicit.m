    %Function will work on Octave, FreeMat, and Matlab
    %Create by Mohammad Tawfik
    %mohammad.tawfik@gmail.com 
    %In assotiation with research paper published on 
    %ResearchGate.Net
    %Author: Mohammad Tawfik
    %Title: Panel Flutter
    %DOI: 10.13140/RG.2.1.1537.6807

    %This code evaluates and plots the stability
    % boundaries for a panel using a beam model
    %The equations use the linear tems only

%Clearing the memory and display
clear all
clc
close all
%Problem Data
NE=10; %number of elements
Length=.50; %beam length
Width=0.02; %beam width
Thickness=0.001; %beam thickness
Modulus=71e9; %Modulus of Elasticity Aluminum (GPa)
Rho=2700; %Density (Kg/m^3)
Alpha=22.5e-6; %Thermal Expansion coefficientt
%Cross-section area
Area=Width*Thickness;
%Second moment of area
Imoment=Width*Thickness*Thickness*Thickness/12;
Le=Length/NE; %Element Length
%Fundamental frequency
Wo=sqrt(Modulus*Imoment/Rho/Area/(Length^4));
%Element matrices
Ke=Modulus*Imoment*[12  ,6*Le   ,-12  ,6*Le; ...
                    6*Le,4*Le*Le,-6*Le,2*Le*Le; ...
                    -12 ,-6*Le  ,12   ,-6*Le; ...
                    6*Le,2*Le*Le,-6*Le,4*Le*Le]/Le/Le/Le;

M=Le*[156  ,22*Le    ,54    ,-13*Le; ...
      22*Le, 4*Le*Le ,13*Le ,-3*Le*Le; ...
      54   ,13*Le    ,156   ,-22*Le; ...
     -13*Le,-3*Le*Le ,-22*Le,4*Le*Le]/420;
Me =Rho*Area*M; %Mass Matrix
Mg =Modulus*Imoment/Wo/(Length^4)*M; %Aerodynamic Damping Matrix
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
%Global stiffness and mass matrix assembly
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
%For a Fixed-Fixed beam the first two and last two
% degree of freedom are fixed
BCs=[1,2,2*NE+1,2*NE+2];
%For simply supported beam, the first and before-last
% degees of freedom are fixed
%BCs=[1,2*NE+1];

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

%Calculating the buckling temperature
% for the static case - without air
TBucklingStatic=min(eig(inv(KgGlobal)*KeGlobal))

%Looping for different values of the dynamic pressure
% to evaluate the buckling temperature at each speed
for ii=1:101
  LambdaBuc(ii)=(ii-1)*2.0;
  TTBuc(ii)=min(eig(inv(KgGlobal)*(KeGlobal+LambdaBuc(ii)*KaGlobal)));
endfor
TTBuc=TTBuc/TBucklingStatic;

%Looping for different temperatures to evaluate
% the critical dynamic pressure
for ii=1:101
  TempFlut=(ii-1)*.017*TBucklingStatic; 
  %The number 0.017 in the above line is
  % totally arbitrary for scaling the graph
  for jj=1:1001
    LL=(jj-1);
    KK=KeGlobal-TempFlut*KgGlobal+LL*KaGlobal;
    MM=MeGlobal+MgGlobal;
    Kappa=eig(inv(MM)*KK);
    IK=imag(Kappa);
    %If any of the imaginary eigenvalues
    % is NOT ZERO then record the critical
    % pressure
    if IK'*IK!=0
      TTFlut(ii)=TempFlut/TBucklingStatic;
      LambdaFlut(ii)=LL;
      break
    endif
  endfor
endfor
figure(1)  
plot(real(TTBuc),LambdaBuc,TTFlut,LambdaFlut)
grid
ylabel("Lambda")
xlabel("Temp/TempCritical")
