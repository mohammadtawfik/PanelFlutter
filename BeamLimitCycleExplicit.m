%This code evaluates and plots the variation
% of the limit cycle amplitude
% for a panel using a beam model

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
%Thermal Stiffness Matrix
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

%Applying the boundary conditions
KeGlobal(BCs,:)=[]; KeGlobal(:,BCs)=[];
KgGlobal(BCs,:)=[]; KgGlobal(:,BCs)=[];
KaGlobal(BCs,:)=[]; KaGlobal(:,BCs)=[];
MeGlobal(BCs,:)=[]; MeGlobal(:,BCs)=[];
MgGlobal(BCs,:)=[]; MgGlobal(:,BCs)=[];

%Creating the complementary-boundary-conditions
% vector. This includes the numbers of the 
% free degrees of freedom
BCsC=1:1:2*(NE+1);
BCsC(BCs)=[];

%Initializing a vector with the latest
% values of ALL the degrees of freedom
Displacements=zeros(2*(NE+1),1);

%***************************
%Calculating the buckling temperature
% for the static case - without air
TBucklingStatic=min(eig(inv(KgGlobal)*KeGlobal))
%***************************


%Looping for different temperatures to
% evaluate the critical dynamic pressure
for ii=1:11
  %Evaluating the temperature
  TempFlut=(ii-1)*0.15*TBucklingStatic; 
  %Storing it in a vector
  % normalized with the basic buckliing temp.
  TTFlut(ii)=TempFlut/TBucklingStatic;
  %The number 0.015 in the above line is
  % totally arbitrary
  %**********************
  %Searching for the flutter pressure
  for jj=1:1001
    %Value of Lamda to test for flutter
    LL=1*(jj-1);
    %Evaluating the total mass
    % and stiffness matrices
    KK=KeGlobal-TempFlut*KgGlobal+LL*KaGlobal;
    MM=MeGlobal+MgGlobal;
    %Obtaining the eigenvalues and vectors
    [vv,Kappa]=eig(inv(MM)*KK);
    %Putting the eigenvalues in a vector
    Kappa=diag(Kappa);
    %Extracting the imaginary values 
    % of the eigenvalues
    IK=imag(Kappa);
    %If any of the imaginary eigenvalues
    % is NOT ZERO then record the critical
    % pressure
    if IK'*IK!=0
      %Sort the eigenvalues
      [Kappa,SortIndex]=sort(Kappa);
      %Obtain the mode-shape associated with
      % the flutter pressure
      Phi=real(vv(:,SortIndex(1)));
      %Normalize the mode shape with the 
      % maximum value
      Phi=Phi/max(abs(Phi));
      %Exit the dynami-pressure loop
      break
    endif
  endfor
  %Store the flutter dynamic pressure
  % as the first point on the graph
  LamdaCycle(1,ii)=LL;
  %Initialize the displacements with the 
  % mode-shape of flutter
  Displacements(BCsC)=Phi;
  %Initialize the limit cycle amplitude vector
  Cvect(1)=0;
  %Loop for different limit cycle amplituudes
  for kk=1:5
    c0=0.5*kk; %The amplitude
    Cvect(kk+1)=c0; %Store the amplitude
    %Calculating the Nonlinear Stiffness
    % matrix for the structure
    N2G=zeros(2*(NE+1),2*(NE+1));
    %looping for the elements
    for jj=1:NE
      %Extracting Node Displacements
      %NOTE: the amplitude is multiplied
      % by the thickness to get a physical
      % value for the displacements
      Wbe=c0*Thickness*Displacements(2*jj-1:2*jj+2);
      %Evaluating the element nonlinear matrix
      N2e=CalcNonLinearBeam(Modulus*Area,TbInv,Le,Width,Wbe);
      %Assembling Non-Linear Stiffness Matrices
      N2G(2*jj-1:2*(jj+1),2*jj-1:2*(jj+1))= ...
        N2G(2*jj-1:2*(jj+1),2*jj-1:2*(jj+1))+N2e;
      %END of Non-Linear Assembly
    end %for elements
    %Applying boundary conditions
    N2G(BCs,:)=[]; N2G(:,BCs)=[];
    %Searching for the dynamic pressure 
    % that will satisfy the equilibrium
    for jj=LL:2*LL
      %Evaluating the total stiffness
      % matric - including the nonlinear matrix
      KK=KeGlobal-TempFlut*KgGlobal+jj*KaGlobal+N2G;
      %Evaluating the total mass matrix
      MM=MeGlobal+MgGlobal;
      %Obtaining the eigenvalues
      Kappa=eig(inv(MM)*KK);
      IK=imag(Kappa);
      %If any of the imaginary eigenvalues
      % is NOT ZERO then record the critical
      % pressure
      if IK'*IK!=0
        %This diaplys the current satus of calculations
        [kk ii jj]
        %Storing the equilibrium
        % dynamic pressure
        LamdaCycle(kk+1,ii)=jj;
        %Exit the dynamic pressure loop
        break
      endif
    endfor
  endfor
endfor
%Ploting the resulting graphs
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
%Storing the results in a text file
save XYZ.txt LamdaCycle -ascii