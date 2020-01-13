%This program calculates the
% post-buckling deflection of the panel
% with Aerodynamic lodaing

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

clear all
close all
clc

%Problem Data
E        =7.1e10; %Modulus of elasticity
Nu       =0.3;    %Poisson's ratio
Alpha    =22.5e-6;%Thermal expansion coefficient
Rho      =2700;   %Density
Thickness=0.002;  %Dimensions in meters
Nx       =10;     %Number of elements
Length   =1;
Width    =0.02;
%***************************
%Calculating basic quantities
%Cross-section area
Area=Width*Thickness;
%Second moment of area
Imoment=Width*Thickness*Thickness*Thickness/12;
%Element length
Le=Length/Nx; %Element Length
%Fundamental frequency
Wo=sqrt(E*Imoment/Rho/Area/(Length^4));
%***************************
%Calculating transformation Matrix
TbInv=inv([1,0 ,0    ,0; ...
           0,1 ,0    ,0; ...
           1,Le,Le*Le,Le*Le*Le; ...
           0,1 ,2*Le ,3*Le*Le]);
%***************************
%Calculating the linear matrices
% for the element
%***************************
%Stiffness Matrix
Kb=E*Imoment*[12  ,6*Le   ,-12  ,6*Le; ...
              6*Le,4*Le*Le,-6*Le,2*Le*Le; ...
              -12 ,-6*Le  ,12   ,-6*Le; ...
              6*Le,2*Le*Le,-6*Le,4*Le*Le]/Le/Le/Le;

M=Le*[156  ,22*Le    ,54    ,-13*Le; ...
      22*Le, 4*Le*Le ,13*Le ,-3*Le*Le; ...
      54   ,13*Le    ,156   ,-22*Le; ...
     -13*Le,-3*Le*Le ,-22*Le,4*Le*Le]/420;
%Mass Matrix
Mb =Rho*Area*M; 
%Aerodynamic Damping Matrix
Mg =E*Imoment/Wo/(Length^4)*M; 
%Therman stiffness matrix
Kt =Alpha*E*Area*[36  ,3*Le   ,-36  ,3*Le   ; ...
                        3*Le,4*Le*Le,-3*Le,-Le*Le ; ...
                        -36 ,-3*Le  ,36   ,-3*Le  ; ...
                        3*Le,-Le*Le ,-3*Le,4*Le*Le]/30/Le; 
%Aerodynamic Stiffness matrix
Aa =E*Imoment/(Length^3)*[-30  ,6*Le ,30   ,-6*Le; ...
                                -6*Le,0    ,6*Le ,-Le*Le; ...
                                -30  ,-6*Le,30   ,6*Le; ...
                                6*Le ,Le*Le,-6*Le,0]/60;
%***************************
%***************************
%Global stiffness and mass matrix assembly
%***************************
%Initializing an empty matrix
KbG=zeros(2*(Nx+1),2*(Nx+1));
KtG=zeros(2*(Nx+1),2*(Nx+1));
AaG=zeros(2*(Nx+1),2*(Nx+1));
MbG=zeros(2*(Nx+1),2*(Nx+1));
MgG=zeros(2*(Nx+1),2*(Nx+1));
%Assembling the global matrix
for ii=1:Nx
  KbG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
    KbG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Kb;
  KtG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
    KtG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Kt;
  AaG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
    AaG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Aa;
  MbG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
    MbG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Mb;
  MgG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))= ...
    MgG(2*ii-1:2*(ii+1),2*ii-1:2*(ii+1))+Mg;
end
%***************************
%The boundary conditions
%***************************
%For a Fixed-Fixed beam the first two
% and last two degree of freedom are fixed
BCs=[1,2,2*Nx+1,2*Nx+2];
%For simply supported beam, the first
% and before-last degees of freedom are fixed
%BCs=[1,2*Nx+1];

%Complementary boundary conditions
BCsC=1:1:2*(Nx+1);
BCsC(BCs)=[];

%Number of free nodes and degrees of freedom
ggg=size(BCs);
NFree=(2*(Nx+1)-ggg(2))/2; %Free nodes
gg=NFree*2;                %Free DOF
%A vector with the number of the 
% degrees of freedom of the displacements
% onlt; ie: without the slopes
BendingDisplDOF=1:2:2*(Nx+1);
%Applying the boundary conditions
KbG(BCs,:)=[]; KbG(:,BCs)=[];
KtG(BCs,:)=[]; KtG(:,BCs)=[];
AaG(BCs,:)=[]; AaG(:,BCs)=[];
MbG(BCs,:)=[]; MbG(:,BCs)=[];
MgG(BCs,:)=[]; MgG(:,BCs)=[];

%Finding the buckling temperature
TInitial=min(eig(inv(KtG)*(KbG)));
%Temperature increase in the temperature loop
TIncrem=0.005;

%Initializing the displacements vector
Displacements=zeros(2*(Nx+1),1); %all displ
WInitial=zeros(gg,1); %displ. withoud boundaries
%To start the nonlinear calculations
% the displacements have to be assumed 
% having non-zero values
% otherwise, the nonlinear terms will be zero!
for ii=1:gg
  WInitial(ii)=Thickness*0.5;
end
%Storing the assumed displacements in the 
% total displacements vector
Displacements(BCsC)=WInitial;
%First element in the temperature vector
TempVect(1)=TInitial/TInitial;

%This vector is contains the values
% of the dynamic pressure used in the 
% calculations - NOTE: you may use ANY values
LamX=[0,40,60,100,120];
%Looping on the values of the dynamic
% pressure. The maximum number needs 
% NOT to exceed the size of the vector above
for LLL=1:1:2;
  Lamda=LamX(LLL) %Getting the value from the vector
  %looping for the temperature
  %We already kknow what is the solution
  % at the first step - essentially zero
  for ii=2:1:40+1
    %Evaluating the temperature
    Temp=TInitial+TIncrem*(ii-1);
    %Storing the temperature scaled to the 
    % basic buckling temperature
    TempVect(ii)=Temp/TInitial;
    %If the solution obtained in the pervious
    % temperature is zero, the nonlinear
    % model will become AND STAY linear!
    % that is why we have to, artificially, 
    % force the values NOT to be zero
    if max(Displacements(BendingDisplDOF))<Thickness/100
      WInitial=zeros(gg,1);
      for jj=1:gg
        WInitial(jj)=Thickness*0.5;
      end
      Displacements(BCsC)=WInitial;
    end
    %calculating the linear part of the stiffness and load
    KL=KbG+Lamda*AaG-Temp*KtG;
    Iteration=0;
    %Initializing the Newton-Raphson
    % incrimental changes (it should
    % reach zero as iterations are performed)
    DelW=WInitial/5;
    %Newton-Raphson iterations
    while max(DelW)>=Thickness/1000000
      Iteration=Iteration+1;
      %Calculating the Nonlinear
      % Stiffness of each element
      
      %Initializing the global
      % nonlinear matrix
      N2G=zeros(2*(Nx+1),2*(Nx+1));
      
      %looping for the elements
      for jj=1:Nx
        %Extracting Node Displacements
        Wbe=Displacements(2*jj-1:2*jj+2);
        %Evaluating the nonlinear matrix
        % using numerical integration
        N2e=CalcNonLinearBeam(E*Area,TbInv,Le,Width,Wbe);
        %Assembling Non-Linear Stiffness Matrices
        N2G(2*jj-1:2*(jj+1),2*jj-1:2*(jj+1))= ...
          N2G(2*jj-1:2*(jj+1),2*jj-1:2*(jj+1))+N2e;
        %END of Non-Linear Assembly
      end %for elements
      %Applying boundary conditions
      N2G(BCs,:)=[]; N2G(:,BCs)=[];
      %Performing the Newton-Raphson procedure
      KTan=KL+N2G;
      KTot=KL+N2G/3;
      PsiI=KTot*WInitial;
      DelW=-inv(KTan)*PsiI;
      WInitial=WInitial+DelW;
      %In-case the calculations diverge
      % the check on the number of iterations
      % will prevent infinite loops
      if Iteration>20
        break
      end
      %Updating Nodes' displacements
      Displacements(BCsC)=WInitial;
    end %while Newton-ruphson Iteration
    %Storing the maximum displacement
    % calculated at this temperature
    wmax(LLL,ii)=max(Displacements(BendingDisplDOF))/Thickness;
  end %for temperature loop
  %Storing the results for exporting to
  % text file at the end
  XYZ(:,LLL)=wmax(LLL,:)';
end %for lamda loop
figure(11)
plot(TempVect,wmax)
grid
xlabel("Temp/TCritical")
ylabel("Wmax/Thickness")

%Adding the temperature to the data stored
XYZ=[TempVect', XYZ];
%Exporting the data to a file
save XYZ.txt XYZ -ascii