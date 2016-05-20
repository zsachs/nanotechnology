% Zachariah Sachs
% CHEM 26701
% 19 May 2013

% Nanotechnology

% Included in this script:
% -Constants and variable names.
% -1. Plotted absorption spectra; printed 'ausp#.png'
% -3. Simulated Au nanoparticle spectra; printed 'simauab.png'

%**************************************************************************
% Au constants and equations.
c=2.9979*10^8; %speed of lightt in m/s

%**************************************************************************
% Import data

abs1=importdata('auab0.txt');

% First column is "Wavelength (nm)"
% Second column is "Absorbance (??)"
La1=abs1(:,1);
ab1=abs1(:,2);

%**************************************************************************
% 1. Plot absorbance spectra
sp1=figure;
plot(La1,ab1);
hold on
xlabel('Wavelength (nm)');
ylabel('Absorbance (??)');
title('Au nanoparticle absorbance 1');
%print(sp1,'-dpng','ausp1');

%**************************************************************************
% 2. Calculate the average number of number of Au atoms in nanoparticle and 
% the concentration of Au nanoparticles in solution.	

rho=19.3; %density in g/cm^3
R=6*10^(-7); %average particle radius in cm
M=196.966569; %atomic mass in g/mol
NA=6.02214129*10^-23; %Avogadro's constant in 1/mol
N=(pi*rho*(2R)^3*NA)/(6*M); %number of atoms per nanoparticle

NT=2.0*10^(-5)*NA; %number of Au atoms added as HAuCl4
V=22*10^(-3); %volume of reaction solution in L
Cm=NT/(N*V*NA); %molar concentration of nanoparticles in mol/L

%**************************************************************************
% 3. Simulate the absorption spectrum of Au nanoparticles in solution.	

Vm=1/(1000*Cm); %molar volume of nanoparticles

epsilonm=1.7588; %unitless dielectric constant for water
epsiloninf=8.95; %unitless high freq dielectric constant for gold
omegaplas=1.30*10^16; %bulk plasmon frequency ih Hz
omegadump=1.95*10^14 %dumping frequency in Hz
omega=(2*pi*c)/L; %angular frequency in Hz
eepsilon=epsiloninf-(omegaplas^2/(omega^2+omegadump^2)); %stuff
eeepsilon=(omegaplas^2*omegadump)/(omega(omega^2+omegadump^2)); %stuffs

Cross=(24*pi*R^3*epsilonm^(3/2)*eeepsilon)/(L*((eepsilon+2epsilonm)^2+eeepsilon^2)); %absorption cross section

Q=Cross/(pi*R^2); %extinction cross section of a spherical particle

E=(3*10^(-3)*Vm*Q)/(4*2.303*R); %extinction coefficient in 1/(M*cm)

for L=0:2:1200;
    sim=E;
end

simul=figure;
plot(L,sim);
title('Simulated absorption spectra of Au nanoparticles');
xlabel('Wavelength (nm)');
ylabel('Absorption');
print(simiul,'-dpng','simauab');

%**************************************************************************
% 4. Based on the fact that the citrate anions cover the surface of each
% Au nanoparticle, explain what keeps the nanoparticles from aggregating.

% The citrate ions form a hydrophilic shell around the neutral 
% nanoparticles. 

%**************************************************************************
% 5. Explain the effect observed when the salt solution is added to the
% solution of Au nanoparticles.

% Addition of sodium chloride drastically increases the cation 
% concentration. Aqueous chloride ions are much smaller and more stable 
% than giant citrate ions, so the sodium binds to it and draws the 
% citrate off of the nanoparticles. I imagine gold, a heavy metal, has a
% poor index of refraction on it's own, so it transmits nothing when it
% aggregates.
