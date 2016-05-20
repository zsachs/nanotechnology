% Zachariah Sachs
% CHEM 26701
% 19 May 2013

% Nanotechnology

% Included in this script:
% -Constants and variable names.
% -1. Plotted spectra (absorption and then emission); printed
% 'cdseasp#.png' and 'cdseesp#.png'.
% -2. CdSe peak wavelength and FWHM; visualized.
% -2. Calculated nanocrystal diameter from Eqn (22) and abs/ems pekas
%  (Average ???).
% -3. Color tabulation from video.
% -4. Calculate $\sigma$ of the probable Gaussian distribution from FWHM in
% question #2.
% -5. Plot absorption vs. emission and fit with something probably linear
% -6. Explain the (linear) fit and the shift it represents.

%**************************************************************************
% CdSe constants
Eg=1.751*1.602*10^-19; %(eV) J gap energy
e=1.602*10^-19; %C electron charge
C=2.9979*10^8; %m/s speed of light
h=6.626*10^-34; %J*s lanck's constant
epcdseilon=10.6; %unitless dielectric constant for bulk CdSe
epsilon0=8.854*10^-12; %C^2*J^-1*m^-1 permittivity of vacuum
m0=9.109534*10^-31; %kg electron rest mass
me=0.13*m0; electron effective mass
mh=0.45*m0; hole effective mass

b=-(1.8*e^2)/(4*pi*epcdseilon*epsilon0); % from solving Equation (22)
c=(h^2/8)*((1/me)+(1/mh));

%**************************************************************************
% Import data

abs1=importdata('abs1.txt');

% First column is "Wavelength (nm)"
% Second column is "Absorbance (??)"
La1=abs1(:,1);
ab1=abs1(:,2);

emi1=importdata('ems1.txt');

% First column is "Wavelength (nm)"
% Second column is "Emission (??)"
Le1=emi1(:,1);
em1=emi1(:,2);
%**************************************************************************
% 1. Plot absorbance spectra
s1a=figure;
plot(La1,ab1);
hold on
xlabel('Wavelength (nm)');
ylabel('Absorbance (??)');
title('CdSe nanoparticle absorbance 1');
%print(s1a,'-dpng','cdseasp1');

%**************************************************************************
% 2. Tabulate absorbance peak wavelengths, estimate FWHM, and compute
% average nanocrystal diameter (from peak wavelength)
pkla1=La1(find(ab1==max(ab1),1)); %nm
plot(pkla1,max(ab1),'or');

hma1=find(ab1>max(ab1)/2,1,'last');
plot(La1(hma1),ab1(hma1),'or');

fwhma1=2*(La1(hma1)-pkla1); %assuming symmetry of peak? in nm

a1=Eg-((h*C)/(pkla1*10^-9)); %energy of peak from wavelength in J

xa1 = zeros(2,1);
da1 = sqrt(b^2 - 4*a1*c);
xa1(1) = ( -b + da1 ) / (2*a1);
xa1(2) = ( -b - da1 ) / (2*a1);  %solve Equation (22) quadratic in R

dma1=2*max(real(xa1)); % avg diameter in m.

% Do I take real part or modulus? Since the complex solutions will be
% conjugate the real part will be the same.

pkla=[pkla1]; in nm
dma=[dma1]; in m
fwhma=[fwhma1]; in nm

%**************************************************************************
% 1. Plot emission spectra
s1e=figure;
plot(Le1,em1);
hold on
xlabel('Wavelength (nm)');
ylabel('Emission (??)');
title('CdSe nanoparticle emission 1');
%print(s1e,'-dpng','cdseesp1');

%**************************************************************************
% 2. Tabulate emission peak wavelengths, estimate FWHM, and compute average
% nanocrystal diameter (from peak wavelength)
pkle1=Le1(find(em1==max(em1(900:end))),1); %nm
plot(pkle1,max(em1(900:end)),'or');

hme1=find(em1>max(em1(900:end))/2,1,'last');
plot(Le1(hme1),em1(hme1),'or');

fwhme1=2*(Le1(hme1)-pkle1); %assuming symmetry of peak? in nm

e1=Eg-((h*C)/(pkle1*10^-9)); %energy of peak from wavelength in J

xe1 = zeros(2,1);
de1 = sqrt(b^2 - 4*e1*c);
xe1(1) = ( -b + de1 ) / (2*e1);
xe1(2) = ( -b - de1 ) / (2*e1); %solve Equation (22) quadratic in R

dme1=2*max(real(xe1)); % avg diameter in m.

% Do I take real part or modulus? Ditto

pkle=[pkle1]; in nm
dme=[dme1]; in m
fwhme=[fwhme1]; in nm

%**************************************************************************
% 3. Tabulate nanocrystal diameters in correspondence with the absorption
% peak positions and a description of the color

% Note that emission peaks are not mentioned in this question.

% Also, for those who did the Nanotechnology lab there was a mistake in 
% question #3 for the semiconductor nano particles. It was pointed out that 
% the previous questions #2 & #3 were nearly the same. Rather, in #3 you 
% should plot some of the data from #2 and fit that. Sorry about the
% mistake.

% Basically, take the values from 2, and from the wavelength suggest a
% visible color. Check with video.

% See video

% The samples range from visibly pale yellow through orange to pale red.
% these are the wavelengths the transparent samples transmit, so the
% complimentary colors, namely violet through green, are absorbed.

%**************************************************************************
% 4. Calculate the distribution in the diameter of the nanocrystals for one
% sample based on the FWHM value and assuming a Gaussian distribution for
% the spectra band.

% The FWHM is $\Gamma=2.354*\sigma$. Solve this for $\sigma$ from the 
% numbers in #2

SDa=fwhma/2.354;
SDe=fwhme/2.354;

%**************************************************************************
% 5. Plot the relation between the absoption peak wavelength and emission
% peak wavelength

% For this I should just be able to make a scatter plot of abs vs emis. I
% might expect a linear trend or something corresponding to a shift.

abem=figure;
plot(pkle,pkla,'.');
title('CdSe nanoparticle Absorbption vs. Emission Wavelength');
xlabel('Emission wavelength (nm)');
ylabel('Absorption Wavelength (nm)');
print(abeam,'-dpng','above');

%**************************************************************************
% 6. Explain the shift.

% Consulting Figure 1 of the procedure, when energy is absorbed to create
% the electron-hole pair, the shift of Coulombic particles changes the 
% internuclear distance, making the energy gap between the excited and 
% ground state smaller, so the energy released by the relaxation will be 
% less, and thus of longer wavelength. One might think energy is lost 
% somewhere, but I think the relaxation of internuclear distance is just 
% slower than the electron-hole plugging because the nuclei are just so 
% massive and so respond slower to the Coulombic shift; the internuclear distance takes to long to reach equilibrium.

% Then again maybe not. Photon is absorbed. Electron -hole separation 
% absorbs some energy. Some of that energy is released when the nuclei
% respond to the charge separation. Plugging the hole releases what is 
% left. The nuclei then repel each other to get back to the start to 
% release just a tad more energy. Absorption is all an energy level shift.
% Emission is absorption minus the internuclear relaxations.
