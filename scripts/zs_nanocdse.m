% Zachariah Sachs
% CHEM 26701
% 19 May 2013

% Nanotechnology

% Included in this script:
% -Constants and variable names.
% -1. and 2. Plotted spectra (absorption and then emission) with CdSe peak
%  wavelength and FWHM marked ; printed 'cdseasp#.png' and 'cdseesp#.png'.
%  Peaks written 'peaks.csv'
% -2. Calculated nanocrystal diameter from Eqn (22) and abs(/ems) peaks
%  (Average ???); written 'diam.csv'.
% -3. Color tabulation from video. FWHM values written 'FWHM.csv'.
% -3. Plot nanocrystal diameter vs. absorption energy with 
%  particle-in-a-box fit; printed 'PiBx.png'
% Fix PIB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% -4. Calculate $\sigma$ of the probable Gaussian distribution from FWHM in
%  question #2. Written 'distrib.csv'.
% -5. Plot absorption vs. emission with linear fit; shift written
%  'shft.csv'
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
me=0.13*m0; %electron effective mass
mh=0.45*m0; %hole effective mass

b=-(1.8*e^2)/(4*pi*epcdseilon*epsilon0); % from solving Equation (22)
c=(h^2/8)*((1/me)+(1/mh));

%**************************************************************************
% Import absorbance data

abs1=importdata('abs1.txt');
abs2=importdata('abs2.txt');
abs3=importdata('abs3.txt');
abs4=importdata('abs4.txt');
abs5=importdata('abs5.txt');
abs6=importdata('abs6.txt');
abs7=importdata('abs7.txt');
abs8=importdata('abs8.txt');
abs9=importdata('abs9.txt');
abs10=importdata('abs10.txt');
abs11=importdata('abs11.txt');
abs12=importdata('abs12.txt');
abs13=importdata('abs13.txt');

% First column is "Wavelength (nm)"
% Second column is "Absorbance (??)"
La1=abs1(:,1);
ab1=abs1(:,2);
La2=abs2(:,1);
ab2=abs2(:,2);
La3=abs3(:,1);
ab3=abs3(:,2);
La4=abs4(:,1);
ab4=abs4(:,2);
La5=abs5(:,1);
ab5=abs5(:,2);
La6=abs6(:,1);
ab6=abs6(:,2);
La7=abs7(:,1);
ab7=abs7(:,2);
La8=abs8(:,1);
ab8=abs8(:,2);
La9=abs9(:,1);
ab9=abs9(:,2);
La10=abs10(:,1);
ab10=abs10(:,2);
La11=abs11(:,1);
ab11=abs11(:,2);
La12=abs12(:,1);
ab12=abs12(:,2);
La13=abs13(:,1);
ab13=abs13(:,2);

% Import emission data

emi1=importdata('ems1.txt');
emi2=importdata('ems2.txt');
emi3=importdata('ems3.txt');
emi4=importdata('ems4.txt');
emi5=importdata('ems5.txt');
emi6=importdata('ems6.txt');
emi7=importdata('ems7.txt');
emi8=importdata('ems8.txt');
emi9=importdata('ems9.txt');
emi10=importdata('ems10.txt');
emi11=importdata('ems11.txt');
emi12=importdata('ems12.txt');
emi13=importdata('ems13.txt');

% First column is "Wavelength (nm)"
% Second column is "Emission (??)"
Le1=emi1(:,1);
em1=emi1(:,2);
Le2=emi2(:,1);
em2=emi2(:,2);
Le3=emi3(:,1);
em3=emi3(:,2);
Le4=emi4(:,1);
em4=emi4(:,2);
Le5=emi5(:,1);
em5=emi5(:,2);
Le6=emi6(:,1);
em6=emi6(:,2);
Le7=emi7(:,1);
em7=emi7(:,2);
Le8=emi8(:,1);
em8=emi8(:,2);
Le9=emi9(:,1);
em9=emi9(:,2);
Le10=emi10(:,1);
em10=emi10(:,2);
Le11=emi11(:,1);
em11=emi11(:,2);
Le12=emi12(:,1);
em12=emi12(:,2);
Le13=emi13(:,1);
em13=emi13(:,2);

%**************************************************************************
% 2. Tabulate absorbance peak wavelengths, estimate FWHM, and compute
% average nanocrystal diameter (from peak wavelength)
pkla1=La1(find(ab1==max(ab1),1)); %nm
hma1=find(ab1>max(ab1)/2,1,'last');
fwhma1=2*(La1(hma1)-pkla1); %assuming symmetry of peak? in nm
a1=Eg-((h*C)/(pkla1*10^-9)); %energy of peak from wavelength in J
xa1 = zeros(2,1);
da1 = sqrt(b^2 - 4*a1*c);
xa1(1) = ( -b + da1 ) / (2*a1);
xa1(2) = ( -b - da1 ) / (2*a1);  %solve Equation (22) quadratic in R
dma1=2*max(real(xa1)); % avg diameter in m.

% Do I take real part or modulus? Since the complex solutions will be
% conjugate the real part will be the same.

pkla2=La2(find(ab2==max(ab2(600:end)),1,'last')); %nm
hma2=find(ab2>max(ab2(600:end))/2,1,'last');
fwhma2=2*(La2(hma2)-pkla2); %assuming symmetry of peak? in nm
a2=Eg-((h*C)/(pkla2*10^-9)); %energy of peak from wavelength in J
xa2 = zeros(2,1);
da2 = sqrt(b^2 - 4*a2*c);
xa2(1) = ( -b + da2 ) / (2*a2);
xa2(2) = ( -b - da2 ) / (2*a2);  %solve Equation (22) quadratic in R
dma2=2*max(real(xa2)); % avg diameter in m.

pkla3=La3(find(ab3==max(ab3(435:end)),1,'last')); %nm
hma3=find(ab3>max(ab3(435:end))/2,1,'last');
fwhma3=2*(La3(hma3)-pkla3); %assuming symmetry of peak? in nm
a3=Eg-((h*C)/(pkla3*10^-9)); %energy of peak from wavelength in J
xa3 = zeros(2,1);
da3 = sqrt(b^2 - 4*a3*c);
xa3(1) = ( -b + da3 ) / (2*a3);
xa3(2) = ( -b - da3 ) / (2*a3);  %solve Equation (22) quadratic in R
dma3=2*max(real(xa3)); % avg diameter in m.

pkla4=La4(find(ab4==max(ab4(485:end)),1,'last')); %nm
hma4=find(ab4>max(ab4(485:end))/2,1,'last');
fwhma4=2*(La4(hma4)-pkla4); %assuming symmetry of peak? in nm
a4=Eg-((h*C)/(pkla4*10^-9)); %energy of peak from wavelength in J
xa4 = zeros(2,1);
da4 = sqrt(b^2 - 4*a4*c);
xa4(1) = ( -b + da4 ) / (2*a4);
xa4(2) = ( -b - da4 ) / (2*a4);  %solve Equation (22) quadratic in R
dma4=2*max(real(xa4)); % avg diameter in m.

pkla5=La5(find(ab5==max(ab5(580:end)),1,'last')); %nm
hma5=find(ab5>max(ab5(580:end))/2,1,'last');
fwhma5=2*(La5(hma5)-pkla5); %assuming symmetry of peak? in nm
a5=Eg-((h*C)/(pkla5*10^-9)); %energy of peak from wavelength in J
xa5 = zeros(2,1);
da5 = sqrt(b^2 - 4*a1*c);
xa5(1) = ( -b + da5 ) / (2*a5);
xa5(2) = ( -b - da5 ) / (2*a5);  %solve Equation (22) quadratic in R
dma5=2*max(real(xa5)); % avg diameter in m.

pkla6=La6(find(ab6==max(ab6(600:end)),1,'last')); %nm
hma6=find(ab6>max(ab6(600:end))/2,1,'last');
fwhma6=2*(La6(hma6)-pkla6); %assuming symmetry of peak? in nm
a6=Eg-((h*C)/(pkla6*10^-9)); %energy of peak from wavelength in J
xa6 = zeros(2,1);
da6 = sqrt(b^2 - 4*a6*c);
xa6(1) = ( -b + da6 ) / (2*a6);
xa6(2) = ( -b - da6 ) / (2*a6);  %solve Equation (22) quadratic in R
dma6=2*max(real(xa6)); % avg diameter in m.

pkla7=La7(find(ab7==max(ab7(650:end)),1,'last')); %nm
hma7=find(ab7>max(ab7(650:end))/2,1,'last');
fwhma7=2*(La7(hma7)-pkla7); %assuming symmetry of peak? in nm
a7=Eg-((h*C)/(pkla7*10^-9)); %energy of peak from wavelength in J
xa7 = zeros(2,1);
da7 = sqrt(b^2 - 4*a7*c);
xa7(1) = ( -b + da7 ) / (2*a7);
xa7(2) = ( -b - da7 ) / (2*a7);  %solve Equation (22) quadratic in R
dma7=2*max(real(xa7)); % avg diameter in m.

pkla8=La8(find(ab8==max(ab8(770:end)),1,'last')); %nm
hma8=find(ab8>max(ab8(770:end))/2,1,'last');
fwhma8=2*(La8(hma8)-pkla8); %assuming symmetry of peak? in nm
a8=Eg-((h*C)/(pkla8*10^-9)); %energy of peak from wavelength in J
xa8 = zeros(2,1);
da8 = sqrt(b^2 - 4*a8*c);
xa8(1) = ( -b + da8 ) / (2*a8);
xa8(2) = ( -b - da8 ) / (2*a8);  %solve Equation (22) quadratic in R
dma8=2*max(real(xa8)); % avg diameter in m.

pkla9=La9(find(ab9==max(ab9(810:end)),1,'last')); %nm
hma9=find(ab9>max(ab9(810:end))/2,1,'last');
fwhma9=2*(La9(hma9)-pkla9); %assuming symmetry of peak? in nm
a9=Eg-((h*C)/(pkla9*10^-9)); %energy of peak from wavelength in J
xa9 = zeros(2,1);
da9 = sqrt(b^2 - 4*a9*c);
xa9(1) = ( -b + da9 ) / (2*a9);
xa9(2) = ( -b - da9 ) / (2*a9);  %solve Equation (22) quadratic in R
dma9=2*max(real(xa9)); % avg diameter in m.

pkla10=La10(find(ab10==max(ab10(750:end)),1,'last')); %nm
hma10=find(ab10>max(ab10(750:end))/2,1,'last');
fwhma10=2*(La10(hma10)-pkla10); %assuming symmetry of peak? in nm
a10=Eg-((h*C)/(pkla10*10^-9)); %energy of peak from wavelength in J
xa10 = zeros(2,1);
da10 = sqrt(b^2 - 4*a10*c);
xa10(1) = ( -b + da10 ) / (2*a10);
xa10(2) = ( -b - da10 ) / (2*a10);  %solve Equation (22) quadratic in R
dma10=2*max(real(xa10)); % avg diameter in m.

pkla11=La11(find(ab11==max(ab11(770:end)),1,'last')); %nm
hma11=find(ab11>max(ab11(770:end))/2,1,'last');
fwhma11=2*(La11(hma11)-pkla11); %assuming symmetry of peak? in nm
a11=Eg-((h*C)/(pkla11*10^-9)); %energy of peak from wavelength in J
xa11 = zeros(2,1);
da11 = sqrt(b^2 - 4*a11*c);
xa11(1) = ( -b + da11 ) / (2*a11);
xa11(2) = ( -b - da11 ) / (2*a11);  %solve Equation (22) quadratic in R
dma11=2*max(real(xa11)); % avg diameter in m.

pkla12=La12(find(ab12==max(ab12(810:end)),1,'last')); %nm
hma12=find(ab12>max(ab12(810:end))/2,1,'last');
fwhma12=2*(La12(hma12)-pkla12); %assuming symmetry of peak? in nm
a12=Eg-((h*C)/(pkla12*10^-9)); %energy of peak from wavelength in J
xa12 = zeros(2,1);
da12 = sqrt(b^2 - 4*a12*c);
xa12(1) = ( -b + da12 ) / (2*a12);
xa12(2) = ( -b - da12 ) / (2*a12);  %solve Equation (22) quadratic in R
dma12=2*max(real(xa12)); % avg diameter in m.

pkla13=La13(find(ab13==max(ab13(910:end)),1,'last')); %nm
hma13=find(ab13>max(ab13(910:end))/2,1,'last');
fwhma13=2*(La13(hma13)-pkla13); %assuming symmetry of peak? in nm
a13=Eg-((h*C)/(pkla13*10^-9)); %energy of peak from wavelength in J
xa13 = zeros(2,1);
da13 = sqrt(b^2 - 4*a13*c);
xa13(1) = ( -b + da13 ) / (2*a13);
xa13(2) = ( -b - da13 ) / (2*a13);  %solve Equation (22) quadratic in R
dma13=2*max(real(xa13)); % avg diameter in m.

pkla=[pkla1,pkla2,pkla3,pkla4,pkla5,pkla6,pkla7,pkla8,pkla9,pkla10,pkla11,pkla12,pkla13]; %in nm
dma=[dma1,dma2,dma3,dma4,dma5,dma6,dma7,dma8,dma9,dma10,dma11,dma12,dma13]; %in m
fwhma=[fwhma1,fwhma2,fwhma3,fwhma4,fwhma5,fwhma6,fwhma7,fwhma8,fwhma9,fwhma10,fwhma11,fwhma12,fwhma13]; %in nm

%**************************************************************************
% 1. Plot absorbance spectra
s1a=figure;
plot(La1,ab1);
hold on
plot(pkla1,max(ab1),'or');
plot(La1(hma1),ab1(hma1),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 1');
hold off
print(s1a,'-dpng','cdseasp1');

s2a=figure;
plot(La2,ab2);
hold on
plot(pkla2,max(ab2(600:end)),'or');
plot(La2(hma2),ab2(hma2),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 2');
hold off
print(s2a,'-dpng','cdseasp2');

s3a=figure;
plot(La3,ab3);
hold on
plot(pkla3,max(ab3(435:end)),'or');
plot(La3(hma3),ab3(hma3),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 3');
hold off
print(s3a,'-dpng','cdseasp3');

s4a=figure;
plot(La4,ab4);
hold on
plot(pkla4,max(ab4(485:end)),'or');
plot(La4(hma4),ab4(hma4),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 4');
hold off
print(s4a,'-dpng','cdseasp4');

s5a=figure;
plot(La5,ab5);
hold on
plot(pkla5,max(ab5(580:end)),'or');
plot(La5(hma5),ab5(hma5),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 5');
hold off
print(s5a,'-dpng','cdseasp5');

s6a=figure;
plot(La6,ab6);
hold on
plot(pkla6,max(ab6(600:end)),'or');
plot(La6(hma6),ab6(hma6),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 6');
hold off
print(s6a,'-dpng','cdseasp6');

s7a=figure;
plot(La7,ab7);
hold on
plot(pkla7,max(ab7(650:end)),'or');
plot(La7(hma7),ab7(hma7),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 7');
hold off
print(s7a,'-dpng','cdseasp7');

s8a=figure;
plot(La8,ab8);
hold on
plot(pkla8,max(ab8(770:end)),'or');
plot(La8(hma8),ab8(hma8),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 8');
hold off
print(s8a,'-dpng','cdseasp8');

s9a=figure;
plot(La9,ab9);
hold on
plot(pkla9,max(ab9(810:end)),'or');
plot(La9(hma9),ab9(hma9),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 9');
hold off
print(s9a,'-dpng','cdseasp9');

s10a=figure;
plot(La10,ab10);
hold on
plot(pkla10,max(ab10(750:end)),'or');
plot(La10(hma10),ab10(hma10),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 10');
hold off
print(s10a,'-dpng','cdseasp10');

s11a=figure;
plot(La11,ab11);
hold on
plot(pkla11,max(ab11(770:end)),'or');
plot(La11(hma11),ab11(hma11),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 11');
hold off
print(s11a,'-dpng','cdseasp11');

s12a=figure;
plot(La12,ab12);
hold on
plot(pkla12,max(ab12(810:end)),'or');
plot(La12(hma12),ab12(hma12),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 12');
hold off
print(s12a,'-dpng','cdseasp12');

s13a=figure;
plot(La13,ab13);
hold on
plot(pkla13,max(ab13(910:end)),'or');
plot(La13(hma13),ab13(hma13),'or');
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('CdSe nanoparticle absorption 13');
hold off
print(s13a,'-dpng','cdseasp13');

%**************************************************************************
% 2. Tabulate emission peak wavelengths, estimate FWHM, and compute average
% nanocrystal diameter (from peak wavelength)
pkle1=Le1(find(em1==max(em1(900:end)),1)); %nm
hme1=find(em1>max(em1(900:end))/2,1,'last');
fwhme1=2*(Le1(hme1)-pkle1); %assuming symmetry of peak? in nm
e1=Eg-((h*C)/(pkle1*10^-9)); %energy of peak from wavelength in J
xe1 = zeros(2,1);
de1 = sqrt(b^2 - 4*e1*c);
xe1(1) = ( -b + de1 ) / (2*e1);
xe1(2) = ( -b - de1 ) / (2*e1); %solve Equation (22) quadratic in R
dme1=2*max(real(xe1)); % avg diameter in m.

% Do I take real part or modulus? Ditto

pkle2=Le2(find(em2==max(em2),1)); %nm
hme2=find(em2>max(em2)/2,1,'last');
fwhme2=2*(Le2(hme2)-pkle2); %assuming symmetry of peak? in nm
e2=Eg-((h*C)/(pkle2*10^-9)); %energy of peak from wavelength in J
xe2 = zeros(2,1);
de2 = sqrt(b^2 - 4*e2*c);
xe2(1) = ( -b + de2 ) / (2*e2);
xe2(2) = ( -b - de2 ) / (2*e2); %solve Equation (22) quadratic in R
dme2=2*max(real(xe2)); % avg diameter in m.

pkle3=Le3(find(em3==max(em3),1)); %nm
hme3=find(em3>max(em3)/2,1,'last');
fwhme3=2*(Le3(hme3)-pkle3); %assuming symmetry of peak? in nm
e3=Eg-((h*C)/(pkle3*10^-9)); %energy of peak from wavelength in J
xe3 = zeros(2,1);
de3 = sqrt(b^2 - 4*e3*c);
xe3(1) = ( -b + de3 ) / (2*e3);
xe3(2) = ( -b - de3 ) / (2*e3); %solve Equation (22) quadratic in R
dme3=2*max(real(xe3)); % avg diameter in m.

pkle4=Le4(find(em4==max(em4),1)); %nm
hme4=find(em4>max(em4)/2,1,'last');
fwhme4=2*(Le4(hme4)-pkle4); %assuming symmetry of peak? in nm
e4=Eg-((h*C)/(pkle4*10^-9)); %energy of peak from wavelength in J
xe4 = zeros(2,1);
de4 = sqrt(b^2 - 4*e4*c);
xe4(1) = ( -b + de4 ) / (2*e4);
xe4(2) = ( -b - de4 ) / (2*e4); %solve Equation (22) quadratic in R
dme4=2*max(real(xe4)); % avg diameter in m.

pkle5=Le5(find(em5==max(em5),1)); %nm
hme5=find(em5>max(em5)/2,1,'last');
fwhme5=2*(Le5(hme5)-pkle5); %assuming symmetry of peak? in nm
e5=Eg-((h*C)/(pkle5*10^-9)); %energy of peak from wavelength in J
xe5 = zeros(2,1);
de5 = sqrt(b^2 - 4*e5*c);
xe5(1) = ( -b + de5 ) / (2*e5);
xe5(2) = ( -b - de5 ) / (2*e5); %solve Equation (22) quadratic in R
dme5=2*max(real(xe5)); % avg diameter in m.

pkle6=Le6(find(em6==max(em6),1)); %nm
hme6=find(em6>max(em6)/2,1,'last');
fwhme6=2*(Le6(hme6)-pkle6); %assuming symmetry of peak? in nm
e6=Eg-((h*C)/(pkle6*10^-9)); %energy of peak from wavelength in J
xe6 = zeros(2,1);
de6 = sqrt(b^2 - 4*e6*c);
xe6(1) = ( -b + de6 ) / (2*e6);
xe6(2) = ( -b - de6 ) / (2*e6); %solve Equation (22) quadratic in R
dme6=2*max(real(xe6)); % avg diameter in m.

pkle7=Le7(find(em7==max(em7),1)); %nm
hme7=find(em7>max(em7)/2,1,'last');
fwhme7=2*(Le7(hme7)-pkle7); %assuming symmetry of peak? in nm
e7=Eg-((h*C)/(pkle7*10^-9)); %energy of peak from wavelength in J
xe7 = zeros(2,1);
de7 = sqrt(b^2 - 4*e7*c);
xe7(1) = ( -b + de7 ) / (2*e7);
xe7(2) = ( -b - de7 ) / (2*e7); %solve Equation (22) quadratic in R
dme7=2*max(real(xe7)); % avg diameter in m.

pkle8=Le8(find(em8==max(em8),1)); %nm
hme8=find(em8>max(em8)/2,1,'last');
fwhme8=2*(Le8(hme8)-pkle8); %assuming symmetry of peak? in nm
e8=Eg-((h*C)/(pkle8*10^-9)); %energy of peak from wavelength in J
xe8 = zeros(2,1);
de8 = sqrt(b^2 - 4*e8*c);
xe8(1) = ( -b + de8 ) / (2*e8);
xe8(2) = ( -b - de8 ) / (2*e8); %solve Equation (22) quadratic in R
dme8=2*max(real(xe8)); % avg diameter in m.

pkle9=Le9(find(em9==max(em9),1)); %nm
hme9=find(em9>max(em9)/2,1,'last');
fwhme9=2*(Le9(hme9)-pkle9); %assuming symmetry of peak? in nm
e9=Eg-((h*C)/(pkle9*10^-9)); %energy of peak from wavelength in J
xe9 = zeros(2,1);
de9 = sqrt(b^2 - 4*e9*c);
xe9(1) = ( -b + de9 ) / (2*e9);
xe9(2) = ( -b - de9 ) / (2*e9); %solve Equation (22) quadratic in R
dme9=2*max(real(xe9)); % avg diameter in m.

pkle10=Le10(find(em10==max(em10),1)); %nm
hme10=find(em10>max(em10)/2,1,'last');
fwhme10=2*(Le10(hme10)-pkle10); %assuming symmetry of peak? in nm
e10=Eg-((h*C)/(pkle10*10^-9)); %energy of peak from wavelength in J
xe10 = zeros(2,1);
de10 = sqrt(b^2 - 4*e10*c);
xe10(1) = ( -b + de10 ) / (2*e10);
xe10(2) = ( -b - de10 ) / (2*e10); %solve Equation (22) quadratic in R
dme10=2*max(real(xe10)); % avg diameter in m.

pkle11=Le11(find(em11==max(em11),1)); %nm
hme11=find(em11>max(em11)/2,1,'last');
fwhme11=2*(Le11(hme11)-pkle11); %assuming symmetry of peak? in nm
e11=Eg-((h*C)/(pkle11*10^-9)); %energy of peak from wavelength in J
xe11 = zeros(2,1);
de11 = sqrt(b^2 - 4*e11*c);
xe11(1) = ( -b + de11 ) / (2*e11);
xe11(2) = ( -b - de11 ) / (2*e11); %solve Equation (22) quadratic in R
dme11=2*max(real(xe11)); % avg diameter in m.

pkle12=Le12(find(em12==max(em12),1)); %nm
hme12=find(em12>max(em12)/2,1,'last');
fwhme12=2*(Le12(hme12)-pkle12); %assuming symmetry of peak? in nm
e12=Eg-((h*C)/(pkle12*10^-9)); %energy of peak from wavelength in J
xe12 = zeros(2,1);
de12 = sqrt(b^2 - 4*e12*c);
xe12(1) = ( -b + de12 ) / (2*e12);
xe12(2) = ( -b - de12 ) / (2*e12); %solve Equation (22) quadratic in R
dme12=2*max(real(xe12)); % avg diameter in m.

pkle13=Le1(find(em13==max(em13),1)); %nm
hme13=find(em13>max(em13)/2,1,'last');
fwhme13=2*(Le13(hme13)-pkle13); %assuming symmetry of peak? in nm
e13=Eg-((h*C)/(pkle13*10^-9)); %energy of peak from wavelength in J
xe13 = zeros(2,1);
de13 = sqrt(b^2 - 4*e13*c);
xe13(1) = ( -b + de13 ) / (2*e13);
xe13(2) = ( -b - de13 ) / (2*e13); %solve Equation (22) quadratic in R
dme13=2*max(real(xe13)); % avg diameter in m.

pkle=[pkle1,pkle2,pkle3,pkle4,pkle5,pkle6,pkle7,pkle8,pkle9,pkle10,pkle11,pkle12,pkle13]; %in nm
dme=[dme1,dme2,dme3,dme4,dme5,dme6,dme7,dme8,dme9,dme10,dme11,dme12,dme13]; %in m
fwhme=[fwhme1,fwhme2,fwhme3,fwhme4,fwhme5,fwhme6,fwhme7,fwhme8,fwhme9,fwhme10,fwhme11,fwhme12,fwhme13]; %in nm

%**************************************************************************
% 1. Plot emission spectra
s1e=figure;
plot(Le1,em1);
hold on
plot(pkle1,max(em1(900:end)),'or');
plot(Le1(hme1),em1(hme1),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 1');
hold off
print(s1e,'-dpng','cdseesp1');

s2e=figure;
plot(Le2,em2);
hold on
plot(pkle2,max(em2),'or');
plot(Le2(hme2),em2(hme2),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 2');
hold off
print(s2e,'-dpng','cdseesp2');

s3e=figure;
plot(Le3,em3);
hold on
plot(pkle3,max(em3),'or');
plot(Le3(hme3),em3(hme3),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 3');
hold off
print(s3e,'-dpng','cdseesp3');

s4e=figure;
plot(Le4,em4);
hold on
plot(pkle4,max(em4),'or');
plot(Le4(hme4),em4(hme4),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 4');
hold off
print(s4e,'-dpng','cdseesp4');

s5e=figure;
plot(Le5,em5);
hold on
plot(pkle5,max(em5),'or');
plot(Le5(hme5),em5(hme5),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 5');
hold off
print(s5e,'-dpng','cdseesp5');

s6e=figure;
plot(Le6,em6);
hold on
plot(pkle6,max(em6),'or');
plot(Le6(hme6),em6(hme6),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 6');
hold off
print(s6e,'-dpng','cdseesp6');

s7e=figure;
plot(Le7,em7);
hold on
plot(pkle7,max(em7),'or');
plot(Le7(hme7),em7(hme7),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 7');
hold off
print(s7e,'-dpng','cdseesp7');

s8e=figure;
plot(Le8,em8);
hold on
plot(pkle8,max(em8),'or');
plot(Le8(hme8),em8(hme8),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 8');
hold off
print(s8e,'-dpng','cdseesp8');

s9e=figure;
plot(Le9,em9);
hold on
plot(pkle9,max(em9),'or');
plot(Le9(hme9),em9(hme9),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 9');
hold off
print(s9e,'-dpng','cdseesp9');

s10e=figure;
plot(Le10,em10);
hold on
plot(pkle10,max(em10),'or');
plot(Le10(hme10),em10(hme10),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 10');
hold off
print(s10e,'-dpng','cdseesp10');

s11e=figure;
plot(Le11,em11);
hold on
plot(pkle11,max(em11),'or');
plot(Le11(hme11),em11(hme11),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 11');
hold off
print(s11e,'-dpng','cdseesp11');

s12e=figure;
plot(Le12,em12);
hold on
plot(pkle12,max(em12),'or');
plot(Le12(hme12),em12(hme12),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 12');
hold off
print(s12e,'-dpng','cdseesp12');

s13e=figure;
plot(Le13,em13);
hold on
plot(pkle13,max(em13),'or');
plot(Le13(hme13),em13(hme13),'or');
xlabel('Wavelength (nm)');
ylabel('Emission');
title('CdSe nanoparticle emission 13');
hold off
print(s13e,'-dpng','cdseesp13');

%**************************************************************************
% 3. Tabulate nanocrystal diameters in correspondence with the absorption
% peak positions and a description of the color

pks=cat(2,pilau.',pkle.'); %nm
csvwrite('peaks',pks);

fwhm=cat(2,fwhma.',fwhme.'); %nm
csvwrite('FWHM',fwhm);

dm=cat(2,vertcat(dma.',mean(dma),std(dma)),vertcat(dme.',mean(dme),std(dme))); %m
csvwrite('diam',dm);

% Note that emission peaks are not mentioned in this question. I just threw
% the emission values in for the halibut

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

% Plot nanocrystal diameters in correspondence with the absorption peak
% positions and a description of the sample color. You should use energy
% (or frequency) units (e.g. eV or cm^-1) not wavelength. Fit a curve
% through the points. Compare the result to what would be expected for a
% particle in a box.

% Fit to particle in a box
x0=[1;m0/2];
epib=@(a,x)(h^2*a(1)^2)./(8*a(2)*x.^2)*(6.24150934*10^18); %eV
z=lsqcurvefit(epib,x0,dma,(h*C*6.24150934*10^18)./(pkla*10^(-9)));

n=z(1)
mc=z(2)

di=min((dma*10^9)):range((dma*10^9))/20:max((dma*10^9));
for i=1:length(di);
    fitted(i)=(h^2*z(1)^2)./(8*z(2)*(di(i)*10^-9).^2)*(6.24150934*10^18);
end

% Plot diameter vs absorption energy with particle-in-a-box fit
pib=figure;
plot((h*C*6.24150934*10^18)./(pkla*10^(-9)),dma*10^9,'.'); %diameter in nm vs peak position in eV
hold on
plot(fitted,di,'r');
title('CdSe nanoparticle diameter versus absorption energy');
xlabel('Absorption energy (eV)');
ylabel('Nanoparticle diameter (nm)');
hold off
print(pib,'-dpng','PiBx');

%Need to properly do the dimensional analysis.

%**************************************************************************
% 4. Calculate the distribution in the diameter of the nanocrystals for one
% sample based on the FWHM value and assuming a Gaussian distribution for
% the spectra band.

% The FWHM is $\Gamma=2.354*\sigma$. Solve this for $\sigma$ from the 
% numbers in #2

SDa=fwhma/2.354;
SDe=fwhme/2.354;

dis=cat(2,SDa.',SDe.'); %nm
csvwrite('distrib',dis);

%**************************************************************************
% 5. Plot the relation between the absoption peak wavelength and emission
% peak wavelength

% For this I just make a scatter plot of abs vs emis and add a linear fit.
% I expect the constant term to correspond to the shift energy.

p1=polyfit(pkle,pkla,1);
pval=polyval(p1,pkle);
shift=p1(2); %nm

abem=figure;
plot(pkle,pkla,'.');
hold on
plot(pkle,pval,'-r');
title('CdSe nanoparticle absorption wavelength versus emission wavelength');
xlabel('Emission wavelength (nm)');
ylabel('Absorption wavelength (nm)');
hold off
print(abem,'-dpng','abvem');

csvwrite('shft',shift);

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
