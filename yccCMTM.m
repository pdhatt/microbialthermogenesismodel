%This document is meant to model arcB heat generation curve
clc
close all
clear all

syms X
%yccC growth Gompertz regression microaerobic
YM = 0.4522;
Y0 = 4.904e-006;
K = 1.274;

%model for growth
OD = YM*(Y0/YM)^(exp(-K*X));

%yccC heat generation graph, sum of 2 gaussian
Amp1 =0.001628;
%Amp2 =0.0009023;
Mean1 =3.382;
%Mean2 =3.490;
SD1 =1.385;
%SD2 =0.6767;

%model for heat generation (HG)
Y1 = Amp1*exp(-0.5*(((X-Mean1)/SD1)^2));
%Y2 = Amp2*exp(-0.5*(((X-Mean2)/SD2)^2));
HG = Y1;


%yccC microaerobic IC [ATP] graph - lognormal curve
A = 12.07;
GeoMean = 4.409;
GeoSD = 1.545;

%model equation is Y=(A/X)*exp(-0.5*(ln(X/GeoMean)/ln(GeoSD))^2)
ATPic = (A/X)*exp(-0.5*(log(X/GeoMean)/log(GeoSD))^2)*(1/1000);




%generating graphs
figure (1)
fplot(OD,[0,10])


%Now we can start calculating the data for our model

%First we need to find the OD at each point which we already have eqn
%Then, we need to convert this to cell #, realizing that the data was
%collected for 10mL of cell sample
cell = OD*8*10^8*10; %number of cells

%We also need to convert [ATP] to moles ATP assuming that 1fL volume per
%cell
molATP = ATPic*10^(-15); %mol ATP/cell

%Also need some standard values for conversion of ATP hydrolysis, 
%this model assumes all heat release is from ATP hydrolysis
dH_ATP = 30543.2; %J/mol
mol = 6.022*10^23; %molecules

%now for each time point, we need to convert ATP molecules -> heat released
ATPheat = molATP*dH_ATP; %J/cell

%Convert this J/cell to account for all cells in volume
heatrelease = ATPheat*cell; %J




%Calculating correction factor
%Our model assumes that all intracellular ATP is lysed as heat, but we need
%a correction factor to convert intracellular ATP --> ATP flux
%This is a rough approximation of the model

%ATP flux was calculated using iJO1366
% ATPflux = 78.818; %mmol/gdcw/h
ATPflux = 2.18939e-17; %molATP/cell/s

%now we need to calculate the peak of our heat generation
h1 = diff(HG,X);

%{ 
Just to check to make sure graph looks correct
figure(5)
fplot(h1,[0,10])
%}

%Now to find zero of h1 to find maximum of HG
x0 = 1;
tpeak = fzero(@(X) -(1876956209499947*exp(-((200*X)/277 - 3382/1385)^2/2)*((40000*X)/76729 - 135280/76729))/1152921504606846976,x0);

%now we have tpeak, so we can plug in and find HGmax
HGfun = @(X) (1876956209499947*exp(-((200*X)/277 - 3382/1385)^2/2))/1152921504606846976;
yccCHG_max = HGfun(tpeak);

%This value is the maximum heat generation for this strain
%We can compare this to the maximum if all ATP flux in the cell was heat
%releasing at the rate of ATP hydrolysis (-30.5kJ/mol)
ATP_max = ATPflux*dH_ATP; %J/cell/s == W/cell

%We can also find out the number of cells at the heat peak, since we
%already calculated the time to heat peak as tpeak
ODfun = @(X) YM*(Y0/YM)^(exp(-K*X));
heatpeakcells = ODfun(tpeak)*8*10^8*10;

%Now we can find the maximum heat generation using these parameters, if all
%ATP flux in the cell is lysed and wasted
maxHG = heatpeakcells*ATP_max;

%now we can find our inefficiency score for this cell (how much of its ATP
%flux is being wasted as heat
yccC_heatinefficiency = yccCHG_max / maxHG


%finding intracellular atp at this time
molATPfun = @(X) (1207*exp(-(162259276829213363391578010288128*log((1000*X)/4409)^2)/61413774137875918409655166395489))/(100000000000000000000*X);
ATP_cmax = molATPfun(tpeak);

%calculating flux/intracellular [ATP]
D = (ATPflux/ATP_cmax)*yccC_heatinefficiency;
%D is our correction factor



%Now we simply multiply our correction factor D in with our heat
%generation, to get our curve, correcting to W for units

correctedheatrelease = (D*heatrelease)*1000;%multiply by factor of 1000
%to get from W to mW
HG2 = 1000*HG;

figure(2)
fplot(correctedheatrelease, [0,10],'black','linestyle','--','linewidth', 2)
hold on
%title("{\it \DeltayccC} Heat Generation")
fplot(HG2, [0,10],'color', [1, 100/255, 0/255], 'linewidth', 2)
set(gca,'FontSize',24)
%legend( "ThermoGen Model","2-Gaussian Equation")
xlim([0 10])
ylim([0 2])
%ylabel("Heat Flow [mW]")
%xlabel("Time [h]")
hold off




%attempt to check error with sum of residuals
residualsq = sqrt((correctedheatrelease - HG)^2);
%{
figure(5)
fplot(residualsq, [0,10])
hold on
title("Residual Error plot")
hold off
%}

%Sum residuals
f1 = @(X) (((1876956209499947.*exp(-((200.*X)./277 - 3382/1385).^2./2))./1152921504606846976 - (5982458096879051460429997944910762549*(1600403944050803/147573952589676412928).^exp(-(637.*X)./500).*exp(-(162259276829213363391578010288128.*log((1000.*X)./4409).^2)/61413774137875918409655166395489))./(773712524553362671811952640000000000000.*X)).^2).^(1/2);
sol = integral(f1,0,10);




