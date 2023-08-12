%This document is meant to model glnL overexpression data
clc
close all
clear all

syms X
%glnL overexpression strain growth Gompertz regression microaerobic
YM = 0.001395;
Y0 = 3.189e-007;
K = 11.223;

%model for growth
OD = YM*(Y0/YM)^(exp(-K*X));


%glnL overexpression strain heat generation graph
Amp1 = 0.001292;
Mean1 = 3.918;
SD1 = 1.702;

%model for heat generation (HG) - SINGLE GAUSSIAN chosen here
HG = Amp1*exp(-0.5*(((X-Mean1)/SD1)^2));


%glnL overexpression strain intracellular (IC) [ATP] graph - lognormal curve
A = 0.02315;
GeoMean = 4.000;
GeoSD = 2.474;

%model equation is Y=(A/X)*exp(-0.5*(ln(X/GeoMean)/ln(GeoSD))^2)
ATPic = (A/X)*exp(-0.5*(log(X/GeoMean)/log(GeoSD))^2)*1/1000;



%generating graphs
figure (1)
fplot(HG,[0,10])


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
tpeak = fzero(@(X) -(5958298335808185*exp(-((500*X)/851 - 1959/851)^2/2)*((250000*X)/724201 - 979500/724201))/4611686018427387904, x0);

%now we have tpeak, so we can plug in and find HGmax
HGfun = @(X) (5958298335808185*exp(-((500*X)/851 - 1959/851)^2/2))/4611686018427387904;
WTHG_max = HGfun(tpeak);

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
WT_heatinefficiency = WTHG_max / maxHG;


%finding intracellular atp at this time
molATPfun = @(X) (463*exp(-(40564819207303340847894502572032*log(X/4)^2)/66570061195607304043194483566881))/(20000000000000000000000*X);
ATP_cmax = molATPfun(tpeak);

%calculating flux/intracellular [ATP]
D = (ATPflux/ATP_cmax)*WT_heatinefficiency;
%D is our correction factor


%Now we simply multiply our correction factor D in with our heat
%generation, to get our curve, correcting to W for units

correctedheatrelease = (D*heatrelease)*1000; %multiply by 1000 to get to mW
HG2 = HG*1000;

figure(2)
fplot(correctedheatrelease, [0,10],'black','linestyle','--','linewidth', 2)
hold on
%title("Wild-Type Heat Generation")
fplot(HG2, [0,10],'black','linewidth', 2)
set(gca,'FontSize',24)
%legend( "ThermoGen Model","2-Gaussian Equation")
xlim([0 10])
%ylim([0 10])
%ylabel("Heat Flow [mW]")
%xlabel("Time [h]")
hold off

%Model standard error measurement
%define HG as a function to get correlation coeff
fHG = HGfun;
fcorrectedheatrelease = @(X) (69699121228378656451759364376730057785539739087071*(2108482682833683/9223372036854775808)^exp(-(11223*X)/1000)*exp(-(40564819207303340847894502572032*log(X/4)^2)/66570061195607304043194483566881))/(13937965749081639463459823920405225941237760000000000*X);
fHGlist = [];
fcorrectedheatreleaselist = [];


%define an array of function outputs
for i = 1.25:0.00277777778:4.5
    fHGlist = [fHGlist, fHG(i)];
    fcorrectedheatreleaselist = [fcorrectedheatreleaselist, fcorrectedheatrelease(i)];
end




%disp(fHGlist)


%print lists of all heat generation outputs from 2 gaussian fitted equation
%(fHGlist) and from MTM model (fcorrectedheatreleaselist)
fprintf('2Gauss HG prediction: [%s]\n', join(string(fHGlist), ','));


fprintf('MTM prediction: [%s]\n', join(string(fcorrectedheatreleaselist), ','));



%% Section 2: AUC Analysis by Integration of Heat Generation Every 30min increment of the IC [ATP] and HG curves to prove coincident peaks

%initialize lists to store HG and ATP AUC (area under the curve) values
HGauc = [];
ATPauc = [];

timelist = [2 : 0.5 : 10];
%Define for loop, looping through # of hours: [0:0.5:10]
for i = timelist
    tlast = i-0.5;
    tcurr = i;
    ATPauc = [ATPauc, int(ATPic,tlast,tcurr)]; %find definite integral for AUC
    HGauc = [HGauc, int(HG,tlast,tcurr)]; %find definite integral for AUC
end

ATPauc_rounded = round(ATPauc, 10); %round values to print as list
HGauc_rounded = round(HGauc, 10); %round values to print as list
 
disp("ATPauc_rounded: ")
disp(ATPauc_rounded)
disp("HGauc_rounded: ")
disp(HGauc_rounded)


