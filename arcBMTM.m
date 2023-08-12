%This document is meant to model arcB knockout strain in E. coli BW25113 heat generation curve
clc
close all
clear all

syms X
%arcB growth Gompertz regression
YM = 0.4094;
Y0 = 0.000007933;
K = 1.21;

OD = YM*(Y0/YM)^(exp(-K*X));

%arcB heat generation graph, sum of 2 gaussian
Amp1 =0.0003852;
Amp2 =0.001576;
Mean1 =2.226;
Mean2 =3.966;
SD1 =0.4421;
SD2 =1.112;

Y1 = Amp1*exp(-0.5*(((X-Mean1)/SD1)^2));
Y2 = Amp2*exp(-0.5*(((X-Mean2)/SD2)^2));
HG = Y1 + Y2;


%arcB microaerobic IC [ATP] graph - lognormal curve
A = 22.16;
GeoMean = 5.559;
GeoSD = 1.841;

%model equation is Y=(A/X)*exp(-0.5*(ln(X/GeoMean)/ln(GeoSD))^2)
ATPic = ((A/X)*exp(-0.5*(log(X/GeoMean)/log(GeoSD))^2))*(1/1000);

mu = diff(OD,X);


%generating graphs
figure (1)
fplot(OD,[0,10])
hold on 
fplot(mu,[0,10])
title("OD600 vs Time[h]")
hold off

figure(2)
fplot(HG,[0,10])
hold on
title("Heat Generated [W] vs Time [h]")
hold off

figure(3)
fplot(ATPic,[0,10])
hold on 
title("IC [ATP] vs Time [h]")
hold off

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


%Now to find zero of h1 to find maximum of HG
x0 = 1;
tpeak = fzero(@(X) - (7268017165041563*exp(-((125*X)/139 - 1983/556)^2/2)*((15625*X)/19321 - 247875/77284))/4611686018427387904 - (7105685817192919*exp(-((10000*X)/4421 - 22260/4421)^2/2)*((100000000*X)/19545241 - 222600000/19545241))/18446744073709551616,x0);

%now we have tpeak, so we can plug in and find HGmax
HGfun = @(X)(7268017165041563*exp(-((125*X)/139 - 1983/556)^2/2))/4611686018427387904 + (7105685817192919*exp(-((10000*X)/4421 - 22260/4421)^2/2))/18446744073709551616;
arcBHG_max = HGfun(tpeak);

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
arcB_heatinefficiency = arcBHG_max / maxHG;


%finding intracellular atp at this time
molATPfun = @(X) (277*exp(-(40564819207303340847894502572032*log((1000*X)/5559)^2)/30218920771581602135127064162089))/(12500000000000000000*X);
ATP_cmax = molATPfun(tpeak);

%calculating flux/intracellular [ATP]
D = (ATPflux/ATP_cmax)*arcB_heatinefficiency;
%D is our correction factor

%Now we simply multiply our correction factor D in with our heat
%generation, to get our curve, correcting to W for units

correctedheatrelease = (D*heatrelease)*1000; %multiplied by 1000 to get from W to mW
HG2 = HG*1000;

figure(4)
fplot(correctedheatrelease, [0,10],'black','linestyle','--', 'linewidth', 2)
hold on
%title("{\it \DeltaarcB} Heat Generation")
fplot(HG2, [0,10],'r','linewidth', 2)
set(gca,'FontSize',24)
%legend( "ThermoGen Model","2-Gaussian Equation")
xlim([0 10])
ylim([0 2])
%ylabel("Heat Flow [mW]")
%xlabel("Time [h]")
hold off


%define HG as a function to get correlation coeff
fHG = HGfun;
fcorrectedheatrelease = @(X) (9631524560795680777368169159096849027*(2859560737405723/147573952589676412928)^exp(-(121*X)/100)*exp(-(40564819207303340847894502572032*log((1000*X)/5559)^2)/30218920771581602135127064162089))/(1208925819614629174706176000000000000*X);
n = 10;

fHGlist = [];
fcorrectedheatreleaselist = [];


%define an array of function outputs
for i = 1.8527777781:0.00277777778:4.498
    fHGlist = [fHGlist, fHG(i)];
    fcorrectedheatreleaselist = [fcorrectedheatreleaselist, fcorrectedheatrelease(i)];
end



%print lists of all heat generation outputs from 2 gaussian fitted equation
%(fHGlist) and from MTM model (fcorrectedheatreleaselist)
fprintf('2Gauss HG prediction: [%s]\n', join(string(fHGlist), ','));

fprintf('MTM prediction: [%s]\n', join(string(fcorrectedheatreleaselist), ','));





%% Section 2: AUC Analysis by Integration of Heat Generation Every 30min increment of the IC [ATP] and HG curves to prove coincident peaks

%initialize lists to store HG and ATP AUC (area under the curve) values
HGauc = [];
ATPauc = [];

timelist = [0.5 : 0.5 : 10];
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


