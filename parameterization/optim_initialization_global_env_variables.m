function []=optim_initialization_global_env_variables(air_temp,RH,initCO2,initCI,PAR)
%Enzyme activation settings
global EAPPDK;
global PRac 
global RedoxEnyAct;
global GsResponse;
EAPPDK=1;%%if EAPPDK=0 PPDK is fully actived; %%if EAPPDK=1 include PPDK activation by PDRP
PRac=1;%%if PRac=0 Rubisco is fully actived; %%if PRac=1 include Rubisco activation by Rca
RedoxEnyAct=1; %%if RedoxEnyAct=0 activities of photosynthetic enzymes are not regulated by light intensity; %%if RedoxEnyAct=1 include light induced enzyme activation 
GsResponse=1; %%if GsResponse=0 Ball berry model no time dependent Gs response ; %%if GsResponse=1 include Gs response, using ki and kd

%Inpute parameter settings

global MeasuredTemp;
MeasuredTemp=air_temp;


global kdcon;
kdcon=1;%=1: constant kd; =0: kd change with light

global Para_mata;
Para_mata=1;%%if Para_mata=1, C4 Metabolic model and Gs model integrated  if Para_mata=0 Steady-state mdoel and gs model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Environmental parameters
global WeatherTemperature;
global Air_O2;
global Air_CO2;
global WeatherRH;
global WeatherWind;
global Radiation_NIR;
global Radiation_LW;
global Radiation_PARo;
global PhiLeaf


WeatherWind=3.5;%m/s %Set as 3.5 to match the bundary layer conductance from licor data
PhiLeaf=0;%Mpa
WeatherTemperature=air_temp;

Air_O2=210;
Air_CO2=initCO2;

WeatherRH=RH; % 65%
Radiation_NIR=0;
Radiation_LW=0;
Radiation_PARo=PAR;

global O2;
O2=Air_O2*1.26/1000;%O2 concentration 

global CI;
CI=initCI*0.4/(3 * 10^4); %The initial intercelluar CO2 Concnetration mmol/L


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%C4 decarboxylation pathway settings
% global C13ratio;% C insotope simulation
% C13ratio=0;%0.01115;

global U;
global V;
U=0;% light partition coefficient
V=0;% light partition coefficient
global Ratio;
Ratio=4; % Enezyme activity variation factor
global pathway_option;
global RatioPPDK;%PPDK in BSC
RatioPPDK=0;
global Pvpr8M;
Pvpr8M=0;%if 0 Glycerate kinase in BSchl; if 1 Glycerate kinase in Mchl;
global Bchl_CP;%Pi content in BS chloroplast
Bchl_CP= 25.0;
global MC_CP;%Pi content in MC
MC_CP=15.0;
global Mchl_CP;%Pi content in MC chloroplast
Mchl_CP=15.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathway_option=0;
%%% 0=Normol NADP-ME type 
%%% 1=Asp+Mal transport and MDH type 
%%% 2=Asp+Mal and PCK type 
%%% 3 Asp+Mal and PCK+MDH type 
%%% 4 Asp and PCK only type
%%% 6 DiT2 mutant
%%% 7 NAD-ME type
%%% 8 NAD-ME+PCK type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reaction rate record
global TIME_M;
global OLD_TIME_M;
global Meta_VEL;
TIME_M=0;
OLD_TIME_M=0;
Meta_VEL=zeros(1,9);

global TIME_N;
global OLD_TIME;

TIME_N=0;
OLD_TIME=0;
global Gs_VEL;
Gs_VEL=zeros(1,9);

global TIME_K;
global OLD_TIME_K;
global MY_VEL;
TIME_K=0;
OLD_TIME_K=0;
MY_VEL=zeros(1,9);

