%%%%%% PRELIMINARY ROCKET DESIGN %%%%%%
%% Constants
g0 = 9.81; % Gravitational acceleration [m/s^2]
R_Earth = 6371e03; % Earth radius
mu = 3.986e14; % Gravitational constant for the Earth
C_D = 0.04; % approx drag coeff, independent of Re

mPL = 1000; % Payload mass [kg]

%% Falcon 9
% First stage
mp1_F9 = 395700; % Propellant mass [kg]
ms1_F9 = 25600; % Empty mass [kg]
D1_F9 = 3.7; % Diameter [m]
Isp_1_F9 = 283; % Specific impulse SL [s]

% Second stage
mp2_F9 = 92670; % Propellant mass [kg]
ms2_F9 = 3900; % Empty mass [kg]
D2_F9 = 3.7; % Diameter [m]
Isp_2_F9 = 348; % Specific impulse VAC [s]

% Speed increment calculations
m01_F9 = mPL + ms2_F9 + mp2_F9 + ms1_F9 + mp1_F9; % Initial mass first stage[kg]
mf1_F9 = m01_F9 - mp1_F9; % Final mass first stage [kg]
m02_F9 = mPL + ms2_F9 + mp2_F9; % Initial mass second stage [kg]
mf2_F9 = m02_F9 - mp2_F9; % Final mass second stage[kg]
Delta_V1_F9 = Isp_1_F9 * g0 * log(m01_F9/mf1_F9);
Delta_V2_F9 = Isp_2_F9 * g0 * log(m01_F9/mf1_F9);

Delta_Vtot_F9 = Delta_V1_F9 + Delta_V2_F9; %Total speed increment due to thrust [m/s]

%% Saturn V
% First stage
mp1_SV = 2149500; % Propellant mass [kg]
ms1_SV = 130570; % Empty mass [kg]
D1_SV = 10.06; % Diameter [m]
Isp_1_SV = 263; % Specific impulse SL [s]

% Second stage
mp2_SV = 451650; % Propellant mass [kg]
ms2_SV = 41590; % Empty mass [kg]
D2_SV = 10.06; % Diameter [m]
Isp_2_SV = 390; % Specific impulse VAC [s]

% Speed increment calculations
m01_SV = mPL + ms2_SV + mp2_SV + ms1_SV + mp1_SV; % Initial mass first stage[kg]
mf1_SV = m01_SV - mp1_SV; % Final mass first stage [kg]
m02_SV = mPL + ms2_SV + mp2_SV; % Initial mass second stage [kg]
mf2_SV = m02_SV - mp2_SV; % Final mass second stage[kg]
Delta_V1_SV = Isp_1_SV * g0 * log(m01_SV/mf1_SV);
Delta_V2_SV = Isp_2_SV * g0 * log(m01_SV/mf1_SV);

Delta_Vtot_SV = Delta_V1_SV + Delta_V2_SV; %Total speed increment due to thrust [m/s]

%% Losses (calculated in launch.m)
Delta_V_air = -80.7963088635336;
Delta_V_grav = -5754.55567737017;

%% THE DELTA V WE CAN OBTAIN AT THE MOST IS:
% Falcon 9
Delta_V_ach_F9 = Delta_Vtot_F9 + Delta_V_air + Delta_V_grav
% Saturn V (2 stages)
Delta_V_ach_SV = Delta_Vtot_SV + Delta_V_air + Delta_V_grav

H_max_F9 = (mu / Delta_V_ach_F9^2 - R_Earth)/1000;
H_max_SV = (mu / Delta_V_ach_SV^2 - R_Earth)/1000;