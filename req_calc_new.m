%%%%%% PRELIMINARY ROCKET DESIGN %%%%%%
%% Constants
g0 = 9.81; % Gravitational acceleration [m/s^2]
R_Earth = 6371e03; % Earth radius
T_Earth = 24*3600; % Period of Earth's rotation [s]
mu = 3.986e14; % Gravitational constant for the Earth
C_D = 0.04; % approx drag coeff, independent of Re

mPL = 1000; % Payload mass [kg]

%% Falcon 9 information:
mPL_F9 = 22800; % Maximum payload [kg]

mp1_F9 = 395700; % Propellant mass [kg]
ms1_F9 = 25600; % Empty mass [kg]
D1_F9 = 3.7; % Diameter [m]
Isp_1_F9 = 283; % Specific impulse SL [s]

mp2_F9 = 92670; % Propellant mass [kg]
ms2_F9 = 3900; % Empty mass [kg]
D2_F9 = 3.7; % Diameter [m]
Isp_2_F9 = 348; % Specific impulse VAC [s]

m01_F9 = mPL_F9 + ms2_F9 + mp2_F9 + ms1_F9 + mp1_F9; % Initial mass first stage[kg]
mf1_F9 = m01_F9 - mp1_F9; % Final mass first stage [kg]
m02_F9 = mPL_F9 + ms2_F9 + mp2_F9; % Initial mass second stage [kg]
mf2_F9 = m02_F9 - mp2_F9; % Final mass second stage[kg]

%% Preliminary mass calculations
m01_M5 = 68742; % Nominal initial mass for the Miura 5

% This includes a nominal payload mass of 900kg to LEO launching from the
% equator. If we launch a heavier mass, the structure will also be heavier
% and more propellant will be required.

% Initial estimation:
m01 = 69000; 

% Mass ratios are based on Falcon 9, as not enough info is available for
% Miura 5
mu1_F9 = mf1_F9_nom / m01_F9_nom;
mu2_F9 = mf2_F9_nom / m02_F9_nom;

% First stage:
mf1 = mu1_F9 * m01;
mp1 = m01-mf1;
ms1 = 4100; % Chosen based on other rockets

% Second stage:
m02 = mf1 - ms1;
mf2 = mu2_F9 * m02;
mp2 = m02-mf2;
ms2 = mf2 - mPL;

% Verification of the total mass
MTOT = mp1 + ms1 + ms2 + mp2 + mPL;
Error = abs(MTOT - m01)/m01 * 100;

%% Engine characteristics
% First option (easy one): Use Space X Merlin for the first stage and
% TEPREL-B

% Second option (more "complex" one): Use both TEPREL-C as the Miura
%% Delta V calculations
Delta_V1 = Isp1 * g0 * log(m01/mf1);
Delta_V2 = Isp2 * g0 * log(m02/mf2);
Delta_Vtot = Delta_V1 + Delta_V1;
