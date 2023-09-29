%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PRELIMINARY ROCKET DESIGN %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
g0 = 9.81;              % Gravitational acceleration [m/s^2]
R_Earth = 6371e03;      % Earth radius
T_Earth = 24*3600;      % Period of Earth's rotation [s]
mu = 3.986e14;          % Gravitational constant for the Earth
C_D = 0.04;             % approx drag coeff, independent of Re

mPL = 1000;             % Payload mass [kg]

%% Falcon 9 information:
mPL_F9 = 22800;         % Maximum payload [kg]

mp1_F9 = 395700;        % Propellant mass [kg]
ms1_F9 = 25600;         % Empty mass [kg]
D1_F9 = 3.7;            % Diameter [m]
Isp_1_F9 = 283;         % Specific impulse SL [s]

mp2_F9 = 92670;         % Propellant mass [kg]
ms2_F9 = 3900;          % Empty mass [kg]
D2_F9 = 3.7;            % Diameter [m]
Isp_2_F9 = 348;         % Specific impulse VAC [s]

m01_F9 = mPL_F9 + ms2_F9 + mp2_F9 + ms1_F9 + mp1_F9;    % Initial mass first stage[kg]
mf1_F9 = m01_F9 - mp1_F9;                               % Final mass first stage [kg]
m02_F9 = mPL_F9 + ms2_F9 + mp2_F9;                      % Initial mass second stage [kg]
mf2_F9 = m02_F9 - mp2_F9;                               % Final mass second stage[kg]

% mu1_F9 = mf1_F9 / m01_F9;
% mu2_F9 = mf2_F9 / m02_F9;
epsilon1_F9 = ms1_F9/(ms1_F9+mp1_F9); 
epsilon2_F9 = ms2_F9/(ms2_F9+mp2_F9); 

SR1_F9 = ms1_F9/m01_F9; %"Structural ratio" first stage
SR2_F9 = ms2_F9/m02_F9; %"Structural ratio" second stage

%% Engine characteristics
% First option (easy one): Use Space X Merlin for the first stage and
% TEPREL-B

% Second option (more "complex" one): Use both TEPREL-C as the Miura
m01_M5 = 68742;                     % Initial mass Miura 5 (from website)
mf1_M5 = 0.167 * m01_M5;            % Final mass first stage
% ms1_M5 = (SR1_F9+0.009) * m01_M5;   % Structural mass first stage
ms1_M5 = (SR1_F9+0.015) * m01_M5;
mp1_M5 = m01_M5 - mf1_M5;           % Propellant mass first stage

m02_M5 = mf1_M5 - ms1_M5;           % Initial mass second stage
mf2_M5 = m02_M5 * (0.18);           % Final mass second stage
ms2_M5 = (SR2_F9+0.0222) * m02_M5;   % Structural mass second stage
mp2_M5 = m02_M5 - mf2_M5;           % Propellant mass second stage

mPL_M5_calc = mf2_M5 - ms2_M5;      % Nominal payload of Miura 5 900kg to LEO

tbo1_M5 = 182;                      % Burning time first stage (from website)
tbo2_M5 = 420;                      % Burning time second stage (from website)
T1_M5_TOT = 950*10^3;               % Total thrust for 5 TEPREL-C regular engines
T1_M5 = T1_M5_TOT/5;                % One TEPREL-C regular engine
T2_M5 = 50*10^3;                    % One TEPREL-C vacuum optimized engine

m_dot1_M5 = (mp1_M5/5)/tbo1_M5;     % Assuming all five rockets are the same and consume the same
m_dot2_M5 = mp2_M5/tbo2_M5;         % In second stage, there's only one engine

Isp1 = (T1_M5/m_dot1_M5)/g0;        % Isp TEPREL-C
Isp2 = (T2_M5/m_dot2_M5)/g0;        % Isp TEPREL-C vacuum

%% Preliminary mass calculations
m01 = 69000;        % Initial mass estimation        

mu1 = 0.17;         % Mass ratio first stage
mu2 = 0.18;         % Mass ratio second stage

% First stage:
mf1 = mu1 * m01;    % Final mass first stage
mp1 = m01-mf1;      % Propellant mass first stage
ms1 = 4200;         % Chosen based on other rockets

% Second stage:
m02 = mf1 - ms1;    % Initial mass second stage
mf2 = mu2 * m02;    % Final mass second stage
mp2 = m02-mf2;      % Propellant mass second stage
ms2 = mf2 - mPL;    % Structural mass second stage

epsilon1 = ms1/(ms1 + mp1);
epsilon2 = ms2/(ms2 + mp2);
% % Verification of the total mass
% MTOT = mp1 + ms1 + ms2 + mp2 + mPL;
% Error = abs(MTOT - m01)/m01 * 100;

%% Delta V calculations
Delta_V1 = Isp1 * g0 * log(m01/mf1);
Delta_V2 = Isp2 * g0 * log(m02/mf2);
Delta_Vtot = Delta_V1 + Delta_V1;
