%%%%%% PRELIMINARY ROCKET DESIGN %%%%%%
%% Constants
g0 = 9.81; % Gravitational acceleration [m/s^2]
R_Earth = 6371e03; % Earth radius
T_Earth = 24*3600; % Period of Earth's rotation [s]
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

%% Losses 
% Calculated in launch_sphere.m
Delta_V_air = -41.2334859977086;
Delta_V_grav = -853.943621262776;

% Preliminary data from slides
Delta_V_grav_p = -1200;
Delta_V_air_p = -300;
%% THE DELTA V WE CAN OBTAIN AT THE MOST IS:
% Falcon 9
Delta_V_ach_F9 = Delta_Vtot_F9 + Delta_V_air + Delta_V_grav;
% Saturn V (2 stages)
Delta_V_ach_SV = Delta_Vtot_SV + Delta_V_air + Delta_V_grav;

%% NECESSARY VELOCITY FOR TARGET ORBIT
H_target = 1000e03; %1000km
R_target = H_target + R_Earth;

V_target = sqrt(mu/R_target);

if (Delta_V_ach_F9 > V_target)
    disp('Direct launch is possible for Falcon 9')
else 
    disp('Direct launch is not possible for Falcon 9')
end 

if (Delta_V_ach_SV > V_target)
    disp('Direct launch is possible for Saturn V')
else 
    disp('Direct launch is not possible for Saturn V')
end

%% MISSION ITERATIONS
H_parking = linspace (100, 600, 26)' * 10^3;
R_parking = H_parking(1) + R_Earth;

Mission_info = zeros(size(H_parking,1) + 1, 7);

V_EarthROT = 2 * pi / (T_Earth) * R_Earth;
R_parking = H_parking(1) + R_Earth;

E0 = -mu/R_Earth;
Ef = Mission_info(1,2)^2/2 - mu/R_parking;
DeltaE = abs(Ef - E0);

% Direct launch
Mission_info(1,1) = 0;                                          % Height of parking orbit
Mission_info(1,2) = V_target;                                   % Speed of parking orbit
Delta_V_launch = roots([0.5, Mission_info(1,2), -DeltaE]);
Mission_info(1,3) = Delta_V_launch(Delta_V_launch>0);           % Delta_V necessary to reach parking orbit
Mission_info(1,4) = 0;                                          % First Hohmann impulse
Mission_info(1,5) = 0;                                          % Second Hohmann impulse
Mission_info(1,6) = sum(Mission_info(1,4:5));                   % Total Delta_V required for the Hohmann transfer
Mission_info(1,7) = sum(Mission_info(1,3:5));                   % Total Delta_V required for the whole mission


% Launches to parking orbit

% Paloma: this is obviously not correct. The benefit to launch first to a
% parking orbit is that the gravity losses are not accounted for in the
% transfers, therefore gravity loss is smaller. I need to change that.

Mission_info(2:size(Mission_info,1),1) = H_parking(:);

for i = 2:size(Mission_info,1)
    R_parking = H_parking(i-1) + R_Earth;
    Mission_info(i,2) = sqrt(mu/R_parking);
   
    E0 = -mu/R_Earth;
    Ef = Mission_info(i,2)^2/2 - mu/R_parking;
    DeltaE = abs(Ef - E0);
    Delta_V_launch = roots([0.5, Mission_info(i,2), -DeltaE]);
    Mission_info(i,3) = Delta_V_launch(Delta_V_launch>0); 
    % Launch estimation (without Earth rotation)
    
    % Hohmann transfer
    a_t = 0.5 * (R_parking + R_Earth + R_target);
    V_A = sqrt(2*mu/R_parking - mu/a_t);
    V_B = sqrt(2*mu/R_target - mu/a_t);
    Mission_info(i,4) = V_A - Mission_info(i,2);
    Mission_info(i,5) = V_B - V_target;
    Mission_info(i,6) = sum(Mission_info(i,4:5));
    Mission_info(i,7) = sum(Mission_info(i,3:5));
end

%% File export
H = Mission_info(:,1)';
V = Mission_info(:,2)';
Delta_V = Mission_info(:,3)';
Delta_Vt1 = Mission_info(:,4)';
Delta_Vt2 = Mission_info(:,5)';
Delta_VtT = Mission_info(:,6)';
Delta_VTOT = Mission_info(:,7)';

A = [H;V;Delta_V;Delta_Vt1; Delta_Vt2; Delta_VtT; Delta_VTOT];
fileID = fopen('Mission_analysis.txt','w');
FS_header = '%-6s %-10s %-10s %-10s %-10s %-10s %-10s\r\n';
fprintf(fileID, FS_header, 'H', 'V', 'Delta_V', 'Delta_Vt1', 'Delta_Vt2', 'Delta_VtT', 'Delta_VTOT');
FS_data = '%-6.0f %-.4e %-.4e %-.4e %-.4e %-.4e %-.4e \r\n';
fprintf(fileID, FS_data, A);

%% Actually useful stuff
% Speed difference
Speed_diff_F9 = V_target/Delta_Vtot_F9;
% 40% less propellant mass and structural mass
mp1 = (1-0.6) * mp1_F9;
mp2 = (1-0.6) * mp2_F9;
ms1 = (1-0.8) * ms1_F9;
ms2 = (1-0.8) * ms2_F9;

mf1 = ms1 + ms2 + mp2 + mPL;
mf2 = mPL + ms2;
m01 = mp1 + ms1 + mp2 + ms2 + mPL;
m02 = ms2 + mp2 + mPL;

% I assume the Isp to be the same of Falcon 9
% Isp1 = Isp_1_F9;
% Isp2 = Isp_2_F9;

%
Isp1 = 235;
Isp2 = 330;

% Calculation of the Delta_V we can get with this
Delta_V1_rocket = Isp1 * g0 * log(m01/mf1);
Delta_V2_rocket = Isp2 * g0 * log(m02/mf2);

Delta_V_rocket = Delta_V1_rocket + Delta_V2_rocket;
%% Nominal characteristics Falcon 9 (again)
mPL_F9 = 22800;

m01_F9_nom = ms1_F9 + ms2_F9 + mp1_F9 + mp2_F9 + mPL_F9;
mf1_F9_nom = ms1_F9 + ms2_F9 + mp2_F9 + mPL_F9;
m02_F9_nom = ms2_F9 + mp2_F9 + mPL_F9;
mf2_F9_nom = ms2_F9 + mPL_F9;

mu1_F9 = mf1_F9_nom / m01_F9_nom;
mu2_F9 = mf2_F9_nom / m02_F9_nom;

SR1_F9 = ms1_F9/m01_F9_nom;
SR2_F9 = ms2_F9/m02_F9_nom;


%% Applied to Miura 5
m01_M5 = 68742 - 900;
mf1_M5 = mu1_F9 * m01_M5;

ms1_M5 = SR1_F9 * m01_M5;
m02_M5 = mf1_M5 - ms1_M5;
mf2_M5 = m02_M5 * mu2_F9;
ms2_M5 = SR2_F9 * m02_M5;
mPL_M5_calc = mf2_M5 - ms2_M5;

mp1_M5 = m01_M5 - mf1_M5;
mp2_M5 = m02_M5 - mf2_M5;

PL_rate = mPL/mPL_M5_calc;
ms1_M5 + ms2_M5 + mp1_M5 + mp2_M5 + mPL

% Now I calculate the rocket parameters again, taking into account that the
% Miura is "supposed" to be able to carry nearly three times the mass that
% is required from us.

mp1 = 0.2 * mp1_M5;
mp2 = 0.2 * mp2_M5;

ms1 = 0.1 * ms1_M5;
ms2 = 0.1 * ms2_M5;

% Calculation of Isp of the two stages bc it is not given to us:
tbo1_M5 = 182;
tbo2_M5 = 420;
T1_M5_TOT = 950*10^3; % Total thrust for 5 TEPREL-C regular engines
T1_M5 = T1_M5_TOT/5;
T2_M5 = 50*10^3; % One TEPREL-C vacuum optimized engine

m_dot1_M5 = (mp1_M5/5)/tbo1_M5; % Assuming all five rockets are the same and consume the same
m_dot2_M5 = mp2_M5/tbo2_M5;

Isp1 = (T1_M5/m_dot1_M5)/g0;
Isp2 = (T2_M5/m_dot2_M5)/g0;
%THIS IS NOT CORRECT. CHECK



