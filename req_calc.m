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
Delta_V2_F9 = Isp_2_F9 * g0 * log(m02_F9/mf2_F9);

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

Mission_info = zeros(size(H_parking,1) + 1, 7);

V_EarthROT = 2 * pi / (T_Earth) * R_Earth;

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

