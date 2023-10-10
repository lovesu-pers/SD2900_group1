%%%% BALSA WOOD GLIDER CALCULATIONS %%%%
%%%% Design limits:
% Max wingspan: 500 mm 
% Max fuselage length: 500 mm
% Max chord length = 100 mm
% Fuselage building material available: 5 x 15 x 1000 mm
% Lifting surfaces building material available: 2.5 x 100 x 1000 mm
close all; clc; clear
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');  
%% Constants and design parameters
Cl_alpha = 0.11 * 180 / pi; % Sectional lift curve slope of flat plate [rad^-1]
eps_slope = 0.1;            % (\partial \varepsilon) / (\partial \alpha)
etaHT = 1;                 % Tail effectiveness: eta_HT = qt/q

alpha = linspace(0,5,50);            % Range of angles of attack calculated [deg]

m_f = 10.4*1e-3;                     % Mass of entire fuselage [kg]
rho_balsa = m_f/(1*5e-3*15e-3)       % Density of the balsa [kg/m^3]

%% Initial reference surface estimation
% Assume the glider has the highest values for fuselage length and
% wingspan. Rectangular planform shape:

b = 400;            % Wingspan [mm]
AR = 5;             % Aspect ratio [-]; should be between 4.5 and 7.5
c = b/AR;           % Mean chord [mm]
S = b^2/AR;         % Wing area (reference area)

xf = 420;           % Dim. of fuselage x (length)
yf = 15;            % y
zf = 5;             % z

ct = 50;            % Dim. of tail x (chord)
bt = 110;           % y 
zt = 2.5;           % z

xw = 100;           % Dim. of wingspan x
yw = b;             % y
zw = zt;            % z

bVT = 30;
cVT = 30;
tVT = 5;

V_t = ct*bt*zt;     % Tail volume
V_w = xw*yw*zw;     % Wing volume
V_f = xf*yf*zf;     % Fuselage volume
V_HT = ct*bt*zt;
V_VT = bVT*cVT*tVT;

% Centroid calculations
xLE_w = 100+5;        % Position of wing LE from nose      
xj_w = 50;          % Position of wing centroid with reference to wing LE
xLE_t = 420-ct+5;        % Position of tail LE from nose  
xj_t = 20;          % Position of tail centroid with reference to tail LE
xLE_VT = xLE_t;       % Position of vertical tail LE from nose 
xj_VT = cVT/2;         % Poistion of vertical tail centroid with reference to vertical tail LE

xi_t = xLE_t + xj_t;    % Position of wing centroid with reference to nose
xi_w = xLE_w + xj_w;    % Position of tail centroid with reference to nose
xi_VT = xLE_VT + xj_VT; % Position of vertical tail centroid with reference to nose
xi_f = xf/2;            % Position of fuselage centroid

x_cg_w = c/2;            % Centre of gravity of the wing
x_cg_HT = ct/2;

%% Initial estimations for horizontal and vertical tail
S_H =ct*bt;          % Horizontal tail area
S_T = bVT*cVT;        % Vertical tail planform area
bt = S_H/ct;            % Span of the tail
ARt = bt^2/S_H;

% Vol_w = S * 2.5;        % Wing volume
% Vol_HT = S_H * 2.5;     % Horizontal tail volume
% Vol_VT = S_T * 2.5;     % Vertical tail volume
% Vol_f = xf*yf*zf;

% Mass addition
M = 0.005;                  % Added mass [kg]
Vol_M = M/rho_balsa * 10^9;     % Equivalent balsa wood volume of the added mass [mm^3]
xi_M = 10;                  % Position of the added mass with reference to the nose
x_M = xi_M - xLE_w;         % Position of the added mass with reference to wing LE

V_TOT = V_HT + V_w + V_f + V_VT + Vol_M;
x_cg = (xi_w*V_w + xi_t * V_HT + xi_f * V_f + xi_VT * V_VT + x_M * Vol_M)/...
       (V_TOT);    % Centre of Mass with respect to nose
x_cg_bar = x_cg / c;    % Dimensionless centre of mass with respect to nose                 

% Centre of mass with respect to wing LE
x_cg_2 = x_cg - xLE_w;
x_cg_2_bar = x_cg_2/c;

%% Calculation of wing lift coefficient of lifting surfaces
CL_alpha_w = Cl_alpha * AR / (2+sqrt(AR^2 + 4));                  % Lift curve slope of wing eq(9-73)
CL_alpha_t = Cl_alpha * ARt / (2+sqrt(ARt^2 + 4));                % Lift curve slope of tail eq(9-73)

%% Wing Contribution
% Equations (24-67), (24-71)
x_ac_w = c/4;                                           % Aerodynamic centre of the wing
x_ac_w_bar = 1/4;                                       % Dimensionless aerodynamic centre of the wing

CL0_w = 0;                                              % Zero lift coefficient of the wing (symmetrical)
CMac_w = 0;                                             % Moment coefficient of the wing about its AC (symmetrical wing)

CM0_w = CMac_w + CL0_w * (x_cg_w/c - x_ac_w_bar)     % Residual pitch moment coefficient of the wing (symmetrical)
CMalpha_w = CL_alpha_w * (x_cg_w/c - x_ac_w_bar)       % Longitudinal Stability derivative of the wing (- for stability)

% CL_w plot (linear part)
CL_w = CL0_w + CL_alpha_w * deg2rad(alpha);
figure(1)
plot(alpha, CL_w);
xlabel('\alpha (deg)')
ylabel('C_{L_w}')
title('Lift curve (linear section) of the wing')

% CM_w plot
CM_w = CM0_w + CMalpha_w * deg2rad(alpha);
figure(2)
plot(alpha, CM_w);
xlabel('\alpha (deg)')
ylabel('C_{M_w}')
title('Pitch moment curve of the wing')

%% Horizontal Tail Contribution
% Equations (24-68) and (24-73)
x_ac_HT = ct/4 + xLE_t - xLE_w;                         % Aerodynamic centre of the tail
x_ac_HT_bar = x_ac_HT / c;

CL0_HT = 0;                                             % Zero lift coefficient of the wing (symmetrical)
CMac_HT = 0;                                            % Moment coefficient of the wing about its AC (symmetrical wing)   

l_HT_bar = (x_ac_HT - x_cg_2)/c;                        % Distance from CG to AC_HT

CM0_HT = CMac_HT + CL0_HT * (l_HT_bar)                 % Residual pitch moment coefficient of the HT (symmetrical)

VHT = S_H*l_HT_bar/S;                                   % Horizontal tail volume
CMalpha_HT = -etaHT * VHT * CL_alpha_t * (1-eps_slope)   % Longitudinal Stability derivative of the HT (- for stability)

% CL_HT plot (linear part)
CL_HT = CL0_HT + CL_alpha_t * deg2rad(alpha);
figure(3)
plot(alpha, CL_HT);
xlabel('\alpha (deg)')
ylabel('C_{L_{HT}}')
title('Lift curve (linear section) of the horizontal tail')

% CM_HT plot
CM_HT = CM0_HT + CMalpha_HT * deg2rad(alpha);
figure(4)
plot(alpha, CM_HT);
xlabel('\alpha (deg)')
ylabel('C_{M_{HT}}')
title('Pitch moment curve of the horizontal tail')

%% Stick-fixed longitudinal stability
i_w = deg2rad(1);                                           % Wing incidence angle
i_t = deg2rad(-3);                                            % Tail incidence angle

% Aerodynamic centre of glider CHECK THE ORIGIN OF THIS!
x_ac_bar = (x_ac_w_bar + CL_alpha_t/CL_alpha_w * etaHT * S_H/S * x_ac_HT_bar * (1-eps_slope))/...
            (1+CL_alpha_t/CL_alpha_w * etaHT * S_H/S * (1-eps_slope)); 

dfus = 2*sqrt((yf/2)^2 + (zf/2)^2);                 % "Diameter" of the fuselage
K_WB = 1 + 0.025 * (dfus/b) - 0.25 * (dfus/b)^2;    % Fuselage correction factor                 

CL0 = CL0_w + etaHT * S_H/S * CL0_HT;                                                   % Zero lift coefficient of glider
CLalpha = K_WB * CL_alpha_w + etaHT * S_H/S * CL_alpha_t * (1-eps_slope);                   % Slope of lift curve of glider
CM0 = CM0_w + etaHT * VHT * CL_alpha_t * (i_w - i_t)                                     % Residual moment coefficient of glider
CMalpha = CL_alpha_w * (x_cg_2_bar - x_ac_bar) - etaHT * VHT * CL_alpha_t * (1-eps_slope);   % Longitudinal stability of glider
CMalpha_true = CMalpha * pi/180

% CL_TOT plot (linear part)
CL = CL0 + CLalpha * deg2rad(alpha);
figure(5)
plot(alpha, CL);
xlabel('\alpha (deg)')
ylabel('C_{L}')
title('Lift curve (linear section) of the glider')

% CM_TOT plot
CM = CM0 + CMalpha * deg2rad(alpha);
figure(6)
plot(alpha, CM);
xlabel('\alpha (deg)')
ylabel('C_{M}')
title('Pitch moment curve of the glider')
% grid on

figure(7)
axis equal
rectangle('position',[xLE_w, -b/2, c,  b])
hold on
rectangle('position',[0, -15/2, xf,  15])
rectangle('position',[xLE_t, -bt/2, ct,  bt])
plot(x_cg,0,'o','MarkerSize',10)
plot(x_ac_bar*c+xLE_w,0,'o','MarkerSize',10)
plot(xi_M,0,'*','MarkerSize',10)
hold off
legend('cg','AC','Dead mass')

SM = x_ac_bar - x_cg_2_bar           % Static margin
SM_2 = -CMalpha/CLalpha
x_cg_2_bar
alpha_trim = rad2deg(-CM0/CMalpha)

CL_trim = CL0 + CLalpha * deg2rad(alpha_trim);
mtot = rho_balsa*V_TOT*1e-9 + M;
Vtrim = sqrt(2*mtot*9.8/(CL_trim*1.12*S*1e-6))