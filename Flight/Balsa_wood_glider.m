%%%% BALSA WOOD GLIDER CALCULATIONS %%%%
%%%% Design limits:
% Max wingspan: 500 mm 
% Max fuselage length: 500 mm
% Max chord length = 100 mm
% Fuselage building material available: 5 x 15 x 1000 mm
% Lifting surfaces building material available: 2.5 x 100 x 1000 mm

%% Constants and design parameters
Cl_alpha = 0.11 * 180 / pi; % Sectional lift curve slope of flat plate [rad^-1]

%% Initial reference surface estimation
% Assume the glider has the highest values for fuselage length and
% wingspan. Rectangular planform shape:

b = 500;            % Wingspan [mm]
AR = 6;             % Aspect ratio [-]; should be between 4.5 and 7.5
c = b/AR;           % Mean chord [mm]
S = b^2/AR;         % Wing area (reference area)

L_f = 500;          % Length of fuselage

%% Initial estimations for horizontal and vertical tail
S_H = 0.2 * S;      % Horizontal tail area
S_T = 0.5 * S_H;    % Vertical tail planform area

