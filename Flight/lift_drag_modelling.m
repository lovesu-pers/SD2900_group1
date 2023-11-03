%%%%% lift and drag modelling %%%%%

% Determine the maximum lift coefficient
% Generally, the maximum lift coefficient for takeoff is about 80% of the
% maximum lift coefficient for landing
% According to the table of typical values of maximum lift coefficient -
% the coefficient value tange of millitary trainers 
% https://www.ae.utexas.edu/~varghesep/class/aircraft/Suggestions.p
CL = ;
CLmax = 1.8; 
% Maximum lift coefficient normally 1.2-1.8
CLmaxl = 2.5; 
% Landing maximum lift coefficient normally 1.6-2.2 fighters 1.6-2.6
CLmaxTO = 0.8*CLmaxl; 
% Maximun lift coefficient for takeoff normally 1.4-2.0

% Determine the zero-lift resistance coefficient
% The equivalent skin friction coefficient method is used to estimate the
% zero-lift resistance formula
Cfe = 0.0040; 
% based on the equivalent skin friction coefficient table - naval flighter,
% taking subsonic speed 0.0040
k = 2.3;
% Due to the estimate stage, the aircraft design plane area cannot yet be
% determined, so statistical values are used to estimate
% wetted area ratio k temporarily take 2.3 in the estimation stage
Cdo = Cfe*k; %zero-lift resistance coefficient
% During takeoff and landing, flaps and landing gear have a greater impact,
% and deltaCdo depends on the subsequent specific selection

% Determine the maximum lift-to-drag ratio
% At subsonic speed, the aircraft drag coefficient
A = 10.2; % bigger maybe 12?
% wing aspect ratio, take 3.48 as the aspect ratio of supersonic flighter
% aircraft is generally 2-4
LE = 60; % wing leading edge sweep angle
e = 4.61*(1-0.045*A^0.68)*(cos(LE))^0.15-3.1;% Oswald efficiency factor
Cd = Cdo+CLmaxTO^2/(pi*A*e); % aircraft drag coefficient
% The lift-to-drag ratio directly on two design factors, wing span and
% wetted area. At subsonic speeds, the foamula for the maximum lift-to-drag
% ratio can be used to estimate the lift-to-drag ratio L/D
LDmax = 0.5*(pi*A*e/Cdo)^0.5;

c = 0.9880; % chord length(m)
Cmax = 2*c/100;
t = c*12/100; % thickness
x_chamber = c*4/10; % x position of the chamber
T = 269; %temperature K
mu = 1.458e-6*(T^1.5/(110.4+T));% Ns/m2 
Mach = 0.263;
V_cruise = Mach*340.2900; % cruise speed m/s 
rho = 0.909; %kg/m3
Re = (rho * V_cruise*c)/mu; %Reynolds numbers
CDi = (CL^2)/(pi*A*e);
