%%%%% lift and drag modelling %%%%%

% Determine the maximum lift coefficient
% Generally, the maximum lift coefficient for takeoff is about 80% of the
% maximum lift coefficient for landing
% According to the table of typical values of maximum lift coefficient -
% the coefficient value tange of millitary trainers 
% https://www.ae.utexas.edu/~varghesep/class/aircraft/Suggestions.p
%need CL or we can use the maximum CLmax for analysis?
CLmax = 1.54; 
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
%Mach = 0.263;
V_cruise = 77;
Mach = V_cruise/340.2900;
%V_cruise = Mach*340.2900; % cruise speed m/s 
rho = 0.909; %kg/m3
Re = (rho * 77*c)/mu; %Reynolds numbers
CDimax = (CLmax^2)/(pi*A*e);
Cf = Cfe; % skin friction coeeficient, same as equivalent one before
CDf = Cf*k; % skin friction drag coefficient


%%%%% 3D

%CLalpha = 2*pi*A/(2+sqrt(A^2+4));
%CLmax0 = 2*(W/S)/(rho*(V_cruise^2));
%cbar = S/b;
%Clalpha = CLalpha/0.9;
%Cmalpha = 0.0002*alpha-0.01198;
%alpha_ZL = -Cl0/Clalpha;
%CL0 = abs(alpha_ZL)*Clalpha;
%CM_aw = CLalpha*(Cmalpha/Clalpha);

beta = (1-(Mach^2))^0.5; %Mach number parameter(Prandtl-Glauert)
Clalpha = 2*pi;
kappa = 2*pi; %ratio of 2-dimension lift curve slope to 2*pi
%In general, most airfoils only approximately dispaly the 2*pi lift slope
%as predicted by thin airfoil theory. That is because airfoils are not
%actually infinitely thin in practice, and will deviate from thin airfoil
%theory by a small amount
LAMBDA_c2 = 0; %sweepback of mid-chord 
%The resulting expression below is valid for noncurved planform shapes and
%Mach is not higher than 0.8
CLalpha = (2*pi*A)/(2+sqrt(((A*beta/kappa)^2)*(1+((tan(LAMBDA_c2)^2)/(beta^2)))+4)); %per rad

W = 1321.88; %5880N = 1321.88Ibf
W_S = 12.8237; %614N/m2 = 12.8237Ibf/ft2
S = W/W_S; %ft2
Vs = 0.952*sqrt(2*W/(rho*S*CLmax)); %Wing loading on stalling speed
CLmax0 = CLmax*0.9; %Max lift coefficient of wing in its clean configuration
