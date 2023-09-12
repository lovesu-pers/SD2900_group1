close all
clear
clc

set(groot,"defaultFigurePosition", [50,50,900,600])
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot, 'DefaultAxesFontSize',15);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'DefaultAxesLineWidth', 1)
set(groot, 'DefaultLineLineWidth', 2)


%% Reference parameters
T0 = 288.15;
p0 = 101325;
rho0 = 1.225;
R_earth = 6371*1e3;
g0 = 9.81;


%%
%  Initial state vector
V_init = 0;
gamma_init = pi/2;
X_init = 0;
H_init = 0;

U0 = [V_init; gamma_init; X_init; H_init];

% Stage diameters
d_SV = 14.6;
d_1 = d_SV;
d_2 = d_SV;


massfrac_v = [0.3061, 0.55*0.3450];  % mf/m0
m0 = [2822,642]*1e3;            % Initial mass
mf = m0.*massfrac_v;            % Final mass
mp = m0 - m0.*massfrac_v;       % Propellant mass
m_empt = m0 - mp;               

T_sv = [1*34500,1.5*6500]*1e3;      % Stage thrust
I_SP_V = [1.1*263, 0.99*390];            % Stage I_SP

% Mass flow for each stage
mdot_v = T_sv./(9.8*I_SP_V);

Delta_V = -g0*I_SP_V.*log(massfrac_v); % Delta V for each stage 
Delta_V_TOT = sum(Delta_V)

t_turn = 5; % Start of gravity turn after launch in sec
t_bo = (m0-m0.*massfrac_v) ./ (T_sv./(g0*I_SP_V));  % Burn time for each stage
t_sep = [2,2];                                      % Separation time


AV = [pi*(d_1/2)^2,pi*(d_1/2)^2];   % Stage cross sections

C_D = 0.035; % approx drag coeff, independent of Re

%% Simulation parameters
dt_max = 1e-2; % Max timestep before grav.turn
tspan_in = [0, t_turn]; % Inital simulation interval befor turn

tspan = [t_turn+dt_max, 5*sum(t_bo)]; % Time interval for rest of launch

% ODE45 options
opts = odeset('MaxStep',dt_max); % Fine mesh for pre turn launch
opts_main = odeset('Reltol',5e-14,'AbsTol',1e-14,'Stats','on'); % Let ODE78 choose step size

% Runs the two simulation sections
tic % Used to check performance of solvers, can be ignored
[t_in,U_in] = ode45(@(t,U) odefun(t,U,C_D,AV,t_bo,mdot_v,m0,T_sv,t_turn,t_sep), ...
    tspan_in, U0, opts);
[t_main,U_main] = ode78(@(t,U) odefun(t,U,C_D,AV,t_bo,mdot_v,m0,T_sv,t_turn,t_sep), ...
    tspan, U_in(end,:), opts_main);
toc % Used to check performance of solvers, can be ignored

t = [t_in; t_main];
U = [U_in; U_main];

% Results
V_R = U(:,1);
gamma = U(:,2); % Radians limited to [0,pi]
gamma(1) = pi/2; % Fix
X_R = U(:,3);
H_R = U(:,4);

X_R0 = (R_earth./(R_earth+H_R));
a_R = gradient(V_R,t);

% Delta V calculations 
delta_V_air = -cumtrapz(t,C_D*AV(1).*atmos(H_R,12).*V_R.^2. ...
    /(2*mass_func(t,t_bo,mdot_v,m0,t_sep))); % Cumulative integral
delta_V_air_tot = -trapz(t, C_D*AV(1).*atmos(H_R,12) ... 
    .*V_R.^2./(2*mass_func(t,t_bo,mdot_v,m0,t_sep)));

delta_V_grav = -cumtrapz(t,gfunc(H_R).*sin(gamma));
delta_V_grav_tot = -trapz(t,gfunc(H_R).*sin(gamma));


rhot = zeros(1,length(H_R)); % Density during launch as function time
mt = zeros(1,length(H_R));   % Rocket mass during launch as function time
Tt = zeros(1,length(H_R));  % Thrust during launch as function time
for j = 1:length(H_R)
    rhot(j) = atmos(H_R(j),12);
    mt(j) = mass_func(t(j),t_bo,mdot_v,m0,t_sep);
    Tt(j) = thrust_func(t(j),t_bo,T_sv,t_sep);
end

% Dynamic pressure during launch as function time
q_R = 0.5*V_R.^2.*rhot';

R = R_earth+H_R;
theta = X_R./(R);

cart_x = R.*sin(theta);
cart_y = R.*cos(theta);
%% Entire launch profile
figure(2)
plot(t,X_R./1000,t,H_R./1000)
xlabel('$t$ [s]')
ylabel('$X, \ H $ [km]')
xline((t_bo(1)))
xline(sum(t_bo),'--')
legend('$X$','$H$','Burn-out stage 1', 'Burn-out stage 2', 'Location','northwest')



%% Misc data
figure(3)
set(gcf,'Position',[100,100,900,600])
subplot(6,1,1)
plot(t, gamma*180/pi)
ylabel('Angle $\gamma$ [$^{\mathrm{o}}$]')
xline((t_bo(1)))
xline(sum(t_bo),'--')


subplot(6,1,2)
plot(t, [0;diff(gamma*180/pi)])
ylabel('$\dot{\gamma}$ [$^{\mathrm{o}}$/s]')
xline((t_bo(1)))
xline(sum(t_bo),'--')


subplot(6,1,3)
plot(t,X_R./1000)
ylabel('$X$ [km]')
xline((t_bo(1)))
xline(sum(t_bo),'--')


subplot(6,1,4)
plot(t,H_R./1000)
ylabel('$H$ [km]')
xline((t_bo(1)))
xline(sum(t_bo),'--')

V_req = sqrt((R_earth+H_R).*gfunc(H_R) )/1e3;
subplot(6,1,5)
plot(t, V_R / 1000)
hold on
plot(t, V_req, '-k')
hold off
ylabel('$V$ [km/s]')

xline((t_bo(1)))
xline(sum(t_bo),'--')

subplot(6,1,6)
plot(t,q_R./1000)
xlabel('$t$ [s]')
ylabel('$q$ [kPa]')
xline((t_bo(1)))
xline(sum(t_bo),'--')

axis tight

%% Profile for first 500 km in arc length
figure(4)
plot(X_R./1000,H_R./1000)
xlabel('$x$ [km]')
ylabel('$y$ [km]')
xlim([0,500])


%% Air density and gravity during launch
figure(5)
subplot(2,1,1)
plot(t, rhot)
ylabel('$\rho$ [kg/m$^3$]')
xline(sum(t_bo))
subplot(2,1,2)
plot(t, gfunc(H_R))
ylabel('$g$ [m/s$^2$]')
xlabel('$t$ [s]')
xline((t_bo(1)))
xline(sum(t_bo),'--')

%% Delta V losses and total expended Delta V
figure(6)
subplot(3,1,1)
plot(t, delta_V_air./1e3)
ylabel('$\Delta V_\mathrm{drag}$ [km/s]')

subplot(3,1,2)
plot(t, delta_V_grav./1e3)
ylabel('$\Delta V_\mathrm{grav}$ [km/s]')

subplot(3,1,3)
plot(t,V_R-(delta_V_grav+delta_V_air))
yline(Delta_V_TOT,'-')
xlabel('$t$ [s]')

%% Mass, thrust to weight, and accel plots

figure(7)
subplot(4,1,1)
plot(t, mt/1e3)
ylabel('$m$ [tons]')
xline((t_bo(1)))
xline(sum(t_bo),'--')
legend('$m(t)$', 'Burn-out stage 1', 'Burn-out stage 2', 'Location','northeast')
subplot(4,1,2)
plot(t, Tt./1e6)
ylabel('$T$ [MN]')
subplot(4,1,3)
plot(t, Tt./(mt.*gfunc(H_R).'))
ylabel('$T/W$ [-]')

subplot(4,1,4)
plot(t, a_R./g0)
ylabel('$a/g_0$ [-]')
xlabel('$t$ [s]')

%% Inital turn plots
figure(8)
subplot(4,1,1)
plot(t, gamma*180/pi)
xlim([0, 2*t_turn])
xline(t_turn)
ylabel('Angle $\gamma$ [$^{\mathrm{o}}$]')

subplot(4,1,2)
plot(t, gradient(gamma,t))
xlim([0, 2*t_turn])
xline(t_turn)
ylabel('$\dot{\gamma}$ [$^{\mathrm{o}}$/s]')

subplot(4,1,3)
plot(t, V_R)
xlim([0, 2*t_turn])
xline(t_turn)
ylabel('$V$ [m/s]')

subplot(4,1,4)
plot(t, H_R)
ylabel('$H$ [m]')
xlim([0, 2*t_turn])
xline(t_turn)
xlabel('$t$ [s]')


%% Orbit in cartesian coords
xx = [cart_x, cart_x]./1e3;
yy = [cart_y, cart_y]./1e3;
zz = -ones(size(xx));

figure(9)
% plot(cart_x./1e3,cart_y./1e3, 'Color', 'r', 'LineWidth',2);

surf(xx,yy,zz,[t,t],'LineWidth',3,'EdgeColor','interp',FaceColor='none')
circ([0,0], R_earth./1e3, 3);




xlabel('$x_{cart}$ [km]')
ylabel('$y_{cart}$ [km]')
legend('Orbit', 'Earth')
axis equal
view(2)
colormap("parula")
axis padded

%%
test = (1/V_R)*(gfunc(H_R)-V_R.^2./(R_earth+H_R)) .* cos(gamma);
figure(10)
plot(t,test)
%%
function m = mass_func(t,t_bov,mdot_v,m0,t_sep)
    if t<t_bov(1)
        m = m0(1) - t*mdot_v(1);
     elseif (t_bov(1)<t) & (t<(+t_bov(1)+t_sep(1)))
        m = m0(2);
    elseif (t_bov(1)+t_sep(1)<t) & (t<(t_bov(2)+t_bov(1)+t_sep(1)))
        m = m0(2) - (t-t_sep(1)-t_bov(1))*mdot_v(2);
    else
        m = m0(2) - (t_bov(2))*mdot_v(2);
    end
    
end

function T = thrust_func(t,t_bo,T_sv,t_sep)
    if t<t_bo(1)
        T = T_sv(1);
    elseif t<(t_bo(1)+t_sep(1))
        T = 0;
    elseif t<(t_bo(2)+t_bo(1)+t_sep(1))
        T = T_sv(2);

%     elseif t>(t_bov(2)+t_bov(1)+t_sep(1))
%         T = 0;
    else
        T = 0;
    end
end

function g =  gfunc(H)
    R_earth = 6371*1e3;
    G = 6.6743e-11;
    M = 5.972e24;
    g = G*M./(R_earth + H).^2;
%     g = 9.81*(R_earth/(R_earth+H)).^2;
end

function dUdt = odefun(t,U,C_D,AV,t_bo,mdot_v,m0v,T_sv,t_turn,t_sep)
    % initiates gravity turn

    
    dUdt = zeros(4,1);
    V = U(1);
    gamma = U(2);

%     if t>1.7*(sum(t_bo+t_sep))
%         gamma = 0;
%     end

    X = U(3);
    H = U(4);

    m = mass_func(t,t_bo,mdot_v,m0v,t_sep);
    
    rho = atmos(H,12);
    g = gfunc(H);
    R_earth = 6371*1000;
    T = thrust_func(t,t_bo,T_sv,t_sep);
    
    if t<t_bo(1)
        A = AV(1);
    elseif (t_bo(1)+t_sep(1)<t) && t < (t_bo(2)+t_bo(1)+t_sep(1))
        A = AV(2);
    else
        A = AV(2);
    end
    D = 0.5*rho*C_D*A*V.^2;

    dUdt(1) = T/m - D/m - g*sin(gamma);
    if t<t_turn
        dUdt(2) = 0;
    elseif t==t_turn
        dUdt(2) = -1.e-1;
    else
        dUdt(2) = -(1/V) * (g - (V^2 / (R_earth + H)) ) * cos(gamma);
    end
    
    
    dUdt(3) = V*cos(gamma);
    dUdt(4) = V*sin(gamma);
    dUdt(isinf(dUdt)|isnan(dUdt)) = 0;
end

function C = circ(cntr, r, dim)
    tau = linspace(0,2*pi);
    x = r*cos(tau)+cntr(1);
    y = r*sin(tau)+cntr(2);
    if dim == 2
        C = patch(x,y,'blue');
    else
        C = patch([x,x],[y,y],zeros(size([x,x])));
    end
end