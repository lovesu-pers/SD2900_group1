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


%%
T0 = 288.15;
p0 = 101325;
rho0 = 1.225;
R_earth = 6371*1e3;
alt = 0:100:90e3;

%%
X_init = 1e-8;
H_init = 1e-8;
Xdot_init = 0;
Hdot_init = 0;
downrange_start = 0;
U0 = [X_init,Xdot_init,H_init,Hdot_init];

g0 = 9.81;

d_SV = 14.6;
d_1 = d_SV*1;
d_2 = d_SV*1;


massfrac_v = [0.3061, 0.3450];
m0 = [2822,500]*1e3;

mf = m0 - m0.*massfrac_v;

m_empt = m0 - mf;

T_sv = 1*[34500,6500]*1e3;
I_SP_V = [263, 390];

mdot_v = T_sv./(9.8*I_SP_V);

Delta_V = -g0*I_SP_V.*log(1./mf);
Delta_V_TOT = sum(Delta_V)

t_turn = 3;
t_bo = (m0-m0.*massfrac_v) ./ (T_sv./(g0*I_SP_V));
t_sep = [2,2];

AV = [pi*(d_1/2)^2,pi*(d_1/2)^2];

C_D = 0.05;

%%
dt_max = 1e-2;
tspan_in = [0, t_turn];
tspan = [t_turn+dt_max, 2*sum(t_bo)];
opts = odeset('MaxStep',dt_max);
opts_main = odeset('MaxStep',1);
% opts_main = odeset('Stats','off') ;

[t_in,U_in] = ode45(@(t,U) odefun(t,U,C_D,AV,t_bo,mdot_v,m0,T_sv,t_turn,t_sep), tspan_in, U0, opts);
[t_main,U_main] = ode45(@(t,U) odefun(t,U,C_D,AV,t_bo,mdot_v,m0,T_sv,t_turn,t_sep), tspan, U_in(end,:), opts_main);
t = [t_in; t_main];
U = [U_in; U_main];

test = U(:,2)./sqrt(U(:,2).^2+U(:,4).^2);
gamma = acos(U(:,2)./sqrt(U(:,2).^2+U(:,4).^2));
X_R = U(:,1);
H_R = U(:,3);
V_R = sqrt(U(:,2).^2 + U(:,4).^2);
X_R0 = (R_earth./(R_earth+H_R));

a_R = gradient(V_R);

delta_V_air = -cumtrapz(t,C_D*AV(1).*atmos(H_R,12).*V_R.^2./(2*mass_func(t,t_bo,mdot_v,m0,t_sep)));
delta_V_air_tot = -trapz(t, C_D*AV(1).*atmos(H_R,12).*V_R.^2./(2*mass_func(t,t_bo,mdot_v,m0,t_sep)));
delta_V_grav = -cumtrapz(t,gfunc(H_R));
delta_V_grav_tot = -trapz(t,gfunc(H_R));

rho_launch = zeros(1,length(H_R));
mt = zeros(1,length(H_R));
Tt = zeros(1,length(H_R));
for j = 1:length(H_R)
    rho_launch(j) = atmos(H_R(j),12);
    mt(j) = mass_func(t(j),t_bo,mdot_v,m0,t_sep);
    Tt(j) = thrust_func(t(j),t_bo,T_sv,t_sep);
    j
end

q_R = 0.5*V_R.^2.*rho_launch';
%%
figure(2)
plot(t,X_R./1000,t,H_R./1000)
xlabel('$t$ [s]')
ylabel('$X, \ H $ [km]')
xline((t_bo(1)))
xline(sum(t_bo),'--')
legend('$X$','$H$','Burn-out stage 1', 'Burn-out stage 2', 'Location','northwest')




figure(3)
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




figure(4)
plot(U(:,1)./1000,H_R./1000)
xlabel('$x$ [km]')
ylabel('$y$ [km]')
xlim([0,500])



figure(5)
subplot(2,1,1)
plot(t, rho_launch)
ylabel('$\rho$ [kg/m$^3$]')
xline(sum(t_bo))
subplot(2,1,2)
plot(t, gfunc(H_R))
ylabel('$g$ [m/s$^2$]')
xlabel('$t$ [s]')
xline((t_bo(1)))
xline(sum(t_bo),'--')

%%
figure(6)
subplot(3,1,1)
plot(t, delta_V_air./1e3)
ylabel('$\Delta V_\mathrm{drag}$ [km/s]')

subplot(3,1,2)
plot(t, delta_V_grav./1e3)
ylabel('$\Delta V_\mathrm{grav}$ [km/s]')

subplot(3,1,3)
plot(t,V_R-(delta_V_grav+delta_V_air))
yline(Delta_V_TOT,'-o')
xlabel('$t$ [s]')

%%

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
plot(t, a_R)
ylabel('$a$ [m/s$^2$]')
xlabel('$t$ [s]')

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
%     if abs(t_turn-t)<1e-5
%         U(2) = 100e-8;
%     end
    if t==t_turn
        U(2) = 1e-1;
        U(1) = 1e-8;
    end
    dUdt = zeros(4,1);
    m = mass_func(t,t_bo,mdot_v,m0v,t_sep);
    cosgamma = U(2)/sqrt(U(2)^2+U(4)^2);
    singamma = U(4)/sqrt(U(2)^2+U(4)^2);
    singamma = fillmissing(singamma,"constant",1);
    cosgamma = fillmissing(cosgamma,"constant",0);
    rho = atmos(U(3),12);
    g = gfunc(U(3));
    R_earth = 6371*1000;
    T = thrust_func(t,t_bo,T_sv,t_sep);
    
    if t<t_bo(1)
        A = AV(1);
    elseif (t_bo(1)+t_sep(1)<t) && t < (t_bo(2)+t_bo(1)+t_sep(1))
        A = AV(2);
    else
        A = AV(2);
    end
    
    if t>t_bo(2)
        T
    end


    dUdt(1) = U(2);
    dUdt(2) = T*cosgamma/m - 0.5*cosgamma/m*rho*C_D*A*(U(2)^2+U(4)^2)- ... 
        U(2)*U(4)/(R_earth + U(3));
    dUdt(3) = U(4);
    dUdt(4) = T*singamma/m - 0.5*rho*singamma/m*C_D*A*(U(2)^2+U(4)^2)- ... 
        (g-U(2)^2/(R_earth + U(3)));
    dUdt = fillmissing(dUdt,"constant",0);
end