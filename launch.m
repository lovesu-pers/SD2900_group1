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

Clrs = linspecer(7);
%%
T0 = 288.15;
p0 = 101325;
rho0 = 1.225;

alt = 0:100:90e3;

% pratio = zeros(length(alt),1);
% 
% for i = 1:length(alt)
%     pratio(i) = atmos(alt(i),10);
% end

% figure(1)
% loglog(pratio,alt);

%%
X_init = 1e-8;
H_init = 10;
Xdot_init = 100e-8;
Hdot_init = 10;
downrange_start = 0;
U0 = [X_init,Xdot_init,H_init,Hdot_init];

g0 = 9.8;

d_SV = 14.6;
d_1 = d_SV*1;
d_2 = d_SV*1;


massfrac_v = [0.3061, 0.3450];
m0 = [2822,642]*1e3;

T_sv = [34500,4400]*1e3;
I_SP_V = [263, 390];

mdot_v = T_sv./(9.8*I_SP_V);

t_bov = (m0-m0.*massfrac_v) ./ (T_sv./(g0*I_SP_V));

AV = [pi*(d_1/2)^2,pi*(d_1/2)^2];

C_D = 0.04;

%%
tspan = [0, 1.2*sum(t_bov)];
% opts = odeset('MaxStep',300e-2);
opts = odeset('Stats','off') ;
[t,U] = ode45(@(t,U) odefun(t,U,C_D,AV,t_bov,mdot_v,m0,T_sv), tspan, U0, opts);
gamma = acos(U(:,2)./sqrt(U(:,2).^2+U(:,4).^2));
X_R = U(:,1);
H_R = U(:,3);
V_R = sqrt(U(:,2).^2 + U(:,4).^2);

delta_V_air = -cumtrapz(C_D*AV(1).*atmos(H_R,12).*V_R.^2./(2*mass_func(t,t_bov,mdot_v,m0)));
delta_V_air_tot = -trapz(t, C_D*AV(1).*atmos(H_R,12).*V_R.^2./(2*mass_func(t,t_bov,mdot_v,m0)));
delta_V_grav = cumtrapz(gfunc(H_R));
delta_V_grav_tot = -trapz(t,gfunc(H_R));

rho_launch = zeros(1,length(H_R));
for j = 1:length(H_R)
    rho_launch(j) = atmos(H_R(j),12);
    j
end
%%
figure(2)
plot(t,X_R./1000,t,H_R./1000)
xlabel('$t$ [s]')
ylabel('$X, \ H $ [km]')
xline((t_bov(1)))
xline(sum(t_bov),'--')
legend('$X$','$H$','Burn-out', 'Location','northwest')




figure(3)
subplot(6,1,1)
plot(t, gamma*180/pi, Color=Clrs(1,:))
ylabel('Angle $\gamma$ [$^{\mathrm{o}}$]')
xline((t_bov(1)))
xline(sum(t_bov),'--')


subplot(6,1,2)
plot(t, [0;diff(gamma*180/pi)], Color=Clrs(2,:))
ylabel('$\dot{\gamma}$ [$^{\mathrm{o}}$/s]')
xline((t_bov(1)))
xline(sum(t_bov),'--')


subplot(6,1,3)
plot(t,X_R./1000, Color=Clrs(3,:))
ylabel('$X$ [km]')
xline((t_bov(1)))
xline(sum(t_bov),'--')


subplot(6,1,4)
plot(t,H_R./1000, Color=Clrs(4,:))
ylabel('$H$ [km]')
xline((t_bov(1)))
xline(sum(t_bov),'--')


subplot(6,1,5)
plot(t, V_R / 1000, Color=Clrs(7,:))
ylabel('$V$ [km/s]')

xline((t_bov(1)))
xline(sum(t_bov),'--')

subplot(6,1,6)
plot(t,0.5*V_R.^2.*rho_launch./1000)
xlabel('$t$ [s]')
ylabel('$q$ [kPa]')
xline((t_bov(1)))
xline(sum(t_bov),'--')




figure(4)
plot(U(:,1)./1000,H_R./1000)
xlabel('$x$ [km]')
ylabel('$y$ [km]')




figure(5)
subplot(2,1,1)
plot(t, rho_launch)
ylabel('$\rho$ [kg/m$^3$]')
xline(sum(t_bov))
subplot(2,1,2)
plot(t, gfunc(H_R))
ylabel('$g$ [m/s$^2$]')
xlabel('$t$ [s]')
xline((t_bov(1)))
xline(sum(t_bov),'--')
%%
figure(6)
subplot(2,1,1)
plot(t, delta_V_air)
ylabel('$\Delta V_\mathrm{drag}$ [m/s]')

subplot(2,1,2)
plot(t, delta_V_grav)
ylabel('$\Delta V_\mathrm{grav}$ [m/s]')
xlabel('$t$ [s]')
%%
function m = mass_func(t,t_bov,mdot_v,m0v)
    if t<t_bov(1)
        m = m0v(1) - t*mdot_v(1);
    elseif (t_bov(1)<t)<t_bov(2)+t_bov(1)
        m = m0v(2) - t*mdot_v(2);
    else
        m = m0v(2) - t_bov(2)*mdot_v(2);
    end
    
end

function T = thrust_func(t,t_bov,T_sv)
    if t<t_bov(1)
        T = T_sv(1);
    elseif (t_bov(1)<t) && t<(t_bov(2)+t_bov(1))
        T = T_sv(2);
    elseif t>(t_bov(2)+t_bov(1))
        T = 0;
    end
end

function g =  gfunc(H)
    R_earth = 6371*1000;
    G = 6.6743e-11;
    M = 5.972e24;
    g = G*M./(R_earth + H).^2;
end

function dUdt = odefun(t,U,C_D,AV,t_bov,mdot_v,m0v,T_sv)
    tspan = 613;

    dUdt = zeros(4,1);
    m = mass_func(t,t_bov,mdot_v,m0v);
    cosgamma = U(2)/sqrt(U(2)^2+U(4)^2);
    singamma = U(4)/sqrt(U(2)^2+U(4)^2);
    rho = atmos(U(3),12);
    g = gfunc(U(3));
    R_earth = 6371*1000;
    T = thrust_func(t,t_bov,T_sv);
    
    if t<t_bov(1)
        A = AV(1);
    elseif (t_bov(1)<t)<t_bov(2)+t_bov(1)
        A = AV(2);
    else
        A = AV(2);
    end
    
    if t>=t_bov(2)+t_bov(1)
        t
    end
    dUdt(1) = U(2);
    dUdt(2) = T*cosgamma/m - 0.5*cosgamma/m*rho*C_D*A*(U(2)^2+U(4)^2)- ... 
        U(2)*U(4)/(R_earth + U(3));
    dUdt(3) = U(4);
    dUdt(4) = T*singamma/m - 0.5*singamma/m*C_D*A*(U(2)^2+U(4)^2)- ... 
        (g-U(2)^2/(R_earth + U(3)));

end