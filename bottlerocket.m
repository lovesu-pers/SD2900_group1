clear
clc
close all

set(groot,"defaultFigurePosition", [200,200,1100,800])
set(groot,'defaulttextinterpreter','none');
set(groot,'defaultLegendInterpreter','none');
set(groot, 'DefaultAxesFontSize',20);
set(groot,'defaultAxesTickLabelInterpreter','none');  
%%
ff = 0.4;
rho = 1000;
m0 = 0.115 + 0.02*3 + 0.050;
p0 = 8.013e5;
patm = 1e5;
gamma = 1.4;
V_tot = 1.5*(1e-1)^3;
V_water = V_tot*ff;
V_air = V_water/ff-V_water;

f = 0:0.01:1;

K = p0*(1-ff)*V_tot;

W3 = (K/(1-gamma)) * ( (V_tot^(1-gamma))-( V_tot-V_tot*ff )^(1-gamma) )

W = (p0*(V_air^gamma).*(1-f).^(gamma))/(1-gamma) .* (V_air.^(1-gamma)- ... 
    V_air.^(1-gamma)*(1-f).^(1-gamma));
W2 = p0*(V_tot-V_tot*ff)/(1-gamma) * ((1-ff)^gamma - (1-ff))

Wm = (1./(m0+f*V_air*rho)).*(p0*(V_air^gamma).*(1-f).^(gamma))/(1-gamma) .* (V_air.^(1-gamma)- ... 
    V_air.^(1-gamma)*(1-f).^(1-gamma));
figure(1)
plot(f,W,f,Wm, LineWidth=2)
xlabel('Volume fraction, f [-]')
ylabel('Work/mass [J/kg]')
%%

rho_water = 1000;
Anozzle = pi*(0.01131)^2;
Ek = 0.8*200;
m = m0;
g = 9.8;
alpha = 40*(pi/180);
v_Ek = sqrt(2*Ek/m)
v_eBern = sqrt(2*(p0-patm)/rho_water)

tdis = (V_water*rho_water)^2 / (v_Ek*m*Anozzle*rho_water) 

v0 = v_Ek*[cos(alpha);sin(alpha)];
x0 = [0;0];
cD = 0.6;
rho = 1.225;
A = pi*(0.085/2)^2; 
mu = cD*rho*A/(2*m);

% v_eP = sqrt(2*(p0-patm)/(rho_water*(1-(Anozzle/A)^2)))

v_eP = sqrt( 2 * (p0*((1-ff)./(1-f)).^gamma - patm) ./ (1-(Anozzle/A)^2)  );

mdot = rho_water*Anozzle*v_eP;

mwater = rho_water*ff*V_tot;
% t_dis = mwater/mdot;

% T_tot = mdot*(v_eP+(p0-patm)*Anozzle/mdot)
% T_tot = mdot.*(v_eP);
% T_tot = 2*Anozzle*(p0*()^gamma)
% v0_T = t_dis*T_tot/m; 
H0 = [x0;v0];
tspan = [0,8];

% figure(3)
% plot(f,T_tot)

[t,H] = ode45(@(t,H) odefun(t,H,mu,g), tspan, H0);

figure(2)
plot(H(:,1), H(:,2), 'LineWidth', 3     )
xlim([0, max(H(:,1))+10])
ylim([0, max(H(:,2))+1])


%%
Nalpha = 10;
alpha_v = (pi/180)*linspace(30,50, Nalpha);

figure(5)
xlabel('Downrange [m]')
ylabel('h [m]')

hold on
xlim([0, max(H(:,1))+10])
ylim([0, max(H(:,2))+1])
for i = 1:Nalpha
    alpha = alpha_v(i);
    v0 = v_Ek*[cos(alpha);sin(alpha)];
    H0 = [x0;v0];
    [t,H] = ode45(@(t,H) odefun(t,H,mu,g), tspan, H0);
    plot(H(:,1), H(:,2), 'LineWidth', 3)
    pause(0.6)
end
hold off
xlim([0, max(H(:,1))+10])
ylim([0, max(H(:,2))+1])


function dHdt = odefun(t,H,mu,g)
dHdt = zeros(4,1);
dHdt(1) = H(3);
dHdt(2) = H(4);
dHdt(3) = -mu*H(3)*norm(H(3:4));
dHdt(4) = -g-mu*H(4)*norm(H(3:4));
end