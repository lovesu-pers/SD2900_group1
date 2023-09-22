%% Header
clear
close all
clc
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot, 'DefaultAxesFontSize',15);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'DefaultAxesLineWidth', 1)
set(groot, 'DefaultLineLineWidth', 2)

%%%% Reqs: Mapping and Aerospace toolbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants and conversion factors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2r = pi/180;
r2d = 180/pi;
RE = 6371e3;
muE = 3.986e5 * (1e3)^3; %m^3/s^2
g0 = 9.80665;
omegaE = [0;0;7.292115855377074e-5;];

earth_REF = referenceSphere('Earth');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Launch site %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kourou, French Guiana
Nji = 100;
deltagrav_ji = zeros(Nji,1);
h_turn_ji = linspace(10,1e3,Nji);
delta_V_grav_tot = zeros(Nji,1);
delta_V_air_tot = zeros(Nji,1);
for ji = 1:Nji

h0 = 2;
lat0 =5.167713*d2r;
long0 =-52.683994*d2r;
r0 = latlong2cart(lat0, long0, h0);
rmag = norm(r0);
V0 = [0;0;0];
VErot = cross([0;0;muE],r0);
turn_azi = (90 - 26.6821)*d2r;
turn_fp = 89.8*d2r;
turnvec = 1*[cos(turn_fp)*cos(turn_azi); ...
        cos(turn_fp)*sin(turn_azi); ...
        sin(turn_fp)];
[X_turn, Y_turn, Z_turn] = enu2ecef(turnvec(1),turnvec(2),turnvec(3), ... 
    lat0,long0,h0,earth_REF,"radians");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROCKET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nstage = 2;
m0     = [2780000,677000,215000];                % Initail/fueled mass kg
massfraction   = [0.2817,0.3663, 1*0.4930];        % mf/m0
mf     = m0.*massfraction;       % Final/empty mass
mprop = m0-mf;
T0      = [33400e3,1.0*4450e3,1.0*2000e3];            % Thrust N
Isp    = [265,390,421];                 % Specific impulse s
d      = [10,10,6.6];                   % Diameter m

altPO  = h_turn_ji(ji);
tsep   =  1;
TW = T0./(g0*m0);

A0     = pi*(d./2).^2; 
mdot0 = T0./(g0*Isp);
tbo = (m0-mf) ./ mdot0;  % Burn time for each stage
DeltaV = -g0*Isp.*log(massfraction);
DeltaVtot = sum(DeltaV);

tstage_index = [1, 0, tbo(1); 
                2, tbo(1)+tsep, tbo(2)+tbo(1)+tsep;
                3, tbo(2)+tbo(1)+2*tsep, tbo(2)+tbo(1)+tbo(3)+2*tsep];       % Stage number, start time, bo-time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE solving       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfin = 0.5*60*60;
U0 = [r0;V0];

opts_turn = odeset('RelTol',10e-10, 'Stats','on', ...
    'Events',@(t,U) turncond(t,U,altPO)); % Let ODE78 choose step size
opts_main = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @steercond );

opts_steer = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events',@(t,U) cruisecond(t,U,tstage_index(end,3)));

opts_cruise = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events',@crashcond);


% Runs pre-turn launch
tic % Used to check performance of solvers, can be ignored
[t_turn,U_turn] = ode45(@(t,U) ode_turn(t,U,m0,mdot0,tstage_index,tbo,T0,A0,altPO), ...
    [0, 10*60], U0, opts_turn);

% Normal gravity turn
timeatturn = t_turn(end);
Vatturn = norm(U_turn(end,4:6));
V_turn = norm(U_turn(end,4:6)) * ([X_turn; Y_turn; Z_turn] - r0);
[t_asc,U_asc] = ode45(@(t,U) ode_main(t,U,m0,mdot0,tstage_index,tbo,T0,A0,altPO), ...
    [t_turn(end), tfin], [U_turn(end,1:3)';V_turn], opts_main);

% Controlled turn
t_steerstart = t_asc(end);
r_asc = U_asc(end, 1:3);
V_asc = U_asc(end, 4:6);
gamma0 = gammafunc(r_asc,V_asc);
[t_steer,U_steer] = ode45(@(t,U) ode_steer(t,U,m0,mdot0,tstage_index,tbo,T0,A0, gamma0,t_steerstart), ...
    [t_steerstart, tfin], [r_asc';V_asc'], opts_steer);

% Stops control and cruises
[t_cruise,U_cruise] = ode45(@(t,U) ode_main(t,U,m0,mdot0,tstage_index,tbo,T0,A0), ...
    [t_steer(end), tfin], U_steer(end,:), opts_cruise);
toc % Used to check performance of solvers, can be ignored

% Adds the solution for the segments
t_main = [t_turn; t_asc;t_steer;t_cruise];
U_main = [U_turn; U_asc; U_steer;U_cruise];

t_calc = [t_turn; t_asc;t_steer];
U_calc = [U_turn; U_asc;U_steer];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessing      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(t_main);
Ncalc = length(t_calc);

% Date used to convert between ECI and ECEF
tt = [2000 1 1 17 15 10];
epoch = datetime(2000,1,1,17,15,10);

dt = zeros(N,1);
for k = 2:N
    dt(k) =  t_main(k) - t_main(k-1);
end

% Saves results in new vectors
r_res = U_main(:,1:3);
V_res = U_main(:,4:6); 
Vmag_res = vecnorm(V_res,2,2); 
h_res = vecnorm(r_res,2,2) - RE;
latlong_res = zeros(N,2);
a_res = gradient(Vmag_res,t_main);

% Max distance from earth CoM
maxr = max(h_res+RE);

gamma = zeros(N,1);  % Flightpath angle
mres = zeros(N,1);   % Mass as function of time

Mres = zeros(N,1);   % Mach 
CDres = zeros(N,1);  % Drag coeff
rhot = zeros(N,1);
gres = zeros(N,1);
A = zeros(N,1);
Tres = zeros(N,1);

for j = 1:N
    % Lat-Long in ECI
    [latlong_res(j,1), latlong_res(j,2)] = cart2latlong(r_res(j,:));
    
    gamma(j) = gammafunc(r_res(j,:),V_res(j,:));
    
    [mres(j),Tres(j),A(j),CDres(j)] = statefunc(t_main(j), ... 
        m0,mdot0,tstage_index,tbo,T0,A0,V_res(j,:),h0);
    
    Mres(j) = norm(V_res(j,:))/atmos(h_res(j),13);

    rhot(j) = atmos(h_res(j),12);

    gres(j) = gfunc(h_res(j)+RE);
end


dragk = zeros(Ncalc,1);
gk = zeros(Ncalc,1);
for k = 1:Ncalc
    gammak = gammafunc(r_res(k,:),V_res(k,:));
    [mk,Tk,Ak,CDk] = statefunc(t_main(k), ... 
        m0,mdot0,tstage_index,tbo,T0,A0,V_res(k,:),h0);
    
    qk = 0.5*Vmag_res(k).^2.*atmos(h_res(k),12);
    

    gk(k) =  gfunc(h_res(k)+RE)*cos(gammak*d2r);
    dragk(k) = CDk.*Ak.*qk./mk;
end
delta_V_air_tot(ji) = -trapz(t_calc,dragk);
delta_V_grav_tot(ji) = -trapz(t_calc,gk);

% Delta V calculations 
% delta_V_air = -cumtrapz(t_main,CDres.*A.*q_R./mres); % Cumulative integral
% delta_V_grav = -cumtrapz(t_main,gres.*cos(gamma));



ji
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
plot(h_turn_ji, delta_V_grav_tot./1e3)
ylabel('$\Delta V$ loss [km/s]' )
xlabel('Turn alt [m]')

%% Functions 
function dUdt = ode_turn(t,U,m0,mdot0,tstage_index,tbo,T0,A0, altPO)
    RE = 6371e3;

    dUdt = zeros(4,1);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    rho = atmos(h,12); 
    
    [m,T,A,CD] = statefunc(t,m0,mdot0,tstage_index,tbo,T0,A0,V,h);
    
    omegaE = [0;0;7.292115855377074e-5];
    muE = 3.986e5 * (1e3)^3; %m^3/s^2


    rhat = r/norm(r);
    Vhat = V/norm(V);
    Vhat(isinf(Vhat)|isnan(Vhat)) = 0;

    g = gfunc(r);
    if t==0
        T=T*rhat;
    else
        T = T*Vhat;
    end
    drag = -0.5*rho * CD * A * (norm(V-cross([0;0;7.292115855377074e-5],r))).^2 /m * Vhat;
    dUdt(1:3) = V ;
    dUdt(4:6) = T/m + drag/m + g;

end

function dUdt = ode_main(t,U,m0,mdot0,tstage_index,tbo,T0,A0, altPO)
    RE = 6371e3;

    dUdt = zeros(4,1);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    rho = atmos(h,12); 
    
    [m,T,A,CD] = statefunc(t,m0,mdot0,tstage_index,tbo,T0,A0,V,h);
    
    omegaE = 0*[0;0;7.292115855377074e-5;];
    muE = 3.986e5 * (1e3)^3; %m^3/s^2


    rhat = r/norm(r);
    Vhat = V/norm(V);
    Vhat(isinf(Vhat)|isnan(Vhat)) = 0;

    g = gfunc(r);
    if t==0
        T=T*rhat;
    else
        T = T*Vhat;
    end
    drag = -0.5*rho * CD * A * (norm(V)).^2 /m * Vhat;
    dUdt(1:3) = V + cross(omegaE,r);
    dUdt(4:6) = T/m + drag/m + g;
end


function dUdt = ode_steer(t,U,m0,mdot0,tstage_index,tbo,T0,A0, gamma0,t0)
    RE = 6371e3;   
    d2r = pi/180;
    dUdt = zeros(4,1);
    r = U(1:3);
    Vold = U(4:6);
    h = norm(r) - RE;
    gammaold = gammafunc(r,Vold)*d2r;
    gammaend = 0.8*pi/2;
    gamma =  gamma0*d2r + (t-t0)*( (gammaend-gamma0*d2r)/(tstage_index(end,3)-t0)  );
    

    delta = gamma - gammaold;
    
    V = vecrot(r,Vold,delta);
    
    
    rho = atmos(h,12); 
    
    [m,T,A,CD] = statefunc(t,m0,mdot0,tstage_index,tbo,T0,A0,V,h);
    
    omegaE = 0*[0;0;7.292115855377074e-5;];
    muE = 3.986e5 * (1e3)^3; %m^3/s^2


    rhat = r/norm(r);
    Vhat = V/norm(V);
    Vhat(isinf(Vhat)|isnan(Vhat)) = 0;
    
    g = gfunc(r);
    if t==0
        T=T*rhat;
    else
        T = T*Vhat;
    end
    drag = -0.5*rho * CD * A * (norm(V)).^2 /m * Vhat;
    dUdt(1:3) = V + cross(omegaE,r);
    dUdt(4:6) = T/m + drag/m + g;
end

function Vnew = vecrot(r,Vold,delta)
    k = cross(r,Vold)/norm(cross(r,Vold));
    Vnew = Vold*cos(delta) + (cross(k,Vold))*sin(delta) + k*dot(k,Vold)*(1-cos(delta));
end

function g =  gfunc(r)
    muE = 3.986e5 * (1e3)^3;
    g = - muE/norm(r).^3 .* r;
end


function gamma = gammafunc(r,V)
    CosTheta = max(min(dot(r,V)/(norm(r)*norm(V)),1),-1);
    gamma = real(acosd(CosTheta));
end

function [m,T,A,CD] = statefunc(t,m0,mdot0,tstage_index,tbo,T0,A0,V,h)
    Nstage = length(m0);
    
    M = norm(V)/atmos(h,13); 
    if M < 0.85
        CD = 0.2;
    else
        CD = 0.11+0.82/M^2-0.55/M^4;
    end
        % CD = 0.5;
    stage = 0;
    for i = 1:Nstage
        if i==Nstage
            stage = i;
            break;
        elseif t<tstage_index(i+1,2)
            stage = i;
            break;
        end
    end
    if t>tstage_index(end,3)
        m = m0(end) - (tbo(end))*mdot0(end);
        T = 0;
        A = A0(end);
    elseif stage==Nstage
        m = m0(stage) - mdot0(stage)*(t-tstage_index(stage,2));
        T = T0(stage);
        A = A0(stage);
    elseif t>tstage_index(stage,3) & t<tstage_index(stage+1,2)
        m = m0(stage+1);
        T = 0;
        A = A0(stage+1);
    else
        m = m0(stage) - mdot0(stage)*(t-tstage_index(stage,2));
        T = T0(stage);
        A = A0(stage);
    end
    
end

function r = latlong2cart(lat,long,h)
    RE = 6371e3;
    R = RE+h;
    r = R*[cos(lat)*cos(long); ...
            cos(lat)*sin(long); ...
            sin(lat)];
end

function [lat, long] = cart2latlong(r)
rmag = norm(r);
lat = asin(r(3)/rmag)*(180/pi);
    if (r(1) > 0) 
        long = atan(r(2)/r(1))*(180/pi);
    elseif (r(2) > 0) 
        long = atan(r(2)/r(1))*(180/pi) + 180;
     else 
        long = atan(r(2)/r(1))*(180/pi) - 180;
    end
end

function [value,isterminal,direction] = turncond(t,U,altPO)
    RE = 6371e3;

    r = U(1:3);

    h = norm(r) - RE;
    
    if h >= altPO
        value = 0;
        isterminal = 1;
        direction = 0;
    else
        value = 1;
        isterminal = 0;
        direction = 0;
    end
end


function [value,isterminal,direction] = crashcond(t,U)
    RE = 6371e3;

    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    if h < 0
        value = 0;
        isterminal = 1;
        direction = 0;
    else
        value = 1;
        isterminal = 0;
        direction = 0;
    end

end

function [value,isterminal,direction] = cruisecond(t,U,stopt)
    RE = 6371e3;

    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    if h < 0
        value = 0;
        isterminal = 1;
        direction = 0;
    % elseif h>500e3
    %     value = 0;
    %     isterminal = 1;
    %     direction = 0;
    elseif t >= 1*stopt
        value = 0;
        isterminal = 1;
        direction = 0;
    else
        value = 1;
        isterminal = 0;
        direction = 0;
    end

end

function [value,isterminal,direction] = steercond(t,U)
    RE = 6371e3;

    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    gamma = gammafunc(r,V);

    if gamma > 50
        value = 0;
        isterminal = 1;
        direction = 0;
    elseif h < 0
        value = 0;
        isterminal = 1;
        direction = 0;
    else
        value = 1;
        isterminal = 0;
        direction = 0;
    end

end