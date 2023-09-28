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
h0 = 2;
lat0 = 5.167713*d2r;
long0 = -52.683994*d2r;
r0 = latlong2cart(lat0, long0, h0);
rmag = norm(r0);
V0 = [0;0;0];

turn_azi = (90 - 26.6821)*d2r;

addrot = 1; % Set to 1 if the velocity from earth's rotation %
            % should be added 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROCKET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nstage = 2;
% m0     = [7.5384e4,7.9259e3];                % Initail/fueled mass kg
m0     = [6.9e4,8.502e3];                % Initail/fueled mass kg
% massfraction   = [0.1678,0.1873];        % mf/m0
massfraction   = [0.18,0.22];        % mf/m0
mf     = m0.*massfraction;       % Final/empty mass
mprop = m0-mf;
% T0      = [813e3, 85e3];     %0.76*       % Thrust N
T0      = [950e3, 50e3];     %0.76*       % Thrust N
% Isp    = [260,350];                 % Specific impulse s
Isp    = [308,343];                 % Specific impulse s
% d      = [3,3];                   % Diameter m
d      = [2,2];                   % Diameter m
tsep   =  1;
tstop = 450;  % Time when stage 2 should stop burning

altPO  = 100;
turn_fp = 89.7*d2r;
turnvec = 1*[cos(turn_fp)*cos(turn_azi); ...
        cos(turn_fp)*sin(turn_azi); ...
        sin(turn_fp)];
[X_turn, Y_turn, Z_turn] = enu2ecef(turnvec(1),turnvec(2),turnvec(3), ... 
    lat0,long0,h0,earth_REF,"radians");

TW = T0./(g0*m0)

A0     = pi*(d./2).^2; 
mdot0 = T0./(g0*Isp);
tbo = (m0-mf) ./ mdot0;  % Burn time for each stage
tbo_tot = sum(tbo); % Total burn time for the whole rocket
i_tbo1 = ceil(tbo(1));
i_tbo2 = i_tbo1 + ceil(tbo(2));
i_tbo = ceil(tbo_tot);

DeltaV = -g0*Isp.*log(massfraction);
DeltaVtot = sum(DeltaV);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE solving       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfin = 1*60*60;
U0 = [r0;V0;m0(1)];



opts_turn = odeset('RelTol',10e-10, 'Stats','on', ...
    'Events',@(t,U) turncond(t,U,altPO)); % Let ODE78 choose step size

opts_stage1 = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @(t,U) stagecond(t,U,mf,1,inf) );

opts_stage2 = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @(t,U) stagecond(t,U,mf,2,tstop) );

opts_cruise = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events',@(t,U) cruisecond(t,U));



% Runs pre-turn launch
tic % Used to check performance of solvers, can be ignored
[t_turn,U_turn] = ode45(@(t,U) ode_turn(t,U,mdot0,1,1,T0,A0), ...
    [0, 10*60], U0, opts_turn);

stage_index = ones(length(t_turn),1);
thrust_index = T0(1)*ones(length(t_turn),1);

% Normal gravity turn stage 1
timeatturn = t_turn(end);
Vatturn = norm(U_turn(end,4:6))
V_turn = norm(U_turn(end,4:6)) * ([X_turn; Y_turn; Z_turn] - r0);
[t_stage1,U_stage1] = ode45(@(t,U) ode_main(t,U,mdot0,1,1,T0,A0), ...
    [t_turn(end), tfin], [U_turn(end,1:3)';V_turn;U_turn(end,7)], opts_stage1);

stage_index = [stage_index;ones(length(t_stage1),1)];
thrust_index = [thrust_index;T0(1)*ones(length(t_stage1),1)];

if norm(U_stage1(end,1:3))<= RE+1
    t_stage2_burn1 = [];
    U_stage2_burn1 = [];
else
    t_stagesep = t_stage1(end);
    [t_stagesep_res,U_stagesep] = ode45(@(t,U) ode_main(t,U,mdot0,2,0,T0,A0), ...
    [t_stagesep, t_stagesep+tsep], [U_stage1(end,1:3)';U_stage1(end,4:6)';m0(2)], opts_stage2);
    stage_index = [stage_index;0*ones(length(t_stagesep_res),1)];
    thrust_index = [thrust_index;0*ones(length(t_stagesep_res),1)];

    t_stage2_start = t_stagesep_res(end);
    
    if norm(U_stagesep(end,1:3))>= RE+1
        [t_stage2_burn1,U_stage2_burn1] = ode45(@(t,U) ode_main(t,U,mdot0,2,1,T0,A0), ...
        [t_stage2_start, inf], [U_stagesep(end,1:3)';U_stagesep(end,4:6)';U_stagesep(end,7)], opts_stage2);
         stage_index = [stage_index;2*ones(length(t_stage2_burn1),1)];
        thrust_index = [thrust_index; T0(2)*ones(length(t_stage2_burn1),1)];
        if norm(U_stagesep(end,1:3)) >= RE+1
            [t_stage2_cruise,U_stage2_cruise] = ode45(@(t,U) ode_main(t,U,mdot0,2,0,T0,A0), ...
            [t_stage2_burn1(end), inf], [U_stage2_burn1(end,1:3)';U_stage2_burn1(end,4:6)'+addrot*cross(omegaE,r0);U_stage2_burn1(end,7)], opts_cruise);
            stage_index = [stage_index;2*ones(length(t_stage2_cruise),1)];
            thrust_index = [thrust_index;0*ones(length(t_stage2_cruise),1)];

            if norm(U_stage2_cruise(end,1:3)) >= RE+1
                 r_trgt = norm(U_stage2_cruise(end,1:3)');
                 V_trgt = sqrt(muE/r_trgt);
                opts_circularize = odeset('RelTol',1e-10, 'MaxStep',1 , ...
                'Stats','on', 'Events',@(t,U) circularizecond(t,U,V_trgt,mf(2)));
                [t_stage2_burn2,U_stage2burn2] = ode45(@(t,U) ode_main(t,U,mdot0,2,1,T0,A0), ...
                [t_stage2_cruise(end), inf], [U_stage2_cruise(end,1:3)'; ...
                U_stage2_cruise(end,4:6)';U_stage2_cruise(end,7)], opts_circularize);

                stage_index = [stage_index;2*ones(length(t_stage2_burn2),1)];
                thrust_index = [thrust_index;T0(2)*ones(length(t_stage2_burn2),1)];
                
                if norm(U_stage2burn2(end,1:3)) > RE + 1
                 opts_stage2.Events = @(t,U) stagecond(t,U,mf,2,inf);
                 [t_stage2_orbit,U_stage2_orbit] = ode45(@(t,U) ode_main(t,U,mdot0,2,0,T0,A0), ...
                [t_stage2_burn2(end), t_stage2_burn2(end)+180*60], [U_stage2burn2(end,1:3)';U_stage2burn2(end,4:6)';U_stage2burn2(end,7)], opts_stage2);
                stage_index = [stage_index;2*ones(length(t_stage2_orbit),1)];
                thrust_index = [thrust_index;0*ones(length(t_stage2_orbit),1)];
                else
                    t_stage2_orbit = [];
                    U_stage2_orbit = [];
                end
            else
                t_stage2_burn2 = [];
                U_stage2burn2 = [];
            end
        else 
            t_stage2_cruise = [];
            U_stage2_cruise = [];
        end
    else
        t_stage2_burn1 = [];
        U_stage2_burn1 = [];
    end
    
   
    
   
end

% Controlled turn

t_main = [t_turn; t_stage1; t_stagesep_res; t_stage2_burn1;t_stage2_cruise;t_stage2_burn2;t_stage2_orbit];
U_main = [U_turn; U_stage1;U_stagesep;U_stage2_burn1;U_stage2_cruise;U_stage2burn2;U_stage2_orbit];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessing      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(t_main);

% Date used to convert between ECI and ECEF
tt = [2000 1 1 17 15 10];
epoch = datetime(2000,1,1,17,15,10);

t_stage1_bo = t_stage1(end);
t_stage2_cutoff = t_stage2_burn1(end);
t_stage2_relight = t_stage2_burn2(1);
t_stage2_bo = t_stage2_burn2(end);
t_events = [t_stage1_bo t_stage2_cutoff t_stage2_relight t_stage2_bo];

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
mres = U_main(:,7);   % Mass as function of time

Mres = zeros(N,1);   % Mach 
CDres = zeros(N,1);  % Drag coeff
rhot = zeros(N,1);
gres = zeros(N,1);
A = zeros(N,1);
Tres = thrust_index;

for j = 1:N
    % Lat-Long in ECI
    [latlong_res(j,1), latlong_res(j,2)] = cart2latlong(r_res(j,:));
    
    gamma(j) = gammafunc(r_res(j,:),V_res(j,:)); 
    
    CDres(j) = CD_func(r_res(j,:),V_res(j,:));
    
    Mres(j) = norm(V_res(j,:))/atmos(h_res(j),13);

    rhot(j) = atmos(h_res(j),12);

    gres(j) = -norm(gfunc(r_res(j,:)));
    j
    if stage_index(j) == 1
        A(j) = A0(1);
    else
        A(j) = A0(2);
    end

end

q_R = 0.5*Vmag_res.^2.*rhot;

% Time span when deltaV losses should be calculated
t_burn = [t_turn; t_stage1; t_stagesep_res; t_stage2_burn1;t_stage2_cruise;t_stage2_burn2];
i_burn = length(t_burn);

% Delta V calculations 
air_loss_int = CDres.*A.*q_R./mres;
g_loss_int = gres.*sin(gamma*d2r);
delta_V_air = -cumtrapz(t_burn,air_loss_int(1:i_burn)); % Cumulative integral
delta_V_air_tot = -trapz(t_burn,air_loss_int(1:i_burn));

delta_V_grav = -cumtrapz(t_burn,g_loss_int(1:i_burn));
delta_V_grav_tot = -trapz(t_burn,g_loss_int(1:i_burn));

% Losses are only accounted for during burn time

% Uses less points for plotting due to performance
nf = 100;                   % Number of frames
sf = floor(N/nf);           % Skip-factor

% Parameters in ECEF
r_ecef = zeros(nf,3);
v_ecef = zeros(nf,3);
latlong_ecef = zeros(nf, 2);
tsf = zeros(nf,1);


for k = 1:nf
    tsf(k) = t_main(sf*k);
    
    % UTC at time t
    utc = epoch + seconds(t_main(k*sf));
    tt = datevec(utc);      % Converts format to 1x6 vector

    lla = eci2lla(r_res(sf*k,:),tt);
    latlong_ecef(k,1) = lla(1);
    latlong_ecef(k,2) = lla(2);
    disp(['ECEF-calc step: ', num2str(k), ' of ', num2str(nf)])
end
%%
mp_curr = U_main(end,7);

DeltaVman1 = 194.2;
DeltaVman2 = 190;

man1 = mp_curr * ( 1 - exp( DeltaVman1 / ( Isp(2) * g0 ) )  )
man2 = (mp_curr+man1) * ( 1 - exp( DeltaVman2 / ( Isp(2) * g0 ) )  )

sparefuel = (mp_curr+man1+man2)-mf(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
subplot(6,1,1)
plot(t_main,Vmag_res/1e3)
ylabel('$|V|$ [km/s]')
xline(t_events)

subplot(6,1,2)
plot(t_main,h_res/1e3)
ylabel('$|h|$ [km]')
xline(t_events)

subplot(6,1,3)
plot(t_main,mres/1e3)
ylabel('$m$ [ton]') 
xline(t_events)

subplot(6,1,4)
plot(t_main,gamma)
ylabel('$\gamma$ [$^o$]')
xline(t_events)

subplot(6,1,5)
plot(t_main, Tres/1e6)
ylabel('$T$ [MN]')
xline(t_events)
subplot(6,1,6)
plot(t_main,a_res/g0)
xlabel('$t$, [s]')
ylabel('$a/g_0$ [-]')
ylim([-1,10])

figure(4)
subplot(2,1,1)
plot(t_main, latlong_res(:,1))
hold on
plot(tsf, latlong_ecef(:,1),'--')
hold off
ylabel('lat [$^o$]')
legend('ECI','ECEF')
subplot(2,1,2)
plot(t_main, latlong_res(:,2))
hold on
plot(tsf, latlong_ecef(:,2),'--')
hold off
xlabel('$t$, [s]')
ylabel('long [$^o$]')

ttemp = (find(t_main<1.1*t_stage2_burn2(end)));
ttemp = t_main(1:ttemp(end));
TT = length(ttemp);
figure(6)
subplot(6,1,1)
plot(ttemp,Vmag_res(1:TT)/1e3)
ylabel('$|V|$ [km/s]')
xline(t_events)

subplot(6,1,2)
plot(ttemp,h_res(1:TT)/1e3)
ylabel('$|h|$ [km]')
xline(t_events)

subplot(6,1,3)
plot(ttemp,mres(1:TT)/1e3)
ylabel('$m$ [ton]') 
xline(t_events)

subplot(6,1,4)
plot(ttemp,gamma(1:TT))
ylabel('$\gamma$ [$^o$]')
xline(t_events)
subplot(6,1,5)
plot(ttemp, Tres(1:TT)/1e6)
ylabel('$T$ [MN]')
xline(t_events)
subplot(6,1,6)
plot(ttemp,a_res(1:TT)/g0)
xlabel('$t$, [s]')
ylabel('$a/g_0$ [-]')
ylim([-1,10])
xline(t_events)

figure(7)
plot(ttemp,q_R(1:TT)/1e3)
xlabel('$t$, [s]')
ylabel('Dynamic pressure $q$ [kPa]')


figure(200)
geoscatter(latlong_ecef(:,1), latlong_ecef(:,2) )

figure(201)
subplot(2,1,1)
plot(t_main, CDres)
ylabel('$C_D$ [-]')
xlim([0,600])
subplot(2,1,2)
plot(t_main, Mres)
ylabel('$M$ [-]')
xlim([0,600])

figure(202)
plot(t_burn, delta_V_grav   )
xlabel('$t$ [s]')
ylabel('Cumulative $\Delta V_{\mathrm{grav}}$ [km/s]')

f1 = figure(1);
ax1 = gca;
ax1.GridLineWidthMode = "auto";
grid on
ax1.Projection = "perspective";
ax1.PlotBoxAspectRatioMode = "manual";
ax1.PlotBoxAspectRatio = [1 1 1]; 
ax1.XLim = 1.01*maxr*[-1, 1];
ax1.YLim = 1.01*maxr*[-1, 1];
ax1.ZLim = 1.01*maxr*[-1, 1];
view([30,30])

hold on
[e1, axes1] = earth_test(t_main(end),RE);
r_rocket = plot3(U_main(:,1),U_main(:,2),U_main(:,3), 'LineWidth',3,'Color','red');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I screwd up something here, it shouldn't become more inefficient with
% time, now it goes fast in the beginning then slows down alot
f2 = figure(2);
close(f2)
f2 = figure(2);
ax2 = gca;
ax2.GridLineWidthMode = "auto";
grid on
ax2.Projection = "perspective";
ax2.PlotBoxAspectRatioMode = "manual";
ax2.PlotBoxAspectRatio = [1 1 1]; 
ax2.XLim = 1.01*maxr*[-1, 1];
ax2.YLim = 1.01*maxr*[-1, 1];
ax2.ZLim = 1.01*maxr*[-1, 1];
view([30,30])

hold on
[e2, axes] = earth_test(0,RE);
r_rocket = plot3(U_main(1,1),U_main(1,2),U_main(1,3), 'LineWidth',3,'Color','red');
hold off

for i = 1:nf
    axes.reset;   
    e2.reset;

    [e2, axes] = earth_test(t_main(sf*i),RE);
    % r_rocket = plot3(U_main(1:i,1),U_main(1:i,2),U_main(1:i,3), 'LineWidth',3,'Color','red');
    r_rocket.XData =  U_main(1:sf*i,1);
    r_rocket.YData = U_main(1:sf*i,2);
    r_rocket.ZData = U_main(1:sf*i,3);
    
    drawnow
    disp(['Movie frame: ', num2str(i), ' of ', num2str(nf)])
end
toc


%% Functions 
function dUdt = ode_turn(t,U,mdot0,stage,thrust_bool,T0,A0)
    RE = 6371e3;

    dUdt = zeros(7,1);
    m = U(7);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    rho = atmos(h,12); 
    
    A = A0(stage);
    mdot = mdot0(stage);

    CD = CD_func(r,V);

    if thrust_bool==0
        T = 0;
        dUdt(7) = 0;
    else
        T = T0(stage);
        dUdt(7) = -mdot;
    end

   
    rhat = r/norm(r);
    Vhat = V/norm(V);
    Vhat(isinf(Vhat)|isnan(Vhat)) = 0;

    g = gfunc(r);
    if t==0
        T=T*rhat;
    else
        T = T*Vhat;
    end
    drag = -0.5*rho * CD * A * norm(V).^2 /m * Vhat;
    dUdt(1:3) = V ;
    dUdt(4:6) = T/m + drag/m + g;

end

function dUdt = ode_main(t,U,mdot0,stage,thrust_bool,T0,A0)
    RE = 6371e3;

    dUdt = zeros(7,1);
    m = U(7);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    rho = atmos(h,12); 
    
    A = A0(stage);
    mdot = mdot0(stage);
    CD = CD_func(r,V);
   if thrust_bool==0
        T = 0;
        dUdt(7) = 0;
    else
        T = T0(stage);
        dUdt(7) = -mdot;
    end
    
    omegaE = 0*[0;0;7.292115855377074e-5;];



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

function dUdt = ode_Toff(t,U,m0,mdot0,tstage_index,tbo,T0,A0,m)
    RE = 6371e3;

    dUdt = zeros(4,1);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    rho = atmos(h,12); 
    
    [~,~,A,CD] = statefunc(t,m0,mdot0,tstage_index,tbo,T0,A0,V,h);
    
    T = 0;

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


function dUdt = ode_steer(t,U,mdot0,stage,thrust_bool,T0,A0, gamma0,t0)
    RE = 6371e3;   
    d2r = pi/180;
    
    dUdt = zeros(7,1);
    
    r = U(1:3);
    Vold = U(4:6);
    m = U(7);
    h = norm(r) - RE;
    gammaold = gammafunc(r,Vold)*d2r;
    gammaend = 0.8*pi/2;
    gamma =  gamma0*d2r + (t-t0)*( (gammaend-gamma0*d2r)/(tstage_index(end,3)-t0)  );
    

    delta = gamma - gammaold;
    
    V = vecrot(r,Vold,delta);
    
    
    rho = atmos(h,12); 
    
    A = A0(stage);
    mdot = mdot0(stage);
    CD = CD_func(r,V);
    if thrust_bool==0
        T = 0;
        dUdt(7) = 0;
    else
        T = T0(stage);
        dUdt(7) = -mdot;
    end

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
    J2 = 1.08262668e-3;

    x = r(1);
    y = r(2);
    z = r(3);

    R = norm(r);
    k_g = -muE/R^3;
    k_xy = k_g*(1+1.5*J2*(1-5*(z/R)^2)/R^2 );
    k_prim = k_g*(1+1.5*J2*(3-5*(z/R)^2)/R^2 );

    g = [k_xy*x; k_xy*y; k_prim*z];
end


function gamma = gammafunc(r,V)
    CosTheta = max(min(dot(r,V)/(norm(r)*norm(V)),1),-1);
    gamma = 90 - real(acosd(CosTheta));
end

function CD = CD_func(r,V)
    RE = 6371e3;
    h = norm(r) - RE; 
    M = norm(V)/atmos(h,13); 
     if M < 0.85
        CD = 0.2;
    else
        CD = 0.23+0.82/M^2-0.55/M^4;
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

function [value,isterminal,direction] = stagecond(t,U,mf0,stage,tstop)
    RE = 6371e3;
    m = U(7);
    r = U(1:3);
    h = norm(r) - RE;
    if m <= mf0(stage)
        value = 0;
        isterminal = 1;
        direction = 0;
    elseif h < 0
        value = 0;
        isterminal = 1;
        direction = 0;
    elseif t>=tstop
        value = 0;
        isterminal = 1;
        direction = 0;
    else
        value = 1;
        isterminal = 0;
        direction = 0;
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

function [value,isterminal,direction] = cruisecond(t,U)
    RE = 6371e3;
    d2r = pi/180;
    r = U(1:3);
    V = U(4:6);
    h = norm(r)-RE;
    gamma = gammafunc(r,V);
    if abs(gamma-2.7) < 0.1
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

function  [value,isterminal,direction] = circularizecond(t,U,V_trgt,mf)
    RE = 6371e3;
    r = U(1:3);
    V = U(4:6);
    h = norm(r)-RE;
    m = U(7);

    if norm(V) >= 1.001*V_trgt
        value = 0;
        isterminal = 1;
        direction = 0;
        disp('Circularizecond target velocity achived')
    elseif m <= mf
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

function [value,isterminal,direction] = steercond(t,U)
    RE = 6371e3;

    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    gamma = gammafunc(r,V);

    if gamma > 70
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