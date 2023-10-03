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
p0 = 1013.25e2;
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
m0     = [1.0823e+05,2.3262e+04];                 % Initial/fueled mass kg
massfraction   = [0.35,0.0708];        % mf/m0
mf     = m0.*massfraction;       % Final/empty mass
mprop  = m0-mf
ms(1)     = mf(1) - m0(2);
ms(2)     = mf(2) - 1000;
T0      = [760e3 + 3*760e3/4, 2*760e3/4];     %0.76*       % Thrust N
Isp    = [308,308];                  % Specific impulse s
d      = [3,3];                   % Diameter m
Ae0    = pi*[2*25e-2/2, 25e-2/2].^2; % Nozzle area
p_e = 0.95*p0;                        % Exhuast pressure                 % Diameter m
tsep   =  1;
tstop  = 200;  % Time when stage 2 should stop burning

stage2_cont = 1;

altPO  = 280;
turn_fp = 89.85*d2r;
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

DeltaV = -g0*Isp.*log(massfraction)
DeltaVtot = sum(DeltaV)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE solving       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfin = 1*60*60;
U0 = [r0;V0;m0(1)];

hstopfreefall = 20e3;
Vcrit = 800;
delta_t = 10;
V_break = 5;
h_break = 2.08e3;

opts_turn = odeset('RelTol',10e-10, 'Stats','on', ...
    'Events',@(t,U) turncond(t,U,altPO)); % Let ODE78 choose step size

opts_stage1 = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @(t,U) stage1cond(t,U,m0,ms,Isp,g0,mprop) ); %gives initial values

opts_reentb = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @(t,U) reentbcond(t,U,mf,m0,ms,Vcrit) );


opts_freefall = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @(t,U) freefallcond(t,U,hstopfreefall) );

opts_break2 = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @(t,U) break2cond(t,U,V_break,h_break) );

opts_land = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @(t,U) landcond(t,U,mf,ms,m0) );


% LAUNCH AND FLIGTH UP
tic % Used to check performance of solvers, can be ignored
%1 Turn
[t_turn,U_turn] = ode45(@(t,U) ode_turn(t,U,Isp,1,1,T0,A0,Ae0,p_e,1), ...
    [0, 10*60], U0, opts_turn);

stage_index = ones(length(t_turn),1);
thrust_index = T0(1)*ones(length(t_turn),1);

%2 Stage One
V_turn = norm(U_turn(end,4:6)) * ([X_turn; Y_turn; Z_turn] - r0);
V_at_turn = norm(U_turn(end,4:6))
[t_stage1,U_stage1] = ode45(@(t,U) ode_main(t,U,Isp,1,1,T0,A0,Ae0,p_e,1), ...
    [t_turn(end), tfin], [U_turn(end,1:3)';V_turn;U_turn(end,7)], opts_stage1);

% plot(t_stage1,vecnorm(U_stage1(:,1:3),2,2)-RE)
% gamma_S1 = gammafunc(U_stage1(:,1:3),U_stage1(:,4:6));
stage_index = [stage_index;ones(length(t_stage1),1)];
thrust_index = [thrust_index;T0(1)*ones(length(t_stage1),1)];



%SEPARATION
mass_return = U_stage1(end,7)-m0(2);

opts_delay = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @(t,U) delay_cond(t,U,t_stage1(end),delta_t) );

% Delay before slowdown burn
[t_delay,U_delay] = ode78(@(t,U) ode_main(t,U,Isp,1,0,T0,A0,Ae0,p_e,2), ...
    [t_stage1(end), inf], [U_stage1(end,1:3)';U_stage1(end,4:6)';mass_return], opts_delay);
stage_index = [stage_index;1*ones(length(t_delay),1)];
thrust_index = [thrust_index;0*ones(length(t_delay),1)];


%REENTRY STARTS
%3 Reentry burn
[t_reentb,U_reentb] = ode78(@(t,U) ode_reentb(t,U,Isp,1,1,T0,A0,Ae0,p_e,1), ...
    [t_delay(end), inf], [U_delay(end,1:3)';U_delay(end,4:6)';U_delay(end,7)], opts_reentb);

stage_index = [stage_index;1*ones(length(t_reentb),1)];
thrust_index = [thrust_index;0.8*T0(1)*ones(length(t_reentb),1)];
plot(t_reentb,vecnorm(U_reentb(:,1:3),2,2)-RE)
%4 Free fall
[t_freefall,U_freefall] = ode78(@(t,U) ode_main(t,U,Isp,1,0,T0,A0,Ae0,p_e,2), ...
    [t_reentb(end), inf], ...
    [U_reentb(end,1:3)';U_reentb(end,4:6)';U_reentb(end,7)], opts_freefall);

stage_index = [stage_index;1*ones(length(t_freefall),1)  ];
thrust_index = [thrust_index;0*ones(length(t_freefall),1)];

%5 Slowdown/brakeburn2
[t_bb2,U_bb2] = ode45(@(t,U) ode_land(t,U,Isp,ms,1,1,T0,A0,Ae0,p_e,0.1), ...
    [t_freefall(end), inf], ...
    [U_freefall(end,1:3)';U_freefall(end,4:6)';U_freefall(end,7)], opts_break2);

stage_index = [stage_index;1*ones(length(t_bb2),1)  ];
thrust_index = [thrust_index;0.5*T0(1)*ones(length(t_bb2),1)];

%6 Landing
[t_land,U_land] = ode45(@(t,U) ode_land(t,U,Isp,ms,1,1,T0,A0,Ae0,p_e,0.5), ...
    [t_bb2(end), inf], [U_bb2(end,1:3)';U_bb2(end,4:6)';U_bb2(end,7)], opts_land);

stage_index = [stage_index;1*ones(length(t_land),1)];
thrust_index = [thrust_index;T0(1)*1*ones(length(t_land),1)];

%TOTAL FLIGHT Stage 1
t_main = [t_turn; t_stage1; t_delay; t_reentb; t_freefall; t_bb2; t_land];
U_main = [U_turn; U_stage1; U_delay; U_reentb; U_freefall; U_bb2; U_land];

%% Stage 2 calc
opts_stage2 = odeset('RelTol',1e-10, 'MaxStep',1 , ...
    'Stats','on', 'Events', @(t,U) stagecond(t,U,mf,2,tstop) );
opts_cruise = odeset('RelTol',1e-10, 'MaxStep',0.5 , ...
    'Stats','on', 'Events',@(t,U) cruisecond(t,U));

if norm(U_stage1(end,1:3))<= RE+1
    t_stage2_burn1 = [];
    U_stage2_burn1 = [];
else
    t_stagesep = t_stage1(end);
    [t_stagesep_res,U_stagesep] = ode45(@(t,U) ode_main(t,U,Isp,2,0,T0,A0,Ae0,p_e,1), ...
    [t_stagesep, t_stagesep+tsep], [U_stage1(end,1:3)';U_stage1(end,4:6)';m0(2)], opts_stage2);
    stage_index2 = [stage_index;0*ones(length(t_stagesep_res),1)];
    thrust_index2 = [thrust_index;0*ones(length(t_stagesep_res),1)];

    t_stage2_start = t_stagesep_res(end);
    
    if norm(U_stagesep(end,1:3))>= RE+1
        [t_stage2_burn1,U_stage2_burn1] = ode45(@(t,U) ode_main(t,U,Isp,2,1,T0,A0,Ae0,p_e,1), ...
        [t_stage2_start, inf], [U_stagesep(end,1:3)';U_stagesep(end,4:6)';U_stagesep(end,7)], opts_stage2);
         stage_index2 = [stage_index2;2*ones(length(t_stage2_burn1),1)];
        thrust_index2 = [thrust_index2; ones(length(t_stage2_burn1),1)];
        if norm(U_stagesep(end,1:3)) >= RE+1
            [t_stage2_cruise,U_stage2_cruise] = ode45(@(t,U) ode_main(t,U,Isp,2,0,T0,A0,Ae0,p_e,1), ...
            [t_stage2_burn1(end), inf], [U_stage2_burn1(end,1:3)';U_stage2_burn1(end,4:6)'+addrot*cross(omegaE,r0);U_stage2_burn1(end,7)], opts_cruise);
            stage_index2 = [stage_index2;2*ones(length(t_stage2_cruise),1)];
            thrust_index2 = [thrust_index2;0*ones(length(t_stage2_cruise),1)];
            
            figure(1)
            plot(t_stage2_cruise,vecnorm(U_stage2_cruise(:,1:3),2,2)/1e3-RE/1e3)

            figure(2)
            gamma_temp = gammafunc(U_stage2_cruise(:,1:3),U_stage2_cruise(:,4:6));
            plot(t_stage2_cruise,gamma_temp)

            figure(3)
            plot(t_stage2_cruise,vecnorm(U_stage2_cruise(:,4:6),2,2)/1e3)
            if norm(U_stage2_cruise(end,1:3)) >= RE+1
                 r_trgt = norm(U_stage2_cruise(end,1:3)');
                 V_trgt = sqrt(muE/r_trgt);
                opts_circularize = odeset('RelTol',1e-10, 'MaxStep',1 , ...
                'Stats','on', 'Events',@(t,U) circularizecond(t,U,V_trgt,mf(2)));
                [t_stage2_burn2,U_stage2burn2] = ode45(@(t,U) ode_main(t,U,Isp,2,1,1*T0,A0,Ae0,p_e,1), ...
                [t_stage2_cruise(end), inf], [U_stage2_cruise(end,1:3)'; ...
                U_stage2_cruise(end,4:6)';U_stage2_cruise(end,7)], opts_circularize);

                stage_index2 = [stage_index2;2*ones(length(t_stage2_burn2),1)];
                thrust_index2 = [thrust_index2;ones(length(t_stage2_burn2),1)];
                
                if norm(U_stage2burn2(end,1:3)) > RE + 1
                 opts_stage2.Events = @(t,U) stagecond(t,U,mf,2,inf);
                 [t_stage2_orbit,U_stage2_orbit] = ode45(@(t,U) ode_main(t,U,Isp,2,0,T0,A0,Ae0,p_e,1), ...
                [t_stage2_burn2(end), t_stage2_burn2(end)+180*60], [U_stage2burn2(end,1:3)';U_stage2burn2(end,4:6)';U_stage2burn2(end,7)], opts_stage2);
                stage_index = [stage_index;2*ones(length(t_stage2_orbit),1)];
                thrust_index2 = [thrust_index2;0*ones(length(t_stage2_orbit),1)];
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
s2_t = [t_turn; t_stage1; t_stagesep_res; t_stage2_burn1;t_stage2_cruise;t_stage2_burn2;t_stage2_orbit];
s2_U = [U_turn; U_stage1;U_stagesep;U_stage2_burn1;U_stage2_cruise;U_stage2burn2;U_stage2_orbit];
s2_Tindex = thrust_index2;
s2_sindex = stage_index2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessing      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(t_main);

% Date used to convert between ECI and ECEF
tt = [2000 1 1 17 15 10];
epoch = datetime(2000,1,1,17,15,10);

t_stage1_bo = t_stage1(end);

t_events = [t_stage1_bo];

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
    
    CDres(j) = CD_func2(r_res(j,:),V_res(j,:));
    
    Mres(j) = norm(V_res(j,:))/atmos(h_res(j),13);

    rhot(j) = atmos(h_res(j),12);

    gres(j) = -norm(gfunc(r_res(j,:)));
    
    Qdot(j) = 1/2*rhot(j)*norm(V_res(j,:))^3*A0(1)*CDres(j,:)/20;


    if stage_index(j) == 1
        A(j) = A0(1);
    else
        A(j) = A0(2);
    end

end

q_R = 0.5*Vmag_res.^2.*rhot;



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
mp_curr = s2_U(end,7);

oev = eci2orb_gooding (muE/(1e3)^3, s2_U(end,1:3)/1e3, s2_U(end,4:6)/1e3)
rA = oev(1)*(1+oev(2))*1e3;
rP = oev(1)*(1-oev(2))*1e3;

DeltaVman1 = 194.2;
DeltaVman2 = 190;

man1 = mp_curr * ( 1 - exp( DeltaVman1 / ( Isp(2) * g0 ) )  );
man2 = (mp_curr+man1) * ( 1 - exp( DeltaVman2 / ( Isp(2) * g0 ) )  );

sparefuel = (mp_curr+man1+man2)-mf(2);

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

ttemp = (find(t_main<500));
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
plot(t_main,q_R/1e3)
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
%%
gammas2 = gammafunc(s2_U(:,1:3), s2_U(:,4:6));
figure(789)
subplot(5,1,1)

plot(s2_t,vecnorm(s2_U(:,4:6),2,2)/1e3)
ylabel('$|V|$ [km/s]')
% xline(t_events)

subplot(5,1,2)
plot(s2_t,vecnorm(s2_U(:,1:3),2,2)/1e3-RE/1e3)
ylabel('$|h|$ [km]')
% xline(t_events)

subplot(5,1,3)
plot(s2_t,s2_U(:,7)/1e3)
ylabel('$m$ [ton]') 
% xline(t_events)

subplot(5,1,4)
plot(s2_t,gammas2)
ylabel('$\gamma$ [$^o$]')
% xline(t_events)

subplot(5,1,5)
plot(s2_t, s2_Tindex/1e6)
ylabel('$T$ [MN]')
% xline(t_events)

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

view(ax1,[89.5611777570063 -6.0387365921562]);
xlim(ax1,[3585720.96991818 4014642.20187433]);

ylim(ax1,[-5190290.94498331 -4761369.71302716]);

zlim(ax1,[407538.30778831 836459.539744459]);
hold on
[e1, axes1] = earth_test(t_main(end),RE);
r_rocket = plot3(U_main(:,1),U_main(:,2),U_main(:,3), 'LineWidth',3,'Color','red');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % I screwd up something here, it shouldn't become more inefficient with
% % time, now it goes fast in the beginning then slows down alot
% f2 = figure(2);
% close(f2)
% f2 = figure(2);
% ax2 = gca;
% ax2.GridLineWidthMode = "auto";
% grid on
% ax2.Projection = "perspective";
% ax2.PlotBoxAspectRatioMode = "manual";
% ax2.PlotBoxAspectRatio = [1 1 1]; 
% ax2.XLim = 1.01*maxr*[-1, 1];
% ax2.YLim = 1.01*maxr*[-1, 1];
% ax2.ZLim = 1.01*maxr*[-1, 1];
% view([30,30])
% 
% hold on
% [e2, axes] = earth_test(0,RE);
% r_rocket = plot3(U_main(1,1),U_main(1,2),U_main(1,3), 'LineWidth',3,'Color','red');
% hold off
% 
% for i = 1:nf
%     axes.reset;   
%     e2.reset;
% 
%     [e2, axes] = earth_test(t_main(sf*i),RE);
%     % r_rocket = plot3(U_main(1:i,1),U_main(1:i,2),U_main(1:i,3), 'LineWidth',3,'Color','red');
%     r_rocket.XData =  U_main(1:sf*i,1);
%     r_rocket.YData = U_main(1:sf*i,2);
%     r_rocket.ZData = U_main(1:sf*i,3);
% 
%     drawnow
%     disp(['Movie frame: ', num2str(i), ' of ', num2str(nf)])
% end
% toc
% 
% 
%% Functions 
function dUdt = ode_turn(t,U,Isp,stage,thrust_bool,T0,A0,Ae0,p_e,Thrust_prop)
    RE = 6371e3;

    dUdt = zeros(7,1);
    m = U(7);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    rho = atmos(h,12); 
    
    A = A0(stage);
    g0 = 9.80665;

    CD = CD_func(r,V);

    if thrust_bool==0
        T = 0;
        dUdt(7) = 0;
    else
        T = Tfunc(h,T0,Ae0,p_e,stage)*Thrust_prop;
        dUdt(7) = -T./(g0*Isp(stage));
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
    drag = -0.5*rho * CD * A * norm(V).^2  * Vhat;
    dUdt(1:3) = V ;
    dUdt(4:6) = T/m + drag/m + g;

end

function dUdt = ode_reentb(t,U,Isp,stage,thrust_bool,T0,A0,Ae0,p_e,Thust_prop)
    RE = 6371e3;

    dUdt = zeros(7,1);
    m = U(7);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    g0 = 9.80665;
    rho = atmos(h,12); 
    
    A = A0(stage);
    
    CD = CD_func2(r,V);
   if thrust_bool==0
        T = 0;
        dUdt(7) = 0;
    else
         T = Tfunc(h,T0,Ae0,p_e,stage)*Thust_prop;
        dUdt(7) = -T./(g0*Isp(stage));
    end
    
    omegaE = 0*[0;0;7.292115855377074e-5;];



    rhat = r/norm(r);
    Vhat = V/norm(V);
    Vhat(isinf(Vhat)|isnan(Vhat)) = 0;

    g = gfunc(r);
    if t==0
        T=T*rhat;
    else
        T = T*-Vhat;
    end
    drag = -0.5*rho * CD * A * (norm(V)).^2  * Vhat;
    dUdt(1:3) = V + cross(omegaE,r);
    dUdt(4:6) = T/m + drag/m + g;
end

function dUdt = ode_land(t,U,Isp,ms,stage,thrust_bool,T0,A0,Ae0,p_e,Thrust_prop)
    RE = 6371e3;
    g0 = 9.80665;
    dUdt = zeros(7,1);
    m = U(7);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    rho = atmos(h,12); 
    
    A = A0(stage);
    
    CD = CD_func2(r,V);

    norm(V)

   if thrust_bool==0
        T = 0;
        dUdt(7) = 0;
   elseif m > ms(stage)
        T = Tfunc(h,T0,Ae0,p_e,stage)*Thrust_prop;
        dUdt(7) = -T./(g0*Isp(stage));
   else
       T = 0;
       dUdt(7) = 0;
    end
    
    omegaE = 0*[0;0;7.292115855377074e-5;];

    rhat = r/norm(r);
    Vhat = V/norm(V);
    Vhat(isinf(Vhat)|isnan(Vhat)) = 0;

    g = gfunc(r);
    if t==0
        T=T*rhat;
    else
        T = T*-Vhat;
    end
    drag = -0.5*rho * CD * A * (norm(V)).^2  * Vhat;
    dUdt(1:3) = V + cross(omegaE,r);
    dUdt(4:6) = T/m + drag/m + g;
end

function dUdt = ode_main(t,U,Isp,stage,thrust_bool,T0,A0,Ae0,p_e,CDmodel)
    RE = 6371e3;
    g0 = 9.80665;
    dUdt = zeros(7,1);
    m = U(7);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    rho = atmos(h,12); 
    
    A = A0(stage);
    if CDmodel == 1
        CD = CD_func(r,V);
    else
        CD = CD_func2(r,V);
    end

   if thrust_bool==0
        T = 0;
        dUdt(7) = 0;
    else
        T = Tfunc(h,T0,Ae0,p_e,stage);
        dUdt(7) = -T./(g0*Isp(stage));
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
    drag = -0.5*rho * CD * A * (norm(V)).^2  * Vhat;
    dUdt(1:3) = V + cross(omegaE,r);
    dUdt(4:6) = T/m + drag/m + g;
end

function dUdt = ode_Toff(t,U,Isp,stage,T0,A0,Ae0,p_e)
    RE = 6371e3;
    g0 = 9.80665;
    dUdt = zeros(7,1);
    m = U(7);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    rho = atmos(h,12); 
    
    A = A0(stage);
    
    CD = CD_func(r,V);
    T = 0; 
    dUdt(7) = 0;
    
    omegaE = 0*[0;0;7.292115855377074e-5;];

    rhat = r/norm(r);
    Vhat = V/norm(V);
    Vhat(isinf(Vhat)|isnan(Vhat)) = 0;

    g = gfunc(r);
   
    drag = -0.5*rho * CD * A * (norm(V)).^2 * Vhat;
    dUdt(1:3) = V + cross(omegaE,r);
    dUdt(4:6) = T/m + drag/m + g;
end


function dUdt = ode_steer(t,U,Isp,stage,thrust_bool,T0,A0,Ae0,p_e, gamma0,t0)
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
    
    CD = CD_func(r,V);
    if thrust_bool==0
        T = 0;
        dUdt(7) = 0;
    else
        T = Tfunc(h,T0,Ae0,p_e,stage);
        dUdt(7) = -T./(g0*Isp(stage));
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
    drag = -0.5*rho * CD * A * (norm(V)).^2  * Vhat;
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
N = size(r,1);
gamma = zeros(N,1);
for i = 1:N
    CosTheta = max(min(dot(r(i,:),V(i,:))/(norm(r(i,:))*norm(V(i,:))),1),-1);
    gamma(i) = 90 - real(acosd(CosTheta));
end
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

function CD = CD_func2(r,V)
    RE = 6371e3;
    h = norm(r) - RE; 
    M = norm(V)/atmos(h,13); 
     if M < 0.85
        CD = 0.6;
    else
        CD = 0.23+0.82/M^2-0.55/M^4;
        CD = 1.8*CD;
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
function [value,isterminal,direction] = stage1cond(t,U,m0,ms,Isp,g0,mprop)
    RE = 6371e3;
    m = U(7);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    
    m1 = m - m0(2);
    
    mp1 = m1 - ms(1);

    mp_req = m1 * ( 1 - exp( norm(V) / ( Isp(1) * g0 ) )  );

    dVmax = -g0*Isp(1)*log( (ms(1)+0.2*mprop(1))/ms(1) );

    if mp1 < 0.11*mprop(1)
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

function [value,isterminal,direction] = reentbcond(t,U,mf0,m0,ms,Vcrit)
    RE = 6371e3;
    m = U(7);
    r = U(1:3);
    V = U(4:6);
    h = norm(r) - RE;
    if m-ms(1) <= 500 %pure fuel mass
        value = 0;
        isterminal = 1;
        direction = 0;
    elseif norm(V) <= Vcrit
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

function [value,isterminal,direction] = landcond(t,U,mf0,ms,m0)
    RE = 6371e3;
    m = U(7);
    r = U(1:3);
    h = norm(r) - RE;
    V = U(4:6);
   if norm(V) < 1
        value = 0;
        isterminal = 1;
        direction = 0;
   elseif h<=0
        value = 0;
        isterminal = 1;
        direction = 0;
   else
        value = 1;
        isterminal = 0;
        direction = 0;
    end
end


function [value,isterminal,direction] = freefallcond(t,U,hstop)
    RE = 6371e3;
    m = U(7);
    r = U(1:3);
    h = norm(r) - RE;
    if h <= hstop
        value = 0;
        isterminal = 1;
        direction = 0;
    elseif h <= 0
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

function [value,isterminal,direction] = delay_cond(t,U,t0,delta_t)
    if (t-t0) > delta_t
           value = 0;
        isterminal = 1;
        direction = 0;
    else
        value = 1;
        isterminal = 0;
        direction = 0;
    end
end


function [value,isterminal,direction]  = break2cond(t,U,V_break,h_break)
    RE = 6371e3;
    m = U(7);
    r = U(1:3);
    h = norm(r) - RE;
    V = norm(U(4:6));
    if h <= h_break
        value = 0;
        isterminal = 1;
        direction = 0;
    elseif V <= V_break
        value = 0;
        isterminal = 1;
        direction = 0;
    elseif h <= 0
        value = 0;
        isterminal = 1;
        direction = 0;
    else
        value = 1;
        isterminal = 0;
        direction = 0;
    end
end 

function T = Tfunc(h,T0,Ae0,p_e,stage)
    p0 = 1013.25e2;
    p = atmos(h,1)*p0;
    T = T0(stage) + (p_e-p)*Ae0(stage);
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

function [value,isterminal,direction] = cruisecond(t,U)
    RE = 6371e3;
    d2r = pi/180;
    r = U(1:3);
    V = U(4:6);
    h = norm(r)-RE;
    gamma = gammafunc(r.',V.');
    if gamma-10 < 0.1
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

    if norm(V) >= 1.027*V_trgt
        value = 0;
        isterminal = 1;
        direction = 0;
        %disp('Circularizecond target velocity achived')
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