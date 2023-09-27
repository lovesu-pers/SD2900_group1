clear
close all
clc
%% Delta-V calculations
d2r = pi/180;
r2d = 180/pi;

mu = 3.986e5 * (1e3)^3;
omega = [0;0;7.292115855377074e-5;];
RE = 6371e3;
g0 = 9.8;

h0 = 2;
lat0 = 5.167713*d2r;
long0 = -52.683994*d2r;
r0 = latlong2cart(lat0, long0, h0);


h_target = 1000e3;
h_park = 287e3;


dV_grav = 2500;
dV_drag = 150;
dV_rot = norm(cross(omega,r0));

dV_park = sqrt(mu/(RE+h_park));
dV_burn1 = sqrt(2*mu/(h_park+RE)-2*mu/(h_target+h_park+2*RE)) ... 
    - sqrt(mu/(RE+h_park));
dV_burn2 = sqrt(mu/(RE+h_target))-sqrt(2*mu/(RE+h_target)-2*mu/(h_target+h_park+2*RE));

dV_tot = -dV_rot + dV_grav+dV_drag+dV_park+dV_burn1+dV_burn2;
%% Expendable rocket optimization
Isp1 = 260; Ve1 = Isp1*g0;
Isp2 = 350; Ve2 = Isp2*g0;
beta = [1;Ve2/Ve1];

Vnondim = dV_tot/(Ve1);

pi1 = 0.1;
pi2 = 0.12;
alpha = [1;1.2];

epsilon1 = 0.07;
epsilon2 = 0.07;


pi1_opt = Nstage(Vnondim,beta,[epsilon1;epsilon2],alpha,pi1)

pi_TOT = (pi1_opt)^2 * prod(alpha)

mPL = 1e3;

m01_opt = mPL/pi_TOT
m02_opt = pi1_opt*m01_opt

mp1_opt = (m01_opt-m02_opt) * (1-epsilon1);
mp2_opt = (m02_opt-mPL) * (1-epsilon2);

ms1_opt = m01_opt - m02_opt - mp1_opt;
ms2_opt = m02_opt-mp2_opt-mPL;

mf1_opt = m01_opt - mp1_opt;
mf2_opt = m02_opt - mp2_opt;

mu1_opt = mf1_opt/m01_opt
mu2_opt = mf2_opt/m02_opt

dV1_opt = -g0*Isp1*log(mu1_opt)
dV2_opt = -g0*Isp2*log(mu2_opt)
dV_tot_opt = dV2_opt+dV1_opt;
%% Newton iteration to find optimum payloadratio of first stage
function p=Nstage(vf,beta,epsilon,alpha,pi_geuss)
    N = size(beta,1);
    p = pi_geuss;
    f = vf;
    tol = 1e-9;
    for k=1:N
        f = f+beta(k)*log(epsilon(k)+alpha(k)*(1-epsilon(k))*p);
    end
    while abs(f)>tol
        f = vf;
        fp = 0;
        for k=1:N
            f = f+beta(k)*log(epsilon(k)+alpha(k)*(1-epsilon(k))*p);
            fp = fp+alpha(k)*beta(k)/(epsilon(k)+alpha(k)*(1-epsilon(k))*p);
        end
        d = -f/fp;
        p = p+d
    end
end

function r = latlong2cart(lat,long,h)
    RE = 6371e3;
    R = RE+h;
    r = R*[cos(lat)*cos(long); ...
            cos(lat)*sin(long); ...
            sin(lat)];
end
