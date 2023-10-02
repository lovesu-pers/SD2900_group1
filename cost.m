t = linspace(0,50);
Nper = 15;
N = t * Nper;
N_ex = N;
N_re = N/100;
%t year;
%N launch total;
%Ner launch per year;

%development;;
Cdevelopment_re = 6000;
Cdevelopment_ex = 1000;

%vehicle
Cvehicle_re0 = 1000;
Cvehicle_re = Cvehicle_re0 * N_re;
Cvehicle_ex0 = 50;
Cvehicle_ex = Cvehicle_ex0 * N_ex;

%total flight opretions per flight;
%Re has more complex systems, more expensive check-out and recovery;
Cflightops_re0 = 5;
Cflightops_ex0 = 2;
Cflightops_re = Cflightops_re0 * N;
Cflightops_ex = Cflightops_ex0 * N;

%launch insurance;
Cinsurance_re0 = 4.16;
Cinsurance_ex0 = 0.86;
Cinsurance_re = Cinsurance_re0 * N;
Cinsurance_ex = Cinsurance_ex0 * N;

%recurring cost of recovery(reusable only);
Crecovery_re = 0.2 * Cflightops_re;

%refurbishment cost(reusable only);
Crefurb_re = 0.02 * Cvehicle_re;


%model the total launch cost as the sum of six individual components;
Claunch_re = Cdevelopment_re + Cvehicle_re + Cflightops_re + Crecovery_re + Crefurb_re + Cinsurance_re;
Claunch_ex = Cdevelopment_ex + Cvehicle_ex + Cflightops_ex + Cinsurance_ex;
%where;
%Claunch is total cost of launch in FY00 dollars(i.e., excluding
%inflation?);
%Cdevelopment is amortization of nonrecurring development cost;
%Cvehicle.re is amortization of vehicle production cost;
%Cvehicle.ex is recurring production cost(theoretical first unit
%cost reduced by learning curve);
%Cflightops is total cost of flight opretions per flight;
%Crecovery is recurring cost of recovery(reusable only);
%Crefurb is refurbishment cost(reusable only);
%Cinsurance is cost of launch insurance;

%cost using million dollars to calculate;

plot (t,Claunch_re,t,Claunch_ex);