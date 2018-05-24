%% Initializaton
% Simulation Parameters
T_vec    = 0:10:5*86400;
lat_gs = 30;
lon_gs = 0;
e_min  = 15;
reltol = 1e-8;
abstol = 1e-9;
% Initialize data vectors

Ni = 50;
N  = 60;
p = primes(N);
T        = zeros(1,N-Ni+1);
P        = zeros(1,N-Ni+1);
F        = zeros(1,N-Ni+1);
inc      = zeros(1,N-Ni+1);
alt      = zeros(1,N-Ni+1);
GMST0    = zeros(1,N-Ni+1);
mean_fit = zeros(1,N-Ni+1);
osc_fit  = zeros(1,N-Ni+1);

datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';

%% Go over all solutions & sim Osculating to get new fitness
for ii = (Ni+1):(N+1)
    if ~any((ii-1)==p)
        load([datafolder '\Walker_Mean_Sol_T_' num2str(ii-1) '.mat']);
        if sol.fit == 20
            T(ii-Ni)        = sol.T;
            P(ii-Ni)        = nan;
            F(ii-Ni)        = nan;
            inc(ii-Ni)      = nan;
            alt(ii-Ni)      = nan;
            GMST0(ii-Ni)    = nan;
            mean_fit(ii-Ni) = nan;
            osc_fit(ii-Ni)  = nan;
        else
            T(ii-Ni)        = sol.T;
            P(ii-Ni)        = sol.P;
            F(ii-Ni)        = sol.F;
            inc(ii-Ni)      = sol.inc;
            alt(ii-Ni)      = sol.alt;
            GMST0(ii-Ni)    = sol.GMST0;
            mean_fit(ii-Ni) = sol.fit;
            osc_fit(ii-Ni)  = WalkerFitness_WGS84_osc([P(ii-Ni),F(ii-Ni),inc(ii-Ni)...
                ,alt(ii-Ni),GMST0(ii-Ni)],T(ii-Ni),T_vec,lat_gs,lon_gs,e_min,reltol,abstol);
        end
    else
        T(ii-Ni)        = ii-1;
        P(ii-Ni)        = nan;
        F(ii-Ni)        = nan;
        inc(ii-Ni)      = nan;
        alt(ii-Ni)      = nan;
        GMST0(ii-Ni)    = nan;
        mean_fit(ii-Ni) = nan;
        osc_fit(ii-Ni)  = nan;
    end
end
%% Plot Results
%Fitness Plot
figure(1)
plot(T,mean_fit,'o',T,osc_fit,'o');
xlabel('# Satellites')
ylabel('Maximum PDOP')
legend('Mean Sim','Osculating Sim')
% inclination
figure(2)
plot(T,inc,'o');
xlabel('# Satellites')
ylabel('Inclination [°]')
% Altitude
figure(3)
plot(T,alt,'o');
xlabel('# Satellites')
ylabel('Altitude [km]')