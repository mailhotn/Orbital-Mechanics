classdef Moon
    %Moon defines a class containing Moon physical characteristics
    % all values are taken from Vallado
    %
    %  ~~~~~~~~~~   Primary Attributes    ~~~~~~~~~~~
    % mu - gravitational parameter (km^3/s^2)
    % Re - equatorial radius (km)
    % J2 - zonal harmonic
    % we - rotation rate (deg/s)
    % Tday - length of sidereal day (s)
    % Tyear - length of year (s) - orbital period
    %
    %  ~~~~~~~~~~   Third-Body Attributes    ~~~~~~~~~~~
    % sma - semimajor axis (km)
    % ecc - eccentricity
    % nMo - mean motion (1/s)
    properties (Constant)
        name = 'Moon';

        % Primary Attributes
        mu = 4902.799;
        Re = 1738.0;
        R = 1738.0;
        J2 = 2.027e-4;
        we = 360/27.32166/86400;
        Tday = 27.32166*86400;
        Tyear = 27.321582*86400;

        % Third-body Atributes (earth primary)
        sma = 384400;
        ecc = 0.05490;
        inc = 5.145396;
        nMo = 2*pi/27.321582/86400;
        mass = 7.3483e22;
    end

    methods (Static)
        function [rMoon, vMoon] = PosJ2000(t)
            %outputs Moon state relative to Earth for given time
            % Generated for one year at J2000 using cftool
            % t is in seconds since epoch, converted to juliandate (days)
            if size(t,1) ~=1
                t =t.';
            end
            t = t/86400;% + 2.4515445e6; % seconds to jdate (since epoch)
            
            
            xFit = 7.1825e3 -3.2679e5*cos(0.23*t) +1.9789e5*sin(0.23*t);
            yFit = -2.7155e4 -1.9623e5*cos(0.23*t) -2.9944e5*sin(0.23*t);
            zFit = -1.1192e4 -4.9359e4*cos(0.2299*t) -1.3319e5*sin(0.2299*t);
            dxFit = -5.5187e-5 +0.5265*cos(0.23*t) +0.8702*sin(0.23*t);
            dyFit = 2.6108e-4 -0.7952*cos(0.23*t) +0.5234*sin(0.23*t);
            dzFit = -1.1777e-4 -0.3538*cos(0.2299*t) +0.1317*sin(0.2299*t);
            rMoon = [xFit; yFit; zFit];
            vMoon = [dxFit; dyFit; dzFit];
        end
    end


end
