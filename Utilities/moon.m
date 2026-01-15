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
    % inc - inclination (deg)
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
            %outputs Sun state relative to Earth for given time
            % Generated for one year at J2000 using cftool
            % t is in seconds since epoch, converted to juliandate (days)
            if size(t,1) ~=1
                t =t.';
            end
            t = t/86400 + 2.4515445e6; % seconds to jdate

            xFit = 3.7885e5*sin(0.23*t-1.0255);
            yFit = 3.591850364918779e+05*sin(24.262202968219331*t+1.715014781889868);
            zFit = 1.425312441205498e+05*sin(24.256011266183258*t+1.479531695405254);
            dxFit = 1.017081200050529*sin(24.262729146994818*t-1.457024036799175);
            dyFit = 0.951978269048581*sin(24.258848362721043*t-2.996524339132248);
            dzFit = 0.377544390751570*sin(24.252373513637487*t+3.049721911685454);
            rMoon = [xFit; yFit; zFit];
            vMoon = [dxFit; dyFit; dzFit];
        end
    end


end
