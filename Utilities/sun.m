classdef Sun
    %Sun defines a class containing Sun physical characteristics
    % all values are taken from Vallado
    %
    %  ~~~~~~~~~~   Attributes    ~~~~~~~~~~~
    % mu - gravitational parameter (km^3/s^2)
    % Re - equatorial radius (km)
    % J2 - zonal harmonic
    %
    %  ~~~~~~~~~~   Third-Body Attributes    ~~~~~~~~~~~
    % sma - semimajor axis (km)
    % ecc - eccentricity
    % inc - inclination (deg)
    % nMo - mean motion (1/s)

    properties (Constant)
        name = "Sun";
        mu = 1.32712428e11;
        Re = 696000;
        R = 696000;

        % Third-body Attributes
        sma = 149598023.0;
        ecc = 0.016708617;
        nMo = 2*pi/(365.2421897*86400);
        mass = 1.9891e30;

    end

    methods (Static)
        function [rSun, vSun] = PosJ2000(t)
            %outputs Moon state relative to Earth for given time
            % Generated for one year at J2000 using cftool
            % t is in seconds since epoch, converted to juliandate (days)
            if size(t,1) ~=1
                t =t.';
            end
            t = t/86400;% + 2.4515445e6; % seconds to jdate (since epoch)


            xFit = -7.5785e5 +2.7427e7*cos(0.0171*t) +1.4681e8*sin(0.0171*t);
            yFit = 2.0973e6 -1.3508e8*cos(0.0170*t) +2.7547e7*sin(0.0170*t);
            zFit = 9.0934e5 -5.8562e7*cos(0.0170*t) +1.1942e7*sin(0.0170*t);
            dxFit = 0.5602 +29.2777*cos(0.0169*t) -6.8244*sin(0.0169*t);
            dyFit = 0.0306 +5.2902*cos(0.0171*t) +26.7240*sin(0.0171*t);
            dzFit = 0.0132 +2.2935*cos(0.0171*t) +11.5862*sin(0.0171*t);
            rSun = [xFit; yFit; zFit];
            vSun = [dxFit; dyFit; dzFit];
        end
    end
end