function [ r_s ] = sun_dir( day )
% Approximate Direction of sun in ECI
RA_sun = 2*pi/365.25*(day-80);
dec_sun = asin(0.39795*cos(2*pi/365.25*(day-173)));
r_s = [cos(dec_sun)*cos(RA_sun); cos(dec_sun)*sin(RA_sun); sin(dec_sun)];
end

