function [ X_ECEF ] = eci2ecef( X_ECI, GMST )
%eci2ecef Converts from ECI to ECEF
%   Rotates through the angle GMST, Velocity is adjusted according to the
%   differentiation of multiple law: d(T*R)/dt = T*dR/dT + dT/dt*R
w_e = 7.2921150e-5; % rad/sec
R = X_ECI(1:3,:);
V = X_ECI(4:6,:);

T = [cosd(GMST), sind(GMST), 0;
    -sind(GMST), cosd(GMST), 0;
              0,          0, 1];

T_dot = w_e*[-sind(GMST),  cosd(GMST), 0;
             -cosd(GMST), -sind(GMST), 0;
                       0,           0, 0];

R1 = T*R;
V1 = T*V + T_dot*R;
X_ECEF = [R1; V1];
end

