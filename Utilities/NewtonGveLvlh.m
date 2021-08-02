function B = NewtonGveLvlh( oE, tol, primary)
%NewtonGveLvlh calculates the GVE of a set of orbital elements.
%  Can be used for an arbitrary number of satellite states, of the same or
%  different satellites
%% Input Arguments
% 
% * oe   - 6xN matrix of orbital element vectors.
%          each 6x1 vector is comprised of [a e inc raan aop M].'
%          all angles are in radians
% * tol  - tolerance for N-R calculation of f
% * primary - primary body structure e.g. output of earth.m
%
%% Output Arguments
%
% * B - 18xN matrix of coefficients such that
%         
%       [da;          [B1  B2  B3;
%        de;           B4  B5  B6;     [u_r;    
%        dinc;  ==     B7  B8  B9;   *  u_s;
%        draan;        B10 B11 B12;     u_w]
%        daop;         B13 B14 B15;
%        dM]       n + B16 B17 B18] 
%                 
%% Algorithm
%
%  1. Calculate true anomaly f with N-R
%  2. Calculate B
% 
%% Handle Input
if nargin < 3
    primary = earth();
end
if nargin < 2
    tol = 1e-14;
end

a    = oE(1,:);
e    = oE(2,:);
inc  = oE(3,:);
% raan = oE(4,:);
aop  = oE(5,:);
M    = oE(6,:);

f    = pi/180*me2ta(180/pi*M,e,tol);
aol  = aop + f;

p = a.*(1-e.^2);
r = p./(1 + e.*cos(f));
b = a.*sqrt(1-e.^2);
mu = primary.mu;
h = sqrt(mu*p);


%% Calculate B
N = size(oE,2);
B = zeros(18,N);

B(1,:)  = 2*a.^2./h.*e.*sin(f);
B(2,:)  = 2*a.^2./h.*p./r;
B(4,:)  = p./h.*sin(f);
B(5,:)  = ((p + r).*cos(f) + r.*e)./h;
B(9,:)  = r./h.*cos(aol);
B(12,:) = r./h./sin(inc).*sin(aol);
B(13,:) = -p./h./e.*cos(f);
B(14,:) = (p + r)./h./e.*sin(f);
B(15,:) = -r./h./sin(inc).*cos(inc).*sin(aol);
B(16,:) = b./a./h./e.*(p.*cos(f) - 2*r.*e);
B(17,:) = -b./a./h./e.*(p + r).*sin(f);
