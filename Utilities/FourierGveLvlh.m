function B = FourierGveLvlh( oE, stop, primary)
%FourierGveLvlh approximates the GVE of a set of orbital elements.
%  Can be used for an arbitrary number of satellite states, of the same or
%  different satellites
%% Input Arguments
% 
% * oe   - 6xN matrix of orbital element vectors.
%          each 6x1 vector is comprised of [a e inc raan aop M].'
%          all angles are in radians
% * stop - criteria to stop evaluating series, two different cases             
%          integer - k at which to cut off the Fourier series
%          else - tolerance - size of Fourier coefficient at which to cut
%                             off series
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
%  1. Set all elements to value of 0 frequency coefficient w/o leading
%     coefficient
%  2. Loop over k until kMax or tolerance reached
%     - Calculate Jk, dJk
%     - Add k frequency elements to each element of B
%  3. Multiply each row of B by leading coefficient.
% 
%% Handle Input
if nargin < 3
    primary = earth();
end
if nargin < 2
    stop = 1e-14;
end
if mod(stop,1) == 0
    kMax = stop;
    tol = -inf;
else
    kMax = inf;
    tol = stop;
end

a    = oE(1,:);
e    = oE(2,:);
inc  = oE(3,:);
% raan = oE(4,:);
aop  = oE(5,:);
M    = oE(6,:);

eta  = sqrt(1-e.^2);
mu = primary.mu;


N = size(oE,2);
B = zeros(18,N);

%% 0 Frequency
B(2,:)  = 1;
B(5,:)  = -3/2*e;
B(9,:)  = -3/2*e.*cos(aop);
B(12,:) = -3/2*e.*sin(aop);
B(13,:) = -e;
B(16,:) = -3*e;

%% Loop k
k = 1;
maxCo = inf;
while (k < kMax) && (maxCo > tol)
    Jk  = besselj(k,k*e);
    dJk = k/2*(besselj(k-1,k*e) - besselj(k+1,k*e));
    sk  = sin(k*M);
    ck  = cos(k*M);

    B(1,:)  = B(1,:)  + 1/k*dJk.*sk;
    B(2,:)  = B(2,:)  + 2*Jk.*ck;
    B(4,:)  = B(4,:)  + dJk.*sk/k;
    B(5,:)  = B(5,:)  + (2*eta.^2./e.*Jk + 2*k^-2.*dJk).*ck;
    B(9,:)  = B(9,:)  + 2*cos(aop).*k^-2.*dJk.*ck - 2/k*eta.*sin(aop)./e.*Jk.*sk;
    B(12,:) = B(12,:) + 2*sin(aop).*k^-2.*dJk.*ck + 2/k*eta.*cos(aop)./e.*Jk.*sk;
    B(13,:) = B(13,:) + 2*eta.^2./e.*Jk.*ck;
    B(14,:) = B(14,:) + (2/k*eta.^2.*dJk + 2/k./e.*Jk).*sk;
    B(16,:) = B(16,:) + (2*eta.^4./e.*Jk + 4*k^-2*e.^2.*dJk).*ck;
    B(17,:) = B(17,:) + (2/k*eta.^2.*dJk + 2/k./e.*Jk).*sk;
    
    k = k + 1;
    maxCo = max([Jk, dJk]);
end

%% Leading Coefficients
B(1,:)  = B(1,:).*4.*a.^(3/2).*e/sqrt(mu);
B(2,:)  = B(2,:).*2.*a.^(3/2).*eta/sqrt(mu);
B(4,:)  = B(4,:).*2.*sqrt(a/mu).*eta.^2;
B(5,:)  = B(5,:).*sqrt(a/mu).*eta;
B(9,:)  = B(9,:).*sqrt(a/mu)./eta;
B(12,:) = B(12,:).*sqrt(a/mu)./eta./sin(inc);
B(13,:) = B(13,:).*-sqrt(a/mu).*eta./e;
B(14,:) = B(14,:).*sqrt(a/mu)./e;
B(15,:) = B(12,:).*-cos(inc); % Not a mistake
B(16,:) = B(16,:).*sqrt(a/mu)./e;
B(17,:) = B(17,:).*-sqrt(a/mu).*eta./e;
