function [xHans,xHans0] = hansenCoeffs(e,n,s,kMax,iMax)

b = (1-sqrt(1-e^2))/e;
cPhi = (1-b^2)/(1+b^2);
xHans = nan(1,kMax);

for k = 1:kMax
    xk = 0;
    mu = k*cPhi;
    for i = 0:iMax
        M = 0;
        p = i+(abs(k-s)+k-s)/2;
        for m = 0:p
            
                M = M + (-1)^m*nchoosek(sym(n-k+p+1),sym(n-k+m+1))*mu^m/factorial(m);
            
        end
        N = 0;
        q = i+(abs(k-s)-(k-s))/2;
        for m = 0:q
            
                N = N + nchoosek(sym(n+k+q+1),sym(n+k+m+1))*mu^m/factorial(m);
            
        end
        xk = xk + M*N*b^(2*i);
    end
    xk = ((1-b^2)^(2*n+3))/((1+b^2)^(n+1))*(-b)^abs(k-s)*xk;
    xHans(k) = xk;
end
xHans0 = (-e)^s*(1+s*sqrt(1-e^2))/(1+sqrt(1-e^2))^s;