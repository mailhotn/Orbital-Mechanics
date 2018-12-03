function [ rho, hBase, rhoBase, hScale] = EarthAtmosphere( r )
%atmo_model Calculates the density of the atmosphere at a given radius
h = r - 6378;
A = [0 25 1.225 7.249;
    25 30 0.038937 6.349;
    30 40 0.017715 6.682;
    40 50 0.003966 7.554;
    50 60 0.001056 8.382;
    60 70 0.00032 7.714;
    70 80 8.76E-05 6.549;
    80 90 1.9E-05 5.799;
    90 100 3.39E-06 5.382;
    100 110 5.29E-07 5.877;
    110 120 9.65E-08 7.263;
    120 130 2.43E-08 9.472;
    130 140 8.47E-09 12.636;
    140 150 3.84E-09 16.149;
    150 180 2.07E-09 22.523;
    180 200 5.46E-10 29.74;
    200 250 2.78E-10 37.105;
    250 300 7.24E-11 45.546;
    300 350 2.41E-11 53.628;
    350 400 9.5E-12 53.298;
    400 450 3.72E-12 58.515;
    450 500 1.58E-12 60.828;
    500 600 6.96E-13 63.822;
    600 700 1.45E-13 71.835;
    700 800 3.61E-14 88.667;
    800 900 1.17E-14 124.64;
    900 1000 5.24E-15 181.05];
ind = (h >= A(:,1)) & (h < A(:,2));
if ~isempty(ind)
    hBase = A(ind,1);
    rhoBase = A(ind,3);
    hScale = A(ind,4);
    rho = rhoBase*exp(-(h-hBase)/hScale);
else
    rho = 0;
    rhoBase = 0;
    hBase = 0;
    hScale = inf;
end