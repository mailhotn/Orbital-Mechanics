function [ rho, hBase, rhoBase, hScale] = EarthAtmosphere( r )
%atmo_model Calculates the density of the atmosphere at a given radius
primary = earth();
h = r-primary.Re;
A = primary.ExpAtmoModel;

ind = (h >= A(:,1)) & (h < A(:,2));
if any(ind)
    hBase = A(ind,1);
    rhoBase = A(ind,3);
    hScale = A(ind,4);
    rho = rhoBase*exp(-(h-hBase)/hScale);
else
    rho = nan;
    rhoBase = nan;
    hBase = nan;
    hScale = nan;
end