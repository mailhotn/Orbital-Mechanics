% wwwwwwwwwwwwwwwwwwwwwww
function map = light_gray
% wwwwwwwwwwwwwwwwwwwwwww
%{
This function creates a color map for displaying the planet as light
gray with a black equator.
r - fraction of red
g - fraction of green
b - fraction of blue
%}
% -----------------------
r = 0.8; g = r; b = r;
map = [r g b;
       0 0 0;
       r g b];
end %light_gray
% wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww