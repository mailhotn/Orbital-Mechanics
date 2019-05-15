function [ H ] = hemid( phi )
%hemid is the hemisphere function with input in degrees
%   See Orbit Design and Management - James R Wertz Chapter 8

if mod(phi,360) < 180
    H = 1;
else
    H = -1;
end
end

