function [ Bxx ] = vec2sparse( V, dimu, indu )
%vec2sparse Creates a matrix in which the elements of a vector are placed
%in blocks along the diagonal such that in each block the corresponding
%element of the diagonal is in the required position.
%This is useful for writing large matrix equations.
%% EXAMPLE: 
% vec2sparse([1,2,3],[3,2],[1,2]) would produce:
% Bxx = [0, 1, 0, 0, 0, 0;
%        0, 0, 0, 0, 0, 0;
%        0, 0, 0, 0, 0, 0;
%        0, 0, 0, 2, 0, 0;
%        0, 0, 0, 0, 0, 0;
%        0, 0, 0, 0, 0, 0;
%        0, 0, 0, 0, 0, 3;
%        0, 0, 0, 0, 0, 0;
%        0, 0, 0, 0, 0, 0];

%% Input Arguments
% V    - length N vector of values
% dimu - 1x2 vector of block dimensions
% indu - 1x2 vector of values index relative to block
%%
if size(dimu,1) == 2
    dimu = dimu.';
end
if size(indu,1) == 2
    indu = indu.';
end
N = length(V);
Bxx = zeros(N*dimu(1),N*dimu(2));
indval = (indu.' + (0:N-1).*dimu.').';
Bxx(indval(:,1),indval(:,2)) = diag(V);
end

