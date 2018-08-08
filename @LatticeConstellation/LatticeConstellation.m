classdef LatticeConstellation < Constellation
    %LatticeConstellation Defines a 3D Lattice Flower Constellation
% Based on the work in J. J. Davis, · Martín, E. Avendaño, and D. Mortari,
% “The 3-D lattice theory of Flower Constellations,” , 2013.
% This implementation creates a constellation with a repeating earth ground
% track.
    
    properties (SetAccess = private)
        % General Architecture
        % nSats (in superclass)
        nPlanes     % N_o
        nAops       % N_w 
        nSatsPerAop % N_so'
        nRepeats    % N_p
        nDays       % N_d
        
        % Phasing
        nC1
        nC2
        nC3
        
        % 
        inc
        ecc
        
        % Initial Conditions
        M1
        raan1
        aop1
        
        % sma
    end
    
    methods
    end
    
end

