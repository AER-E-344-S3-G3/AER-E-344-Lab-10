function [downstream_dist,nozzle_area] = getnozzleparams
%GETNOZZLEPARAMS Get areas and distances of ISU's de Laval nozzle.
%   [downstream_dist,nozzle_area] = GETNOZZLEPARAMS() returns the distance 
%   downstream of the throat and the area corresponding to each of the
%   pressure taps connected to ISU's de Laval nozzle.
%
%   Units:
%       downstream_dist [in]
%       nozzle_area [in^2]
downstream_dist = [-4.00,-1.50,-0.30,-0.18,0.00,0.15,0.30,0.45,0.60, ...
    0.75,0.90,1.05,1.20,1.35,1.45]; % [in]
nozzle_area = [0.800,0.529,0.480,0.478,0.476,0.497,0.518,0.539,0.560, ...
    0.581,0.599,0.616,0.627,0.632,0.634]; % [in^2]
end

