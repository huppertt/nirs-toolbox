function [sources,detectors] = fiducials_FD_triangle(fiducials,mesh)

% [sources,detectors] = fiducials_identity(fiducials,mesh)
%
% finds source/detector locations based on fiducial locations
% and the mesh

% this 

f = fiducials;

sd1 = place_sd_line(f(1,:),f(2,:),[8.89 8.89+12.7*[1:7]]);
sd2 = place_sd_line(f(3,:),f(4,:),[8.89 8.89+12.7*[1:3]]);
sd3 = place_sd_line(f(5,:),f(6,:),[8.89 8.89+12.7*[1:3]]);
sd = [sd1; sd2; sd3];



sources = sd;
detectors = sd;