function [sources,detectors] = fiducials_identity(fiducials,mesh)

% [sources,detectors] = fiducials_identity(fiducials,mesh)
%
% finds source/detector locations based on fiducial locations
% and the mesh

sources = fiducials;
detectors = fiducials;