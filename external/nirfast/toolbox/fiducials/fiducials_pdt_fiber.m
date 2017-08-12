function [sources,detectors] = fiducials_pdt_fiber(fiducials,mesh)

% [sources,detectors] = fiducials_pdt_fiber(fiducials,mesh)
%
% finds source/detector locations based on fiducial locations
% and the mesh


dist = 1; % node distance (mm)

x1 = fiducials(1,1);
y1 = fiducials(1,2);
z1 = fiducials(1,3);

x2 = fiducials(2,1);
y2 = fiducials(2,2);
z2 = fiducials(2,3);

len = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2); % fiber length (mm)

xc = (x2+x1)/2;
yc = (y2+y1)/2;
zc = (z2+z1)/2;

v1 = sqrt((x1-xc)^2 + (y1-yc)^2 + (z1-zc)^2);
v2 = sqrt((x2-xc)^2 + (y2-yc)^2 + (z2-zc)^2);

x1 = xc + (x1-xc)*len/(2*v1);
y1 = yc + (y1-yc)*len/(2*v1);
z1 = zc + (z1-zc)*len/(2*v1);

x2 = xc + (x2-xc)*len/(2*v2);
y2 = yc + (y2-yc)*len/(2*v2);
z2 = zc + (z2-zc)*len/(2*v2);

numpoints = floor(sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)/dist);
dx = (x2-x1)/(numpoints-1);
dy = (y2-y1)/(numpoints-1);
dz = (z2-z1)/(numpoints-1);

sources = zeros(numpoints,3);
if dx == 0
    sources(:,1) = x1;
else
    sources(:,1) = (x1:dx:x2)';
end
if dy == 0
    sources(:,2) = y1;
else
    sources(:,2) = (y1:dy:y2)';
end
if dz == 0
    sources(:,3) = z1;
else
    sources(:,3) = (z1:dz:z2)';
end

detectors = sources(1,:);