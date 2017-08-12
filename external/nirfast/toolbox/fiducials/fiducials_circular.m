function [sources,detectors] = circular_SD(fiducials,mesh)

ri = 15;
i = 1;
for thi=0:2*3.1415/180:2*3.141
    [x(i),y(i)] = pol2cart(thi,ri);
    i = i + 1;
end

sources = [x' y'];
detectors = sources;