function sd_list = place_sd_line(p1,p2,spacing)
% inputs:
% sources will go from point 1 to point 2.
% spacing is a vector of distances from point 1 to each source.

p1p2 = p2 - p1;
% directions are:
    % x = x1 + at
    % y = y1 + bt
    % z = z1 + ct
    % where t=0 is starting point and t=1 is ending point.
    
dist = sqrt((p2(1)-p1(1))^2 + (p2(2)-p1(2))^2 + (p2(3)-p1(3))^2);    
spacing = spacing/dist;

i = 1;
for t = spacing    
    x = p1(1) + p1p2(1)*t;
    y = p1(2) + p1p2(2)*t;
    z = p1(3) + p1p2(3)*t;
    s = [x y z];
    sd_list(i,:) = s;
    i = i+1;
end