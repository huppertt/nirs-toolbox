function n = secorder(Hd)
%SECORDER Returns the order of each section.

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

s = Hd.sosMatrix;
nsec = nsections(Hd);

% Second order sections by default
n = 2*ones(nsec,1);

for i=1:nsec,
    if (s(i,3)==0 && s(i,6)==0),
        % First order sections
        n(i) = 1;
        if (s(i,2)==0 && s(i,5)==0),
            % Zero th order sections
            n(i) = 0;
        end
    end
end


% [EOF]
