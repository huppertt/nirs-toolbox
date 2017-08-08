function update(h)
%UPDATE   

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

% Delete filterquantizer
for i=1:length(h.privfq),
    classname{i} = class(h.privfq(i));
end
h.privfq(strmatch(class(h.filterquantizer),classname)) = [];

% Set the current filterquantizer
h.filterquantizer = update(h.filterquantizer);
set_ncoeffs(h.filterquantizer, naddp1(h));

% Add filterquantizer
h.privfq = [h.privfq; h.filterquantizer];

% [EOF]
