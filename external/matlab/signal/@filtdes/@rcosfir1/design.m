function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Set up design params
N = get(d,'order');

% Get frequency specs, they have been prenormalized
Fc = get(d,'Fc');

win = generatewindow(d);

tm = get(d,'TransitionMode');

R = get(d,tm);    

dt = get(d,'DesignType');

dt_opts = set(d,'DesignType');

% Translate design type when necessary
if strcmpi(dt,dt_opts{2})
    dt = 'sqrt';
end

b = firrcos(N,Fc,R,2,tm,dt,[],win);

% Construct object
Hd = dfilt.dffir(b);



