function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get specs
N = get(d,'order');
F = get(d,'FrequencyVector'); % Freqs have been prenormalized
E = get(d,'FrequencyEdges'); % Freqs have been prenormalized
G = get(d,'GroupDelayVector');
W = get(d,'WeightVector');
R = get(d,'maxRadius');
P = [get(d,'initPnorm'), get(d,'Pnorm')];
D = get(d,'DensityFactor');
ID = get(d,'initDen');

% Setup optional initial denominator argument
optarg = {};
if ~isempty(ID),
    optarg = {ID};
end

% Design filter
[b,a,err,s] = iirgrpdelay(N,F,E,G,W,R,P,{D},optarg{:});

Hd = dfilt.df2sos(s);

