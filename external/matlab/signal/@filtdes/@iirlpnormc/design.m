function varargout = design(d)
%Design  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.


args = setupdesignparams(d);

R = get(d,'maxRadius');

args = {args{1:6},R,args{7:end}};

[b,a,err,s,g] = iirlpnormc(args{:});

% Construct object
h = dfilt.df2sos(s,g);

% Add the dynamic property to store the mask commands
p = schema.prop(h, 'MaskInfo', 'MATLAB array');
set(p, 'Visible', 'Off');
set(h, 'MaskInfo', maskinfo(d));

if nargout,
    
    varargout = {h};
else,
    drawmasknresp(d, h);
end

setfdesign(h, getfdesign(d.ResponseTypeSpecs));
