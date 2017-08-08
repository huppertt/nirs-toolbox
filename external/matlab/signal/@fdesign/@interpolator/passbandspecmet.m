function flag = passbandspecmet(Hf,Hd,ng)
%PASSBANDSPECMET Check whether passband response is within spec.
%   This should be a private method.

%   Author(s): R. Losada
%   Copyright 2009 The MathWorks, Inc.

fd = get(Hf, 'CurrentFDesign');

flag = passbandspecmet(fd,Hd,ng);


% [EOF]
