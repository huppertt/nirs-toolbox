function [F,A,P,nfpts] = super_validatespecs(this)
%SUPER_VALIDATESPECS   Validate the specs

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

% Get amplitudes and frequencies 
F = this.Frequencies;
H = this.FreqResponse;
nfpts = length(F);

if nfpts~=length(H),
    error(message('signal:fspecs:abstractsbarbmagnphase:super_validatespecs:InvalidSpecifications'))
end

% Force row vectors
A = abs(H);
P = angle(H);
F = F(:).';
A = A(:).';
P = P(:).';

% [EOF]
