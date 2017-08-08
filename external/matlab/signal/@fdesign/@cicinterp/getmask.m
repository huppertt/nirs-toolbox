function [F, A] = getmask(this, fcns, R, specs)
%GETMASK   Get the mask.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

% If the specs were not passed in or are [], use the design specifications.
if nargin < 4 || isempty(specs)
    specs = getspecs(this.CurrentSpecs);
end

% In order to compute the nominal gain of CICs, we need R, M & N.  N is
% computed in the @fdfmethod/abstractdesigncic/abstract_design.m method.
fm = getfmethod(this,'multisection');

M = this.DifferentialDelay;
Fp = specs.Fpass;
Aa = specs.Astop;
N = abstract_design(fm,R,M,Fp,Aa,'interp');

% Call shared code to compute mask info.
[F,A] = abstract_cicmask(this,fcns,R,M,N,Fp,Aa);


% [EOF]
