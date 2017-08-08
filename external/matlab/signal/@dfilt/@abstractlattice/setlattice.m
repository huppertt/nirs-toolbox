function lattice = setlattice(Hd, lattice)
%SETLATTICE Overloaded set on the Lattice property.
  
%   Author: V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

% Always make coeffs a row
lattice = set_coeffs(Hd, lattice);

% Store lattice as reference and check data type
set(Hd,'reflattice',lattice);

ncoeffs = Hd.ncoeffs;
oldlength=0;
if ~isempty(ncoeffs), oldlength = ncoeffs(1); end
newlength = length(lattice);
ncoeffs(1) = newlength;
Hd.ncoeffs = ncoeffs;

% Update the ncoeffs property of the plugin.
set_ncoeffs(Hd.filterquantizer, ncoeffs);

if isempty(oldlength) || newlength~=oldlength,
    reset(Hd);
end

% Quantize the coefficients
quantizecoeffs(Hd);

% Hold an empty to not duplicate storage
lattice = [];
