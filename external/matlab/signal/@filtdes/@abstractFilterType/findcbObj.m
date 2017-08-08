function cbt = findcbObj(h,className)
%FINDCBOBJ Find object corresponding to callback.
%
%   Inputs:
%       h - handle to this object
%       className - Name of the class of the object we intend to find


%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

specObjs = get(h,'specobjs');

cbt = find(specObjs,'-isa', className);

