function faultTolStruct = createFaultTolStruct(chkComplexObj)
%

%createFaultTolStruct Create fault tolerance structure
%
%   faultTolStruct = createFaultTolStruct(chkComplexObj) creates a fault
%   tolerance structure. chkComplexObj indicates whether a complex value
%   should be treated as undefined or not. The least squares algorithms do
%   not consider complex values to be undefined. This call should be
%   performed before the main iterative loop.

%   Copyright 2011 The MathWorks, Inc.

faultTolStruct = struct(...
    'undefObj', false, ...
    'undefValue', '', ...
    'currTrialWellDefined', true, ...
    'chkComplexObj', chkComplexObj);
