function faultTolStruct = updateFaultTolStruct(faultTolStruct, ...
    newObjVal, updateUndefValue)
%

%updateFaultTolStruct Update fault tolerance structure
%
%   Optimization Toolbox algorithms with fault tolerance support handle
%   undefined points by either trying several trial points in one iteration
%   (e.g. levenbergMarquardt) or trying one point per iteration (e.g.
%   sfminbx, sfminle). Both types of algorithms are supported by
%   maintaining a fault tolerance structure, faultTolStruct.
%
%   faultTolStruct holds fault tolerance state information for the current
%   "iteration" of an algorithm. From the perspective of faultTolStruct, an
%   "iteration" is completed at the point in the algorithm when the
%   iterative display is updated. faultTolStruct has the following fields:
%
%     undefObj:    Boolean. true if objective value has evaluated to
%                  inf/NaN/complex at any point during the current
%                  iteration, false otherwise.
%     undefValue:  if faultTolStruct.undefObj is true,
%                  faultTolStruct.undefValue is either 'inf', 'NaN',
%                  or 'complex', dependent on how the objective value
%                  is undefined. If faultTolStruct.undefObj is false,
%                  undefValue is an empty string.
%     currTrialWellDefined : Boolean. true if the current trial point in
%                            the current iteration is well defined.
%     chkComplexObj: indicates whether a complex value should be
%                    treated as undefined
%
%   faultTolStruct = updateFaultTolStruct(faultTolStruct, newObjVal,
%   updateUndefValue) updates the fault tolerance structure given the
%   current algorithm state. faultTolStruct.undefValue is only updated in
%   subsequent updates if updateUndefValue is true. This is normally the
%   case if the user has set the display to 'iter'. This call should be
%   performed after every evaluation of the user function(s) in the main
%   iterative loop.

%   Copyright 2011 The MathWorks, Inc.

% Always update whether the current trial point is well defined or not.
% The second part of the AND statement allows the caller to optionally
% check whether a function value is complex.
faultTolStruct.currTrialWellDefined = isfinite(newObjVal) && ...
    (~faultTolStruct.chkComplexObj || isreal(newObjVal));

% If the current trial point is undefined then we update the remaining
% fields of the fault tolerance structure.
if ~faultTolStruct.currTrialWellDefined
    % Check objective function value
    faultTolStruct.undefObj = ~faultTolStruct.currTrialWellDefined;
    if updateUndefValue && faultTolStruct.undefObj
        if isnan(newObjVal)
            faultTolStruct.undefValue = 'NaN';
        elseif ~isfinite(newObjVal)
            faultTolStruct.undefValue = 'Inf';
        else % Result is complex
            faultTolStruct.undefValue = 'complex';
        end
    end
end



