function checkoptionsize(option, optionsize, numvars, numfunelts)
%

%CHECKOPTIONSIZE Verify problem size dependent options 
%   CHECKOPTIONSIZE('OPTION', OPTIONSIZE, NUMVARS) verifies that the
%   specified option is of the correct size. Valid values for 'OPTION' are
%   'TypicalX' and 'HessPattern'. OPTIONSIZE is the size of the specified
%   OPTION. NUMVARS specifies the number of free variables (normally
%   numel(X0)).
%
%   CHECKOPTIONSIZE('JacobPattern', OPTIONSIZE, NUMVARS, NUMFUNELTS)
%   verifies that the option 'JacobPattern' is of the correct size.
%   NUMFUNELTS is the number of elements in the vector (or matrix) returned
%   by FUN.
%
%   This is a helper function for the Optimization Toolbox.

%   Copyright 2008-2011 The MathWorks, Inc.

% Perform size checks
switch lower(option)
    case 'typicalx'
        if prod(optionsize) ~= numvars
            ME = MException('optimlib:checkoptionsize:InvalidSizeOfTypicalX', ...
                getString(message('optimlib:commonMsgs:InvalidSizeOfTypicalX')));
            throwAsCaller(ME);
        end
    case 'hesspattern'
        expsize = [numvars numvars];
        if ~isequal(optionsize, expsize)
            ME = MException('optimlib:checkoptionsize:InvalidSizeOfHessPattern', ...
                getString(message('optimlib:commonMsgs:InvalidSizeOfHessPattern',numvars,numvars)));
            throwAsCaller(ME);
        end
    case 'jacobpattern'
        expsize = [numfunelts numvars];
        if ~isequal(optionsize, expsize)
            ME = MException('optimlib:checkoptionsize:InvalidSizeOfJacobPattern', ...
                getString(message('optimlib:commonMsgs:InvalidSizeOfJacobPattern',numfunelts,numvars)));
            throwAsCaller(ME);
        end
end

