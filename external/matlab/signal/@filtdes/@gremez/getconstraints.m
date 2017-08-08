function varargout = getconstraints(d)
%GETCONSTRAINTS Convert the error approximation vector
%   GETCONSTRAINTS(D) Converts the ErrorBands property to the format that
%   GREMEZ expects.  If 3 outputs are requested GETCONSTRAINTS will return
%   the parameter name, value and description to be used by GENMCODE.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

iae = get(d, 'ErrorBands');
if ~isempty(iae)
    
    % There must be more than one error band
    if any(diff(iae)) && length(iae) == nbands(d.ResponseTypeSpecs, d),
        iae = cellstr(num2str(iae'))';
        for indx = 1:length(iae),
            iae{indx} = ['e' iae{indx}];
        end
    else
        iae = [];
    end
end

if nargout == 1,
    varargout = {iae};
else
    if ~isempty(iae), iae = genmcodeutils('formatcellstr', iae); end
    varargout = {'IAEs', iae, 'Independent Approximation Errors'};
end

% [EOF]
