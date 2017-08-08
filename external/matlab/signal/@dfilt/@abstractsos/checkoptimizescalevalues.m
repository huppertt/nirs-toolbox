function variables = checkoptimizescalevalues(this,variables)
%CHECKOPTIMIZESCALEVALUES check if optimize scale values is possible

%   Copyright 2009-2010 The MathWorks, Inc.


issvnoteq2one = this.issvnoteq2one;
if this.OptimizeScaleValues && ~all(issvnoteq2one),
    % Unit scale values cannot be skipped when specified through a port
    warning(message('signal:dfilt:abstractsos:checkoptimizescalevalues:UnitScaleValues'));
    that = copy(this);
    that.OptimizeScaleValues = false;
    g = that.privScaleValues.';
    variables{3} = g;
end


% [EOF]
