function varargout = set2int(this,coeffWL,inWL)
%SET2INT Scales the coefficients to integer numbers.
%   SET2INT(Hd) scales the coefficients to integer numbers and set the
%   coefficients and input fraction length to zero.
%
%   SET2INT(Hd, COEFFWL) uses the number of bits defined by COEFFWL to
%   represent the coefficients.
%
%   SET2INT(Hd, COEFFWL, INWL) uses the number of bits defined by COEFFWL
%   to represent the coefficients and the number of bits defined by INWL to
%   represent the input.
%
%   G = SET2INT(...) returns the gain G introduced in the filter by scaling
%   the coefficients to integer numbers. G is always a power of 2.
%
%   EXAMPLE: 
%
%   % Part 1: Comparing the step response of a filter in fractional and integer modes
%             b = rcosdesign(0.25,25,4);
%             h = dfilt.dffirt(b);
%             h.Arithmetic = 'fixed';
%             h.InputFracLength = 0; % Integer inputs
%             x = ones(100,1);
%             yfrac = filter(h,x); % Fractional mode
%
%             g = set2int(h);
%             yint = filter(h,x);  % Integer mode
%             % Rescale integer output
%             ysint1 = double(yint)/g;
%             % Verify that the scaled integer output is equal to the 
%             % fractional output
%             max(abs(ysint1-double(yfrac)))
% 
%   % Part 2: Reinterpreting the output binary   
%             % Another way to put the input and the output on the same scale
%             % consists in weighing the MSBs equally. 
%             WL = yint.WordLength;
%             FL = yint.FractionLength + log2(g);
%             ysint2 = fi(zeros(size(yint)),true,WL,FL);
%             ysint2.bin = yint.bin;
%             max(abs(double(ysint2)-double(yfrac)))

%   Copyright 2004-2013 The MathWorks, Inc.

narginchk(1,3);
g = 1;

if ~isfield(get(this),'Arithmetic'),
    error(message('signal:dfilt:dffirt:set2int:FixedPointLicense', 'SET2INT'));
end

if ~strcmpi(this.Arithmetic, 'fixed'),
    error(message('signal:dfilt:dffirt:set2int:FixedPointFcn', 'Arithmetic', 'fixed'));
end

if ~strcmpi(this.FilterInternals,'FullPrecision'),
    error(message('signal:dfilt:dffirt:set2int:FixedPointFcn', 'FilterInternals', 'FullPrecision'));
end

if nargin<3, inWL = this.InputWordLength; end
if nargin<2, coeffWL = this.CoeffWordLength; end

Href = reffilter(this);
b = Href.Numerator;
bfi = fi(b,this.Signed,coeffWL);
coeffFL = bfi.FractionLength;
g = 2^coeffFL;
this.Numerator = b*g;
this.Arithmetic = 'fixed';
this.CoeffWordLength = coeffWL;
this.CoeffAutoScale = 0;
this.NumFracLength = 0;
this.InputWordLength = inWL;
this.InputFracLength = 0;

if nargout>0, varargout{1} = g; end


% [EOF]
