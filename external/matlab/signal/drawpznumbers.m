function h = drawpznumbers(p, hax, varargin)
%DRAWPZNUMBERS Draw the PZNumbers
%   DRAWPZNUMBERS(PZ, HAX) Draw the multiplicity of PZ to the axes HAX.  

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 2, hax = newplot; end

fuzz = diff(get(hax, 'XLim'))/80;
h    = [];

[r,c]=size(p);
if (r>1) && (c>1), PEE=p;
else               PEE=p(:); c = min(r,c); end;

for which_col = 1:c,      % for each column of ZEE ...
    p = PEE(:,which_col);
    
    % Remove NaN's from p (they wouldn't be plotted and they'll break mpoles)
    p(isnan(p)) = [];
    
    [m, p] = sigmpoles(p);
    
    for indx = 1:length(m)
        if m(indx) ~= 1,
            h(end+1) = text(real(p(indx))+fuzz, imag(p(indx)), num2str(m(indx)), ...
                'Parent', hax, 'Vertical', 'Bottom', varargin{:});
        end
    end
end

% [EOF]
