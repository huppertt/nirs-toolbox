function b = base_isspecmet(this, Hd, varargin)
%BASE_ISSPECMET   Returns true if the specification is met.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

m = getmeasurements(this, Hd, Hd.getfmethod);

b = true;

for indx = 1:length(varargin)
    if isprop(this, varargin{indx}{1}) && ~isempty(m.(varargin{indx}{1}))
        switch varargin{indx}{2}
            case '>'
                if m.(varargin{indx}{1}) < this.(varargin{indx}{1})
                    b = false;
                end
            case '<'
                if m.(varargin{indx}{1}) > this.(varargin{indx}{1})
                    b = false;
                end
            case {'=', '=='}
                if m.(varargin{indx}{1}) ~= this.(varargin{indx}{1})
                    b = false;
                end
        end
    end
end

% [EOF]
