function Hd2 = convert(Hd,newstruct)
%CONVERT Convert structure of DFILT object.
%   CONVERT(H,NEWSTRUCT) converts DFILT object H to the structure defined by
%   string NEWSTRUCT.
%
%   EXAMPLE:
%           Hd1 = dfilt.df2t;
%           Hd2 = convert(Hd1,'df1');
%           % returns Hd2 as a direct-form 1 discrete-time filter.
%  
%   See also DFILT.

%   Copyright 1988-2008 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'))

if ~ischar(newstruct)
  error(message('signal:dfilt:singleton:convert:MustBeAString'))
end

if ~ismethod(Hd,['to',newstruct])
  error(message('signal:dfilt:singleton:convert:NotSupported', newstruct));
end  

p = findprop(Hd, 'arithmetic');

% If the property is not there or is private, we are in "double" mode.
if ~isempty(p) && strcmpi(p.AccessFlags.PublicGet, 'On')
    switch lower(Hd.Arithmetic)
        case 'single'
            warning(message('signal:dfilt:singleton:convert:unquantizingSingle'));
            Hd = reffilter(Hd);
        case 'fixed'
            warning(message('signal:dfilt:singleton:convert:unquantizingFixed'));
            Hd = reffilter(Hd);
        otherwise
            % NO OP, double.
    end
end

Hd2 = feval(['to',newstruct], Hd);
setfdesign(Hd2,getfdesign(Hd)); % Carry over fdesign obj
setfmethod(Hd2,getfmethod(Hd));

if ~isa(Hd2,'dfilt.singleton')
  error(message('signal:dfilt:singleton:convert:DFILTErr'));
end
