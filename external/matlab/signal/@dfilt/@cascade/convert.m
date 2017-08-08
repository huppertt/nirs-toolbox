function Hd2 = convert(Hd,newstruct)
%CONVERT Convert structure of DFILT object.
%  Hd2 = CONVERT(Hd,NEWSTRUCT) converts all sections of multisection DFILT
%  object Hd to the structure defined by string NEWSTRUCT.
%
%  EXAMPLE:
%          Hd1 = dfilt.cascade(dfilt.df1,dfilt.df2);
%          Hd2 = convert(Hd1,'df1')
%          % returns Hd2 with all sections as direct-form 1 discrete-time filters.
%  
%   See also DFILT.

%   Copyright 1988-2008 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'))

Hd2 = convert(Hd.Stage(1),newstruct);
for k=2:nstages(Hd)
  Hd2 = cascade(Hd2, convert(Hd.Stage(k),newstruct));
end
