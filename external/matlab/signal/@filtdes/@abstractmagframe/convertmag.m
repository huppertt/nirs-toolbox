function convertmag(h,d,oldvalsStrs,newvalsStrs,passOrStop,convertto)
%CONVERTMAGFIR Convert mag specs.
%
%   Inputs:
%       d          - handle to design method object
%       oldvalStrs - cell of strings of old value props
%       newvalStrs - cell of strings of new value props
%       passOrStop - cell of strings indicating if value 
%                    is passband or stopband
%       convertto  - function handle to conversion method
%
% This should be a private method.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Get old values
oldvals = get(d,oldvalsStrs);
      
for n = 1:length(newvalsStrs),
    newvals(n) = feval(convertto,h,oldvals{n},passOrStop{n});
    newdatatype{n} = 'udouble';
    descript{n} = 'magspec';
end

% Call the modify properties method of the design object
modifyProps(d,oldvalsStrs,newvalsStrs,newvals,newdatatype,descript);



