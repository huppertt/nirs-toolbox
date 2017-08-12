function a = readXPTFile(~,xptfile,args)
%READFILE Read in a SAS export file and create a dataset array.

%   Copyright 2011 The MathWorks, Inc.


pnames = {'readvarnames' 'readobsnames'};
dflts =  {          false          false};
[readvarnames,readobsnames] ...
    = dataset.parseArgs(pnames, dflts, args{:});
readobsnames = onOff2Logical(readobsnames,'ReadObsNames');
readvarnames = onOff2Logical(readvarnames,'ReadVarNames');

% Var names are always present in the XPT file, and are always read, and are
% always valid (max 8 char alphanumeric).
if(readvarnames)
   warning(message('stats:dataset:dataset:XPTReadVarNotSupported'));
end

a = xptread(xptfile,'ReadObsNames',readobsnames);
a = table2dataset(a);
a = set(a,'DimNames',{'Observations' 'Variables'});