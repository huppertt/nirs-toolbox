function modifyProps(d,oldvalsStrs,newvalsStrs,newvals,newdatatype,descript)
%MODIFYPROPS Delete specified old properties and create new ones.
%
%   Note: Properties to be modified cannot have listeners
%
%   Inputs:
%       d - Container object
%       oldvalsStrs - cell of strings of old value props
%       newvalsStrs - cell of strings of new value props
%       newvals     - array of new values
%       newdatatype - cell of data types for new props
%       descript    - cell of descriptions for new props

%   Author(s): R. Losada
%   Copyright 1988-2009 The MathWorks, Inc.

% Find existing dynamic properties
p = get(d,'dynamicProps');

% Create new dynamic properties
% do not put the values in until later so that we avoid the problem of
% setting an invalid value into the first prop, erroring and then the
% second prop doesn exist after that
for n = 1:length(newvals)
    pnew(n) = schema.prop(d,newvalsStrs{n},newdatatype{n});
    set(d,newvalsStrs{n});
    pnew(n).Description = descript{n};
    
    % Retain handles to properties that will remain
    p = find(p,'-not','name',oldvalsStrs{n});
    
    % Delete old properties
    delete(findprop(d,oldvalsStrs{n}));
end

% Contain new and remaining dynamic props
p = [p;pnew(:)];
set(d,'dynamicProps',p);

%now put the new values into the dynamic properties
for n = 1:length(newvals)
    %suppress command line error generation.  suppressing the command line
    %error because I have already alerted the user to the errors with
    %better verbaged dialogs.  no need for both dialogs and command line
    %errors
     try
        set(d,newvalsStrs{n},newvals(n));
     catch  %#ok<CTCH>
     end
end
