function functionNameClashCheck(optionName,optionValue,internalName,exceptionID)
%

% functionNameClashCheck Helper function that checks for name clash.
%
% If user sets option optionName (string) to an optionValue (function handle or
% string) that is the name of an internal helper function, this helper function 
% throws, as caller, an error with ID exceptionID (string) about a possible name 
% clash.

%   Copyright 2009-2011 The MathWorks, Inc.

% Flag that indicates if user-supplied function is of a supported type:
% function handle or a string
if isa(optionValue,'function_handle')
    validFcnInput = true; 
    optionValueString = func2str(optionValue); % needed for message    
elseif ischar(optionValue)
    validFcnInput = true; 
    optionValueString = optionValue;           % needed for message 
else
    validFcnInput = false; 
end
if validFcnInput && strcmpi(optionValueString,internalName)
    ME = MException(exceptionID, ...
        getString(message('optimlib:commonMsgs:FcnNameClash',optionValueString,optionName)));
    throwAsCaller(ME);
end

        