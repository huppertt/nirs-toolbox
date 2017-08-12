function [fieldValue,ME] = getNumericOrStringFieldValue(fieldName, ...
    allowedString,equivNumValue,dataTypeWithArticleForErrMsg,options,defaultopt)
%

%getNumericOrStringFieldValue reads in field value from options structure.
% This function processes option fields that can take either a numeric or
% a special string value.
% It reads in value of field fieldValue. If value is a string, it compares 
% it to allowedString and converts it to equivalent numeric value. If 
% string is not allowedString, it returns non-empty MATLAB exception
% object.
% Note: the input dataTypeWithArticleForErrMsg is the data type of the
% field along with its article for error printing purposes. E.g.
% 'an integer' or 'a real'.

ME = []; % return empty MATLAB exception object unless error occurs
% Read in value of field
fieldValue = optimget(options,fieldName,defaultopt,'fast');
if ischar(fieldValue)
    if strcmpi(fieldValue,allowedString)
        fieldValue = equivNumValue;
    else
        errId = ['optimlib:getNumericOrStringFieldValue:Invalid' fieldName 'Option'];
        ME = MException(errId,getString(message('optimlib:commonMsgs:InvalidOptionValue', ...
            fieldName,dataTypeWithArticleForErrMsg)));
    end
end
