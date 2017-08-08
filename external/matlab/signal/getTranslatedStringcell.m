function str = getTranslatedStringcell(scopestr,str)
%This function is for internal use only. It may be removed in a future
%release.

%GETTRANSLATEDSTRING Get translated string
%   TSTR = getTranslatedString(STR)

%   Copyright 2011 The MathWorks, Inc.

if iscell(str)
    for i = 1:length(str)
        subString = str{i};
        if iscell(subString)
            for ii = 1:length(subString)
                if ischar(subString{ii})
                str{i}{ii} =  ...
                    getTranslatedString(scopestr,  subString{ii});
                end
            end
        else
            if ischar(subString)
            str{i} = ...
                getTranslatedString(scopestr,  subString);
            end
        end
    end
else
    str = getTranslatedString(scopestr,  str);
end

 % [EOF]
