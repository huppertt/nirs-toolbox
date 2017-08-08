function tstr = getTranslatedString(scopestr,str)
%This function is for internal use only. It may be removed in a future
%release.

%GETTRANSLATEDSTRING Get translated string
%   TSTR = getTranslatedString(STR)

%   Copyright 2011 The MathWorks, Inc.

try
    %    tstr = getString(message(sprintf('signal:siggui:renderedlabel:%s',...
    tstr = getString(message(sprintf('%s:%s',scopestr, ...
        regexprep(str,'(\W+)',''))));
catch ME %#ok<NASGU>
    %ADDCATALOG many ui are using this function to be translated. 
    % see for example matlab/toolbox/signal/sigtools/@siggui/@siggui/rendercontrols.m
    % xlate was removed so labels which previously were tranlated are not anymore.
    % This needs to be fixed.
    tstr = str;
    %fprintf('\n<StringToTranslate>%s<sep>%s</StringToTranslate>\n',scopestr,tstr);
end


% [EOF]
