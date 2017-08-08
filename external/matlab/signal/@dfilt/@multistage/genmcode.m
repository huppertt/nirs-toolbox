function str = genmcode(Hm, objname, place)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 2, objname = 'Hd'; end
if nargin < 3, place   = [];   end

subname = sprintf('%ss', objname);

str = {sprintf('%% Create the 1st stage%s', getstage(place)), ...
        getstagestr(Hm, 1, subname, [place 1]), ...
        sprintf('%s = %s(%s);', objname, class(Hm), subname)};

for indx = 2:length(Hm.Stage)
    
    str = {str{:}, ...
            '', ...
            sprintf('%% Add the %s stage%s', genmcodeutils('formatnumplace', indx),...
                getstage(place)), ...
            getstagestr(Hm, indx, subname, [place indx]), ...
            sprintf('addstage(%s, %s);', objname, subname)};
end

str = sprintf('%s\n', str{:});
str(end) = [];

% ----------------------------------------------------------------------
function s = getstage(place)

s = '';
for indx = length(place):-1:1
    s = sprintf('%s of the %s stage', s, ...
        genmcodeutils('formatnumplace', place(indx)));
end

if length(s) > 70,
    s(65:end) = [];
    s = [s ' ...'];
end
% -------------------------------------------------------------------------
function str = getstagestr(Hm, indx, subname, place)

str = genmcode(Hm.Stage(indx), subname, place);

if isa(str, 'sigcodegen.mcodebuffer')
    str = string(str);
end

% [EOF]
