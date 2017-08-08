function displayMethodsProperties(obj)
%DISPLAYMETHODSPROPERTIES Display links to object methods and properties.
%   DISPLAYMETHODSPROPERTIES(OBJ) displays links to object methods and
%   properties. 

%   Copyright 2012 The MathWorks, Inc.

mc = metaclass(obj);
hotlinks = feature('hotlinks');

if hotlinks
    fprintf('\n  <a href="matlab: properties(''%s'')">%s</a>, ', ...
        mc.Name,'Properties');
    fprintf('<a href="matlab: methods(''%s'')">%s</a>\n\n', ...
        mc.Name,'Methods');
else
    fprintf('\n');
end

end
