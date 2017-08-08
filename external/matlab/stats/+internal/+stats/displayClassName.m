function displayClassName(obj,displayPackage)
%DISPLAYCLASSNAME Display class name.
%   DISPLAYCLASSNAME(OBJ) displays the class name for object OBJ.
%
%   DISPLAYCLASSNAME(OBJ,true) also displays the package name.

%   Copyright 2012 The MathWorks, Inc.

if nargin<2
    displayPackage = false;
end

mc = metaclass(obj);
hotlinks = feature('hotlinks');

if hotlinks
    fprintf('  <a href="matlab: helpPopup %s">%s</a>\n', ...
        mc.Name, mc.Name);
else
    fprintf('  %s\n', mc.Name);
end

if displayPackage
    if ~isempty(mc.ContainingPackage)
        fprintf('  %s: %s\n\n','Package', mc.ContainingPackage.Name);
    else
        fprintf('\n');
    end
end

end
