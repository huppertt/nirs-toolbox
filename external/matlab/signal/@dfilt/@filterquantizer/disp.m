function disp(this, spacing)
%DISP Object display.
  
%   Author: R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

if ~isempty(fieldnames(this)),
    siguddutils('dispstr', this, spacing);
end
