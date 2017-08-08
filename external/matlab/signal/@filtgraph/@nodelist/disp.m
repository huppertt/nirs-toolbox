function disp(O)
%DISP Object display.
  
%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

N = length(O);
msg = sprintf( ' %s \n Length: %d.\n',class(O),N);
disp(msg);