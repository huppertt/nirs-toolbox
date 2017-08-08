function hF = createhdlfilter(this)
%CREATEHDLFILTER Returns the corresponding hdlfiltercomp for HDL Code
%generation.

%   Copyright 2007 The MathWorks, Inc.

error(message('signal:dfilt:basefilter:createhdlfilter:NotHdlable', class( this )));
                   

% [EOF]
