function n = order(Hb)
%ORDER Filter order.
%   ORDER(Hb) returns the order of filter Hb.
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

n = base_num(reffilter(Hb), 'thisorder');
