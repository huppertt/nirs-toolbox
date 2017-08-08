function str = getinfoheader(h)
%GETINFOHEADER Return the header to the info

%   Copyright 1988-2013 The MathWorks, Inc.

% Setup the title string
if isreal(h),
  typestr = 'real';
else
  typestr = 'complex';
end

if isfromdesignfilt(h)
  if isfir(h),
    str = ['FIR Digital Filter (', typestr,')'];
  else
    str = ['IIR Digital Filter (', typestr,')'];
  end
  
else
  if isfir(h),
    typestr = ['FIR Filter (', typestr,')'];
  else
    typestr = ['IIR Filter (', typestr,')'];
  end
  str = ['Discrete-Time ',typestr];  
end

% [EOF]
