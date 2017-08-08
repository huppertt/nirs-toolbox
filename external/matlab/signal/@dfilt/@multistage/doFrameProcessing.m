function flag = doFrameProcessing(Hd)
%DOFRAMEPROCESSING Returns true if frame processing if supported by realizemdl()

%   Copyright 2010 The MathWorks, Inc.

nsections = length(Hd.Stage); 
flag = zeros(1,nsections);
for k=1:nsections, 
   flag(k) = doFrameProcessing(Hd.Stage(k));
end 
flag = all(flag);

% [EOF]
