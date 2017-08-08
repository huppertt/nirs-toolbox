function Hnew = interpmergedelays(this,Hvec)
%INTERPMERGEDELAYS   

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

Hnew(1) = Hvec(1);
count = 2;
for k = 2:length(Hvec),
    if strcmpi(class(Hvec(k)),'dfilt.delay') && ...
            strcmpi(class(Hvec(k-1)),'mfilt.firinterp'),
        Hnew(count) = mfilt.firinterp(Hvec(k-1).InterpolationFactor,...
            [zeros(1,Hvec(k).Latency),Hvec(k-1).Numerator]); 
        count = count + 1;
    elseif  strcmpi(class(Hvec(k)),'mfilt.firinterp') && ...
            (k == length(Hvec) || ~strcmpi(class(Hvec(k+1)),'dfilt.delay'))
        Hnew(count) = Hvec(k);
        count = count + 1;
    end    
end


% [EOF]
