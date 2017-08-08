function Hnew = decimmergedelays(this,Hvec)
%DECIMMERGEDELAYS   

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

count = 1;
for k = 1:length(Hvec)-1,
    if strcmpi(class(Hvec(k)),'dfilt.delay') && ...
            strcmpi(class(Hvec(k+1)),'mfilt.firdecim'),
        Hnew(count) = mfilt.firdecim(Hvec(k+1).DecimationFactor,...
            [zeros(1,Hvec(k).Latency),Hvec(k+1).Numerator]); 
        count = count + 1;
    elseif  strcmpi(class(Hvec(k)),'mfilt.firdecim') && ...
            (k == 1 || ~strcmpi(class(Hvec(k-1)),'dfilt.delay'))
        Hnew(count) = Hvec(k);
        count = count + 1;
    end    
end
Hnew(count) = Hvec(end);

% [EOF]
