function h = decimdelaynoble(this,Hd,dfactor,delay,h,count)
%DECIMDELAYNOBLE   

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

if nargin < 4,
    delay = Hd.Latency;
end

if nargin < 6,
    count = 1;
end

d1 = rem(delay,dfactor);

d2 = delay - d1;

% Use noble identity to compute rightmost delay
d2 = d2/dfactor;

if d2 ~= 0,
    h(count) = dfilt.delay(d2);
    count = count + 1;
end

if d1 == 0,
    % We are done
    h(count) = mfilt.firdecim(dfactor,1);
else
    % Fact: d1 < dfactor
    f = factor(dfactor);
    if isprime(dfactor) || all(gcd(f,d1)==1),
        % We are done
        h(count) = mfilt.firdecim(dfactor,[zeros(1,d1),1]);
    else
        idx = 1;
        created_filt = false;
        while length(f(1:end-idx)) >= 1,
            pr = prod(f(1:end-idx));
            if pr <= d1,
                created_filt = true;
                reddfactor = dfactor/pr;
                h(count) = mfilt.firdecim(reddfactor,1);
                count = count + 1;
                h = decimdelaynoble(this,Hd,pr,d1,h,count); 
                break;
            end
            idx = idx + 1;                
        end
        if ~created_filt,
            % Couldn't find further simplification
              h(count) = mfilt.firdecim(dfactor,[zeros(1,d1),1]);
        end
    end
end


% [EOF]
