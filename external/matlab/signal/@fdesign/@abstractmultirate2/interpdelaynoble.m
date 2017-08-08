function h = interpdelaynoble(this,Hd,ifactor,delay,h,count)
%INTERPDELAYNOBLE   

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

if nargin < 4,
    delay = Hd.Latency;
end

if nargin < 6,
    count = 1;
end

d1 = rem(delay,ifactor);

d2 = delay - d1;

% Use noble identity to compute rightmost delay
d2 = d2/ifactor;

if d2 ~= 0,
    h(count) = dfilt.delay(d2);
    count = count + 1;
end

if d1 == 0,
    % We are done
    h(count) = mfilt.firinterp(ifactor,ifactor); % Numerator is a scalar with gain ifactor
else
    % Fact: d1 < ifactor
    f = factor(ifactor);
    if isprime(ifactor) || all(gcd(f,d1)==1),
        % We are done
        h(count) = mfilt.firinterp(ifactor,[zeros(1,d1),ifactor]);
    else
        idx = 1;
        created_filt = false;
        while length(f(1:end-idx)) >= 1,
            pr = prod(f(1:end-idx));            
            if pr <= d1,
                created_filt = true;
                redifactor = ifactor/pr;
                h(count) = mfilt.firinterp(redifactor,redifactor);
                count = count + 1;
                h = interpdelaynoble(this,Hd,pr,d1,h,count); 
                break;
            end
            idx = idx + 1;                
        end
        if ~created_filt,
            % Couldn't find further simplification
             h(count) = mfilt.firinterp(ifactor,[zeros(1,d1),ifactor]);
        end
    end
end

% [EOF]
