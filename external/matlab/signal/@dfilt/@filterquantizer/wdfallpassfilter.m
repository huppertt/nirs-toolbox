function [y,zf] = wdfallpassfilter(this,c,x,zi)
%WDFALLPASSFILTER  
%
%   Reference: M. Lutovac, D. Tosic, B. Evans. Filter Design for Signal
%   Processing using MATLAB and Mathematica. Prentice Hall, 2001.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.


if isempty(c),
    % Special case a wire
    y = x;
    zf = [];
else
    % Preallocate output
    y = zeros(size(x));

    % Preallocate cell to store function handles
    fh = cell(size(c));

    for n = 1:length(c),
        if c(n) > 0 && c(n) <= 0.5,
            fh{n} = @wdfallpassfilter1a;
        elseif c(n) >= -1 && c(n) < -0.5,
            fh{n} = @wdfallpassfilter1b;
        elseif c(n) <= 1 && c(n) > 0.5,
            fh{n} = @wdfallpassfilter1c;
        elseif c(n) < 0 && c(n) >= -0.5,
            fh{n} = @wdfallpassfilter1d;
        elseif c(n) == 0,
            fh{n} = @wdfallpassfilter1e;
        end
    end

    xb = zi;
    for k = 1:size(x,1),
        % Filter in reverse order so the output of one stage is the input of
        % the previous
        xv = [x(k,:);xb];
        for n = length(c)+1:-1:2,
            [xv(n-1,:),xv(n,:)] = feval(fh{n-1},c(n-1),xv(n-1,:),xv(n,:));
        end
        y(k,:) = xv(1,:);
        xb = xv(2:end,:);
    end
    zf = xb;
end


%----------------------------------------------------
function [y,yb] = wdfallpassfilter1a(c,x,xb)
% Case c > 0 && c <= 0.5; Figure 8.51 in reference

a = c;

v  = x - xb;
y  = xb + a*v;
yb = y + v;
%----------------------------------------------------
function [y,yb] = wdfallpassfilter1b(c,x,xb)
% Case c >= -1 && c < -0.5; Figure 8.52 in reference

a = 1 + c;

v  = x - xb;
yb = xb + a*v;
y  = yb - v;

%----------------------------------------------------
function [y,yb] = wdfallpassfilter1c(c,x,xb)
% Case c <= 1 && c > 0.5; Figure 8.48 in reference

a = 1 - c;

v  = x + xb;
yb = a*v - xb;
y  = v - yb;

%----------------------------------------------------
function [y,yb] = wdfallpassfilter1d(c,x,xb)
% Case c < 0 && c >= -0.5; Figure 8.49 in reference

a = -c;

v  = xb - x;
y  = a*v + xb;
yb = y - v;
%----------------------------------------------------
function [y,yb] = wdfallpassfilter1e(c,x,xb)
% Case c = 0  Figure 8.50 in reference


y  = xb;
yb = x;



% [EOF]
