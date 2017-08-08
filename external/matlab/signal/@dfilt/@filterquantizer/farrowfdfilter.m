function [y,z] = farrowfdfilter(this,C,x,d,z)
%FARROWFDFILTER   

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

x = quantizeinput(this,x);
d = quantizeinput(this,d);
nd = size(C,1);
nx = size(x,1);

if ~isempty(z)
    z = [z(1,:);z];
end

y = quantizeinput(this,zeros(size(x)));
for m=1:nx
    % ---------------------------------------------------------------------
    % Step 1: Bank of FIR Filters c(n,l+1) in harris
    % ---------------------------------------------------------------------
    % Load input
    z(1,:)=x(m,:);
    % Apply Differentiating Filters (they share the same states)
    p=C*z;
    % Update States
    z(2:end,:)=z(1:end-1,:);
    
    % ---------------------------------------------------------------------
    % Step 2: Interpolation based on the fractional delay d
    % ---------------------------------------------------------------------
    y(m,:) = d*p(1,:);
    for i=2:nd-1,
        y(m,:) = d*(p(i,:)+y(m,:));
    end
    y(m,:) = y(m,:)+p(nd,:);
end

z = z(2:end,:);



% [EOF]
