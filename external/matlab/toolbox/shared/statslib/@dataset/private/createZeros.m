function y = createZeros(sz,x)
% Create zeros(sz,class(x)), even when x is not built-in
%
% Assumes numeric, and that assigning 0 works
if isobject(x)
    y = x([]); y(1) = 0; y = repmat(y,sz);
else
    y = zeros(sz,class(x));
end
