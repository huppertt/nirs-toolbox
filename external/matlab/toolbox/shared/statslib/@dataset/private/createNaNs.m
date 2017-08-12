function y = createNaNs(sz,x)
% Create NaN(sz,class(x)), even when x is not built-in
%
% Assumes floating point, and that assigning NaN works
if isobject(x)
    y = x([]); y(1) = NaN; y = repmat(y,sz);
else
    y = NaN(sz,class(x));
end
