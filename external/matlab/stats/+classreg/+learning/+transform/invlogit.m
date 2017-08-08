function out = invlogit(in)

%   Copyright 2010 The MathWorks, Inc.


is0 = in==0;
is1 = in==1;
out = zeros(size(in));
out(is0) = -Inf;
out(is1) = +Inf;
todo = ~is0 & ~is1;
out(todo) = log(in(todo)./(1-in(todo)));
end
