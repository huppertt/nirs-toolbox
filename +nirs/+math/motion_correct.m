function ymoco=motion_correct(y,Pmax)

n = length(y);

Xf = nirs.math.lagmatrix(y, 1:Pmax);
Xb = nirs.math.lagmatrix(flipud(y), 1:Pmax);

X = [ones(2*n,1) [Xf; Xb]];
[coef, res] = nirs.math.stepwise(X, [y; flipud(y)]);
    
f=[1; -coef(2:end)];

c=convmtx(f,length(y));
c=c(1:length(y),:);
yf = c*y;

%e = y-pinv(c)*yf;

r=yf;
s = mad(r, 0) / 0.6745;
r = r/s/4.685;
w = (1 - r.^2) .* (r < 1 & r > -1);
ymoco=inv(c)*(w.*yf);


%ymoco=pinv(c)*(w.*yf)+e;