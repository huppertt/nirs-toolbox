function loss = mse(Y,Yfit,W)

%   Copyright 2010 The MathWorks, Inc.


notNaN = ~isnan(Yfit);
loss = sum(W(notNaN) .* (Y(notNaN)-Yfit(notNaN)).^2) / sum(W(notNaN));
end
