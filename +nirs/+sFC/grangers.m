function [G,F,dfe1,dfe2,p]=grangers(data)
% This function runs the Granger's causality models

%TODO- write this


[G, F, df1, df2, p] = nirs.math.mvgc(Y, Pmax,includeZeroLag);

[G, F,dfe1,dfe2, P] = nirs.math.robust_mvgc(data(i).data,pMax,includeZeroLag);
