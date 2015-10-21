function imQ = BenjaminiHochberg(p,FDR)
% Returns the Benjamini-Hochberg adjusted pvalues

s = size(p);   
[p, I] = sort(p(:));

% number of hypothesis tests
m = length(p);

imQ=[1:m]'/m*FDR;

imQ(I) = imQ;
imQ = reshape(imQ, s);

end