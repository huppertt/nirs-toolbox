function q = BenjaminiHochberg(p)
% Returns the Benjamini-Hochberg adjusted pvalues

s = size(p);   
[p, I] = sort(p,'descend');

% number of hypothesis tests
m = length(p);

imQ=p.*m./[1:m]';
imQ=min(imQ,1);

q=ones(s);
q(I) = imQ;
q = reshape(q, s);

end