function q = BenjaminiHochberg(p)
% Returns the Benjamini-Hochberg adjusted pvalues
porig=p;
s = size(p);   
[p, I] = sort(p,'ascend');

% number of hypothesis tests
m = length(p);

imQ=p.*m./[1:m]';

for i=length(imQ):-1:2; if(imQ(i)<imQ(i-1)); imQ(i-1)=imQ(i); end; end;

imQ=min(imQ,1);

q=ones(s);
q(I) = imQ;
q = reshape(q, s);

end