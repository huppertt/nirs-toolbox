function newsubs = irenumber(t, sz, range)
%RENUMBER indices for sptensor subsasgn
%
%  See also SPTENSOR/SUBSASGN

nz = nnz(t);
if (nz == 0)
    newsubs = [];
    return;
end

newsubs = t.subs;
for i = 1 : numel(range)
    r = range{i};
    if r == ':'
        continue;
    elseif numel(r) == 1
        newsubs = [newsubs(:,1:i-1), r*ones(nz,1), newsubs(:,i:end)];
    else
        newsubs(:,i) = r(newsubs(:,i));
    end
end

