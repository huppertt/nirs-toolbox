function M = mode(x)

% M = mode(x)
%
% computes the mathematical mode of x
%
% x is the input vector
% M is the resulting mode


sp_flag = issparse(x);
sx = size(x);
if 1>length(sx)
    sx = [sx, ones(1,1-length(sx))];
end

sizem = sx;
sizem(1) = 1;

% preallocate
if sp_flag
    M = sparse(sizem(1),sizem(2));
else
    M = zeros(sizem,class(x));
end

% remove empty arrays
if isempty(x)
    if docell
        C(:) = {M(1:0)};
    end
    if prod(sizem)>0
        M(:) = NaN;
        if dofreq
            F(:) = 0;
        end
    end
    return
end

% Convert data to operate along columns of a 2-d array
x = permute(x,[1, (1:1-1), (1+1:length(sx))]);
x = reshape(x,[sx(1),prod(sizem)]);
[nrows,ncols] = size(x);

for j=1:ncols
    v = sort(x(:,j));                       
    start = find([1; v(1:end-1)~=v(2:end)]);
    freq = [start(2:end);nrows+1] - start;
    [maxfreq,firstloc] = max(freq);
    M(j) = v(start(firstloc));
end
