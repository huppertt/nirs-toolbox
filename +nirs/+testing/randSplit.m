function idx = randSplit(N, ngroups)
%% randSplit - Returns the groupings to randomly split N items into ngroups
%
% Args: 
%     N       - number of items
%     ngroups - number of groups
%     
% Returns:
%     idx - vector with entries between 1 to N indicating the group identity
%
% Notes:
%     Splits into groups as evenly as possible (i.e. 10 to groups of 3/3/4)
    
    p = randperm(N)';
    idx = zeros(size(p));
    
    r = randi(ngroups);
    for i = 1:ngroups
        lst = mod(i+r,ngroups)+1:ngroups:N;
        idx(p(lst)) = i;
    end
end

