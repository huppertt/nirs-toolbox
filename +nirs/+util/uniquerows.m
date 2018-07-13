function [uvars,uind,uind2] = uniquerows( vars )
% Equivalent to [C,IX,IC] = unique( X , 'rows', 'stable') but can handle cell arrays
try
    [uvars,uind,uind2] = unique( vars , 'rows' , 'stable' );
catch
    N = height(vars);
    hashes = cell(N,1);
    opt.Method = 'SHA-512';
    opt.Format = 'hex';
    for i = 1:N
        hashes{i} = DataHash(vars(i,:),opt);
    end
    [~,uind,uind2] = unique( hashes , 'stable' );
    uvars = vars(uind,:);
end
end