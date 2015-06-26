% this function generates new names based on a linear tranformation
% T of the variables
function newNames = transformNames( obj, T )
    names = obj.conditions;
    for i = 1:size(T,1)
        newNames{i} = '';
        for j = 1:size(T,2)
            c = T(i,j);
            if c == 1
                newNames{i} = [newNames{i} '+' names{j}];
            elseif c == -1
                newNames{i} = [newNames{i} '-' names{j}];
            elseif c > 0
                newNames{i} = [newNames{i} '+' num2str(c) names{j}];
            elseif c < 0
                newNames{i} = [newNames{i} num2str(c) names{j}];
            end
        end

        if newNames{i}(1) == '+'
            newNames{i} = newNames{i}(2:end);
        end
    end

    newNames = newNames(:);
end