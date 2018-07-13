function d=fixnan(d)

ifFailReplaceWith = 1;

lst = isnan(d);

if any(lst(:))
    try
        for j = 1:size(d,2)
            if(all(lst(:,j)))
                d(lst(:,j),j) =ifFailReplaceWith;
            elseif any(lst(:,j))
                % interpolation
                l = lst(:,j);
                d(l,j) = interp1(t(~l), d(~l,j), t(l),'linear','extrap');
                data(i).data = d;
            end
        end
        
    catch
        % just replace with white noise
        d(lst) = ifFailReplaceWith;
    end
    
end