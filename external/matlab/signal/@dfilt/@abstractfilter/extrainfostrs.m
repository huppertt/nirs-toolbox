function strs = extrainfostrs(this)
%EXTRAINFOSTRS   Return the extra info strings.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if isquantized(this)
    str = {evalc('disp(this.filterquantizer)')};
    indx = 1;
    while indx <= length(str),
        
        % Find the first newline feed
        idx = min(strfind(str{indx}, char(10)));
        if isempty(idx),
            indx = indx + 1;
        else

            % If there is a newline feed convert it into two entries of a cell array.
            str = {str{1:indx-1} str{indx}(1:idx-1) [] str{indx}(idx+1:end) str{indx+1:end}};
            indx = indx + 1;
        end
    end
    strs = {sprintf('\n'), str{:}};

else
    strs = {};
end

% [EOF]
