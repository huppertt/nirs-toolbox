function [y,zf] = allpassfilter(this,a,x,zi)
%ALLPASSFILTER

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.




if isempty(a),
    % Special case a wire
    y = x;
    zf = [];
else
    % Preallocate output
    y = zeros(size(x));

    % Breakup the states in top and bottom
    ztop    = [zi(length(a):-1:1,:)];
    zbottom = zi(length(a)+1:end,:);

    for k = 1:size(x,1),

        % Isolate final top state
        ztf = ztop(1,:);

        ztop    = [ztop(2:end,:);x(k,:)];

        % Compute output
        y(k,:) = ztf + a*(ztop-zbottom);

        % Update states
        zbottom = [y(k,:);zbottom(1:end-1,:)];
    end

    % Reassemble states
    zf = [ztop(end:-1:1,:);zbottom];
end



    % [EOF]
