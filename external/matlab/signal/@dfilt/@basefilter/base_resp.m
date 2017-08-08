function [h, w] = base_resp(Hb, fcn, varargin)
%BASE_RESP Return the response specified
%   BASE_RESP(Hd, FCN) Returns the response specified by the function
%   handle FCN for the filters Hd.

%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% This should be private

h = [];

for indx = 1:length(Hb),
    Hd = dispatch(Hb(indx));
    
    for jndx = 1:length(Hd),
        [ht, w] = feval(fcn, Hd(jndx), varargin{:});
        
        if isempty(h),
            h = ht;
        else
            
            % If the number of columns == 1, expand in the columns, otherwise
            % expand in the rows
            if size(ht, 2) == 1, h(:, end+1) = ht;
            else,                h(end+1, :) = ht; end
        end
    end
end

% [EOF]
