function result = parallel_init()

% parallel_init()
%
% Initializes the parallel workers if the toolbox is available
%
% result is 1 if succesful, 0 otherwise

if (exist('matlabpool')) % Start labs if necessary. 
    sz = matlabpool('size'); 
    if (sz ==0) 
        matlabpool('open'); 
    end
    % Check we got some now. 
    sz = matlabpool('size'); 
    if (sz ==0) 
        error('Failed to open parallel workers'); 
        result = 0;
    else
        fprintf('Running on %d workers\n', sz); 
        result = 1;
    end
else
    result = 0;
end