function [errid,errmsg] = dct_validate_inputs(myfun, size_x, class_x, isreal_x, N_in)
%Embedded MATLAB Library Function

% Validate inputs to toolbox/signal/eml/dct and idct.  

% Copyright 2009 The MathWorks, Inc.
errid  = '';  %#ok
errmsg = '';  %#ok

%%
% Catch any errors that the builtin function would have thrown
n_outputs = 1;
n_variable_inputs = 1;
[errid,errmsg] = eml_validate_toolbox_inputs(...
    myfun, n_outputs, n_variable_inputs, size_x, class_x, isreal_x);
if ~isempty(errmsg)
    return
end
%%
% Check limitations specific to the Embedded MATLAB version.
if any(size_x==0)
    errid = ['signal:',myfun,':emptyInput'];
    errmsg = 'Input must not be empty.';
    return
end

if length(size_x)~=2
    errid  = ['signal:',myfun,':nDimensionalInput'];
    errmsg = 'Input must be a vector or 2-dimensional matrix.';
    return
end

if nargin>4
    N = N_in;
else
    if size_x(1) ==1
        N = size_x(2);
    else
        N = size_x(1);
    end
end
r = log2(N);
if r ~= floor(r)
    errid  = ['signal:',myfun,':sizeMustBePower2'];
    errmsg = 'Length of transform dimension must be a power of 2.';
    return
end
