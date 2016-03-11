function disp(X,name)
%DISP Command window display of a sparse tensor.
%
%   DISP(X) displays the tensor without printing its name.
%
%   DISP(X,NAME) displays the tensor with the given name.
%
%   See also SPTENSOR, SPTENSOR/DISPLAY.
%
%MATLAB Tensor Toolbox.
%Copyright 2006, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2006) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: disp.m,v 1.15 2006/10/24 16:12:44 tgkolda Exp $

% Extract the number of nonzeros and number of dimensions
nz = nnz(X);
n = ndims(X);

if ~exist('name','var')
    name = 'ans';
end

if (nz == 0)
    fprintf('%s is an all-zero sparse tensor of size %s\n',...
        name, tt_size2str(X.size));
    return;
else
    fprintf('%s is a sparse tensor of size %s with %d nonzeros\n',...
        name, tt_size2str(X.size), nz);
end

% Stop insane printouts
if (nz > 10000)
    r = input('Are you sure you want to print all nonzeros? (Y/N) ','s');
    if upper(r) ~= 'Y', return, end;
end

% preallocate
output = cell(nz,1);
%%
spc = floor(log10(max(X.subs,[],1)))+1;
if numel(spc) == 1
    fmt = ['\t(%' num2str(spc(1)) 'd)\t%g'];
else
    fmt = ['\t(%' num2str(spc(1)) 'd,'];
    for i = 2:numel(spc)-1
        fmt = [fmt '%' num2str(spc(i)) 'd,'];
    end
    fmt = [fmt '%' num2str(spc(i)) 'd)\t%g'];
end
%%      
for i = 1:nz
    output{i} = sprintf(fmt,X.subs(i,:),X.vals(i));
end
fprintf('%s\n',output{:});

%     function y = fmt(s,v)
%         % nested function has access to n from disp function workspace
%         if n > 1
%             y = [sprintf('\t(') sprintf('%d,',s(1:n-1))...
%                 sprintf('%d) ',s(n)) sprintf('\t%g',v)];
%         else
%             y = sprintf('\t(%d) \t%f',s,v);
%         end
%     end

end
