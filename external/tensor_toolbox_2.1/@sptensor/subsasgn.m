function t = subsasgn(t,s,rhs)
%SUBSASGN Subscripted assignment for sparse tensor.
%
%   We can assign elements to a sptensor in three ways.
%
%   Case 1: X(R1,R2,...,RN) = Y, in which case we replace the
%   rectangular subtensor (or single element) specified by the ranges
%   R1,...,RN with Y. The right-hand-side can be a scalar or an
%   sptensor. 
%
%   Case 2: X(S) = V, where S is a p x n array of subscripts and V is
%   a scalar value or a vector containing p values.
%
%   Linear indexing is not supported for sparse tensors.
%
%   Examples
%   X = sptensor([30 40 20]) %<-- Create an emtpy 30 x 40 x 20 sptensor
%   X(30,40,20) = 7 %<-- Assign a single element to be 7
%   X([1,1,1;2,2,2]) = 1 %<-- Assign a list of elements to the same value
%   X(11:20,11:20,11:20) = sptenrand([10,10,10],10) %<-- subtensor!
%   X(31,41,21) = 7 %<-- grows the size of the tensor
%   X(111:120,111:120,111:120) = sptenrand([10,10,10],10) %<-- grows
%   X(1,1,1,1) = 4 %<-- increases the number of dimensions from 3 to 4
%
%   X = sptensor([30]) %<-- empty one-dimensional tensor
%   X([4:6]) = 1 %<-- set subtensor to ones (does not increase dimension)
%   X([10;12;14]) = (4:6)'  %<-- set three elements
%   X(31) = 7 %<-- grow the first dimension
%   X(1,1) = 0 %<-- add a dimension, but no nonzeros
%
%   See also SPTENSOR, TENSOR/SUBSASGN.
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
% $Id: subsasgn.m,v 1.13 2006/11/26 17:59:01 tgkolda Exp $


switch s.type

    case '.'
        error(['Cannot change field ', s.subs, ' directly.']);

    case '()'

        % Do nothing if both subscripts and RHS are empty
        if isempty(s.subs{1}) && isempty(t.vals)
            return;
        end

        % Figure out if we are doing a subtensor or a list of subscripts...
        type = 'error';       
        if ndims(t) == 1 
            % Multiple subscripts *or* row vector -> subtensor
            if (numel(s.subs) > 1) ...
                    || (ndims(s.subs{1}) == 2 && size(s.subs{1},1) == 1)
                type = 'subtensor';
            % Column vector *or* matrix of subscripts -> subscripts
            elseif (ndims(s.subs{1}) == 2)
                type = 'subscripts';
            end
        else
            % Multiple subscripts -> subtensor
            if numel(s.subs) >= ndims(t)
                type = 'subtensor';
            % Matrix of subscripts -> subscripts
            elseif (ndims(s.subs{1}) == 2) && (size(s.subs{1},2) >= ndims(t))
                type = 'subscripts';
            end
        end

        
        % Case I: Replace a sub-tensor
        if isequal(type,'subtensor')                     

            % Resize the tensor!
            for n = 1:ndims(t)
                if (s.subs{n} == ':')
                    newsz(1,n) = t.size(n);
                else
                    newsz(1,n) = max([t.size(n) s.subs{n}]);
                end
            end
            for n = ndims(t)+1:numel(s.subs)
                newsz(1,n) = max(s.subs{n});
            end
            t.size = newsz;

            % Resize existing subscripts if needed
            if ~isempty(t.subs) && (size(t.size,2) > size(t.subs,2))
                t.subs(:,end+1:size(t.size,2)) = 1;
            end

            if isa(rhs,'sptensor')

                % Delete what currently occupies the specified range
                rmloc = subdims(s.subs,t);
                kploc = setdiff(1:nnz(t),rmloc);
                newsubs = t.subs(kploc,:);
                newvals = t.vals(kploc);

                % Renumber the subscripts
                addsubs = irenumber(rhs, t.size, s.subs);
                t.subs = [newsubs; addsubs];
                t.vals = [newvals; rhs.vals];
                
            elseif numel(rhs) == 1 && rhs == 0

                % Delete what currently occupies the specified range
                rmloc = subdims(s.subs,t);
                kploc = setdiff(1:nnz(t),rmloc);
                t.subs = t.subs(kploc,:);
                t.vals = t.vals(kploc);
               
            elseif numel(rhs) == 1
                
                % Determine number of dimensions (may be larger than
                % current number)
                N = numel(s.subs); 
                
                % Figure out how many indices are in each dimension
                nssubs = zeros(N,1);
                for n = 1:N
                    if s.subs{n} == ':'
                        s.subs{n} = 1:size(t,n);
                    end
                    nssubs(n) = numel(s.subs{n});
                end
                
                % Preallocate (discover any memory issues here!)
                addsubs = zeros(prod(nssubs),N);

                % Generate appropriately sized ones vectors.
                o = cell(N,1);
                for n = 1:N
                    o{n} = ones(nssubs(n),1);
                end

                % Generate each column of the subscripts in turn
                for n = 1:N
                    i = o;
                    i{n} = s.subs{n}';
                    addsubs(:,n) = khatrirao(i);
                end
                
                if ~isempty(t.subs)
                    % replace existing values
                    [junk,loc] = intersect(t.subs,addsubs,'rows');
                    t.vals(loc) = rhs;
                    % pare down list of subscripts to add
                    addsubs = setdiff(addsubs,t.subs,'rows');
                end
                t.subs = [t.subs; addsubs];
                t.vals = [t.vals; rhs*ones(size(addsubs,1),1)];

            else
                error('Invalid RHS')
            end

            return;

        elseif isequal(type,'subscripts')

            % Case II: Replacing values at specified indices

            newsubs = [s.subs{1}];

            % Error check on subscripts
            if (ndims(newsubs) ~= 2) || (size(newsubs,2) < ndims(t))
                error('Invalid subscripts');
            end
            
            % Check for expanding the order
            if size(newsubs,2) > ndims(t)
                t.size(end+1:size(newsubs,2)) = 1;
                if ~isempty(t.subs)
                    t.subs(:,end+1:size(newsubs,2)) = 1;
                end
            end

            % Copy rhs to newvals
            newvals = rhs;

            % Error check the RHS is a column vector. We do not bother to
            % handle any other type of RHS with the sparse tensor.
            if (ndims(newvals) ~= 2) || (size(newvals,2) ~= 1)
                error('Invalid RHS --- must be a column vector');
            end

            % Determine number of nonzeros being inserted. (This is
            % determined by the number of subscripts. Later we will check
            % to see that it matches the size of the RHS.)
            newnnz = size(newsubs,1);

            % Error check on size of newvals
            if numel(newvals) == 1

                % Special case where newvals is a single element to be
                % assigned to multiple RHS. Fix to be correct size.
                newvals = newvals * ones(newnnz,1);

            elseif size(newvals,1) ~= newnnz

                % Sizes don't match!
                error('Number of subscripts and number of values do not match!');

            end

            % Remove duplicates & print warning if any duplicates were
            % removed.
            [newsubs,idx] = unique(newsubs,'rows');
            if size(newsubs,1) ~= newnnz
                warning('Duplicate assignments discarded.');
            end
            newvals = newvals(idx);

            % Find which subscripts already exist and their locations
            [tf,loc] = ismember(newsubs,t.subs,'rows');

            % Split into three groups for processing:
            %
            % Group A: Elements that already exist and need to be changed
            % Group B: Elements that already exist and need to be removed
            % Group C: Elements that do not exist and need to be added
            %
            % Note that we are ignoring any new zero elements, because
            % those obviously do not need to be added. Also, it's
            % important to process Group A before Group B because the
            % processing of Group B may change the locations of the
            % remaining elements.

            idxa = find((tf .* newvals) ~= 0);
            idxb = find((tf .* ~newvals) ~= 0);
            idxc = find((~tf .* newvals) ~= 0);

            % Process Group A: Changing values
            if ~isempty(idxa)
                t.vals(loc(idxa)) = newvals(idxa);
            end

            % Process Group B: Removing values
            if ~isempty(idxb)
                removesubs = loc(idxb);
                keepsubs = setdiff(1:nnz(t),removesubs);
                t.subs = t.subs(keepsubs,:);
                t.vals = t.vals(keepsubs);
            end

            % Process Group C: Adding new, nonzero values
            if ~isempty(idxc)
                t.subs = [t.subs; newsubs(idxc,:)];
                t.vals = [t.vals; newvals(idxc)];
            end

            % Resize the tensor!
            for n = 1:length(t.size)
                smax = max(newsubs(:,n));
                t.size(n) = max(t.size(n), smax);
            end

            return;

        end
        
        error('Invalid call to sptensor/subsasgn');

    case '{}'
        error('Subscript cell reference not supported for sptensor.');

    otherwise
        error('Incorrect indexing into sptensor.')

end


