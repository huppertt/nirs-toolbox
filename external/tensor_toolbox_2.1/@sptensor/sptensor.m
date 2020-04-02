function t = sptensor(varargin)
%SPTENSOR Create a sparse tensor.
%
%   X = SPTENSOR(SUBS, VALS, SZ, FUN) uses the rows of SUBS and VALS
%   to generate a sparse tensor X of size SZ = [m1 m2 ... mn]. SUBS is
%   an p x n array specifying the subscripts of the values to be
%   inserted into S. The k-th row of SUBS specifies the subscripts for
%   the k-th value in VALS. The values are accumulated at repeated
%   subscripts using the function FUN, which is specified by a
%   function handle.
%   
%   There are several simplifications of this four argument call.
%
%   X = SPTENSOR(SUBS,VALS,SZ) uses FUN=@SUM.
%
%   X = SPTENSOR(SUBS,VALS) uses SM = max(SUBS,[],1).
%
%   X = SPTENSOR(SZ) abbreviates X = SPTENSOR([],[],SZ).
%
%   X = SPTENSOR(Y) copies/converts Y if it is an sptensor, an sptenmat, or
%   a dense tensor or MDA (the zeros are squeezed out), or a sparse matrix. 
%   Note that a row-vector, integer MDA is interpreted as a size
%   (see previously constructor).
%
%   S = SPTENSOR is the empty constructor.
%
%   The argument VALS may be scalar, which is expanded to be the
%   same length as SUBS, i.e., it is equivalent to VALS*(p,1).
%
%   Examples
%   subs = [1 1 1; 1 1 3; 2 2 2; 4 4 4; 1 1 1; 1 1 1]
%   vals = [0.5; 1.5; 2.5; 3.5; 4.5; 5.5]
%   siz = [4 4 4];
%   X = sptensor(subs,vals,siz) %<-- sparse 4x4x4, repeats summed
%   X = sptensor(subs,1,siz) %<-- scalar 2nd argument
%   X = sptensor(subs,vals,siz,@max) %<-- max for accumulation
%   myfun = @(x) sum(x) / 3;
%   X = sptensor(subs,vals,siz,myfun) %<-- custom accumulation
%
%   See also SPTENRAND, TENSOR, SPTENMAT, ACCUMARRAY
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
% $Id: sptensor.m,v 1.13 2006/11/26 00:16:28 tgkolda Exp $

% EMPTY Constructor
if nargin == 0
    t.subs = [];
    t.vals = [];
    t.size = [];
    t = class(t,'sptensor');
    return;
end

% SINGLE ARGUMENT
if (nargin == 1)

    source = varargin{1};

    switch(class(source))
        
        case 'sptensor',                % COPY CONSTRUCTOR
            t.subs = source.subs;
            t.vals = source.vals;
            t.size = source.size;
            t = class(t, 'sptensor');
            return;
        
        case 'sptenmat',                % CONVERT SPTENMAT
            
            % Extract the tensor size and order
            siz = source.tsize;
    
            % Convert the 2d-subscipts into nd-subscripts
            asubs = source.subs;
            if ~isempty(source.rdims)
                subs(:,source.rdims) = ...
                    tt_ind2sub(siz(source.rdims),asubs(:,1));
            end
            if ~isempty(source.cdims)
                subs(:,source.cdims) = ...
                    tt_ind2sub(siz(source.cdims),asubs(:,2));
            end
            
            % Copy the values (which do not need to be modified)
            vals = source.vals;
            
            % Store everything
            t.subs = subs;
            t.vals = vals;
            t.size = siz;
            t = class(t, 'sptensor');
            return;
            
        case 'tensor',                  % CONVERT TENSOR
            t.subs = find(source);      % find nonzeros
            t.vals = source(t.subs,'extract');    % extract nonzeros
            t.size = size(source);      
            t = class(t, 'sptensor');
            return;
                      
        case {'numeric','logical','double'},                  

            % Case 1: SPARSE MATRIX
            if issparse(source)     
                [i,j,s] = find(source);
                siz = size(source);
                t.subs = [i, j];
                t.vals = s;
                t.size = siz;
                t = class(t,'sptensor');
                return;                
            end
            
            % Case 2: SPECIFYING THE SIZE
            if ndims(source) == 2 && size(source,1) == 1 && ...
                    isequal(source, floor(source))
                t.subs = [];
                t.vals = [];
                t.size = source;
                t = class(t, 'sptensor');
                return;
            end
            
            % Case 3: An MDA
            t = sptensor(tensor(source));
            return;
           
    end % switch

end % nargin == 1


% CONVERT A SET OF INPUTS
if (nargin == 2) || (nargin == 3) || (nargin == 4)
    
    % Extract the subscripts, values, and sizes
    subs = varargin{1};
    vals = varargin{2};
    
    if (nargin > 2)
        siz = varargin{3};
    else
        siz = max(subs,[],1);
    end
    
    if (nargin == 4)
	fun = varargin{4};
    else
	fun = @sum;
    end

    if isempty(subs)
        newsubs = [];
        newvals = [];
    else
        % Identify only the unique indices
        [newsubs,junk,loc] = unique(subs,'rows');

        % Sum the corresponding values
        newvals = accumarray(loc,vals,[size(newsubs,1) 1],fun);
    end
    
    % Find the nonzero indices of the new values
    nzidx = find(newvals);
    newsubs = newsubs(nzidx,:);
    newvals = newvals(nzidx);
    
    % Store everything
    t.subs = newsubs;
    t.vals = newvals;
    t.size = siz;

    % Create the tensor
    t = class(t, 'sptensor');

    return;
end

error('Unsupported use of function SPTENSOR.');

