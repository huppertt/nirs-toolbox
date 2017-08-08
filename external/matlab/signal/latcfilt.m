function [F,G,Zf] = latcfilt(varargin) 
%LATCFILT Lattice and lattice-ladder filter implementation.
%   [F,G] = LATCFILT(K,X) filters X with the FIR lattice coefficients
%   in vector K.  F is the forward lattice filter result, and G is
%   the backward lattice filter result. If K is such that ABS(K) <= 1,
%   F corresponds to the minimum-phase output of the filter and G
%   corresponds to the maximum-phase output.
%
%   If K and X are vectors, the result is a (signal) vector.
%   Matrix arguments are permitted under the following rules:
%   - If X is a matrix and K is a vector, each column of X is processed
%     through the lattice filter specified by K.
%   - If X is a vector and K is a matrix, each column of K is used to
%     filter X, and a signal matrix is returned.
%   - If X and K are both matrices with the same # of columns, then the
%     i-th column of K is used to filter the i-th column of X.  A
%     signal matrix is returned.
%
%   [F,G,Zf] = LATCFILT(K,X,'ic',Zi) gives access to the initial and final
%   conditions, Zi and Zf, of the lattice states.  Zi must be a vector
%   with length(K) entries.
%
%   [F,G] = LATCFILT(K,1,X) filters X with the all-pole or allpass
%   lattice coefficients in vector K.  F is the all-pole lattice filter
%   result and G is the allpass lattice filter result.  K must be a vector,
%   while X may be a signal matrix.
%
%   [F,G,Zf] = LATCFILT(K,1,X,'ic',Zi) gives access to the initial
%   and final conditions, Zi and Zf, of the lattice states.  Zi must be a
%   vector with length(K) entries.
%
%   [F,G] = LATCFILT(K,V,X) filters X with the ARMA lattice-ladder
%   structure described by the lattice coefficients K and ladder
%   coefficients V.  F is the result of adding all outputs from the
%   ladder coefficients and G is the output of the allpass section of
%   the structure.  K and V must be vectors, while X may be a signal
%   matrix.
%
%   [F,G,Zf] = LATCFILT(K,V,X,'ic',Zi) gives access to the initial and
%   final conditions, Zi and Zf, of the lattice states.  Zi must be a
%   vector with length(K) entries.
%   
%   [F,G,Zf] = LATCFILT(K,V,X,'ic',Zi,DIM) processes the data X along
%   dimension DIM.  If you specify DIM, the filter K must be a vector.  You
%   must enter all six input arguments in order.  To omit an argument,
%   specify it as an empty vector [].  Zf returns the final conditions in
%   columns regardless of the shape of input data matrix X.
%   
%   % Example:
%   %   Filter data with an FIR lattice filter. 
%
%   x=randn(512,1);                     % Create data
%   [f,g]=latcfilt([1/2 1],x);          % Reflection coefficients for 
%                                       % 3-point Moving Average filter
%   subplot(2,1,1);                  
%   plot(f);                            % Max-phase lattice filter output
%   title('Max-phase lattice filter output')
%   subplot(2,1,2);        
%   plot(g);                            % Min-phase lattice filter output
%   title('Min-phase lattice filter output')     
%
%   See also FILTER, TF2LATC, LATC2TF.

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(2,6);

if nargin < 6
    [F,G,Zf] = latcfiltmex(varargin{:});
else
    [K,V,X,ic,Zi,DIM] = deal(varargin{:});
    
    if isempty(K)||isempty(X)
        error(message('signal:latcfilt:invalidinput'));
    end
    
    if ~isempty(DIM)
        if ~isvector(K)
            error(message('signal:latcfilt:inputnotsupported'));
        end
        if (DIM < 1) || (DIM > ndims(X))
            error(message('signal:latcfilt:invaliddim', ndims( X )));
        end
        
        s = size(X);

        [X,perm,nshifts] = shiftdata(X,DIM);
        s_shift = size(X); % New size
        X = reshape(X,size(X,1),[]); % Force into 2-D
    end

    if isempty(V)
        if isempty(ic) && isempty(Zi)
            [F,G,Zf] = latcfiltmex(K,X);
        else
            [F,G,Zf] = latcfiltmex(K,X,ic,Zi);
        end
    else
        if isempty(ic) && isempty(Zi)
            [F,G,Zf] = latcfiltmex(K,V,X);
        else
            [F,G,Zf] = latcfiltmex(K,V,X,ic,Zi);
        end
    end

    if ~isempty(DIM)
        F = reshape(F,s_shift);
        F = unshiftdata(F,perm,nshifts);
        F = reshape(F,s);

        G = reshape(G,s_shift);
        G = unshiftdata(G,perm,nshifts);
        G = reshape(G,s);
    end

end



% [EOF] latcfilt.m
