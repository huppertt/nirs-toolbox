function [P,N] = seqperiod(X)  %#ok<STOUT,INUSD>
%SEQPERIOD Find minimum-length repeating sequence in a vector.
% 
%  P = SEQPERIOD(X) returns the index P of the sequence of samples
%  X(1:P) which is found to repeat (possibly multiple times) in
%  X(P+1:end).  P is the sample period of the repetitive sequence.
%  No intervening samples may be present between repetitions.  An
%  incomplete repetition is permitted at the end of X.  If no
%  repetition is found, the entire sequence X is returned as the
%  minimum-length sequence and hence P=length(X).
%
%  [P,N] = SEQPERIOD(X) returns the number of repetitions N of the
%  sequence X(1:P) in X.  N will always be >= 1 and may be non-
%  integer valued.
%
%  If X is a matrix or N-D array, the sequence period is determined
%  along the first non-singleton dimension of X.
%
%   % Example 1:
%   %   Define data and find the minimum-length repeating sequence in it.
%
%   x = repmat([32,43,54],1,4)      % Defining data 
%   p = seqperiod(x)                % Minimum-length repeating sequence
%
%   % Example 2:
%   %   Find the period of each of the column-subsequences of the matrix.
%
%   x = [4 0 1 6; 
%        2 0 2 7; 
%        4 0 1 5; 
%        2 0 5 6];
%   p = seqperiod(x)

%  Author: D. Orofino
%  Copyright 1988-2012 The MathWorks, Inc.

% The following comment, MATLAB compiler pragma, is necessary to avoid
% compiling this file instead of linking against the MEX-file.  Don't
% remove.
%# mex

error(message('signal:seqperiod:NotSupported'));

% [EOF] seqperiod.m
