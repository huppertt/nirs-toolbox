function varargout = buffer(varargin) %#ok
%BUFFER Buffer a signal vector into a matrix of data frames.
%  Y = BUFFER(X,N) partitions signal vector X into nonoverlapping data
%  segments (frames) of length N.  Each data frame occupies one
%  column in the output matrix, resulting in a matrix with N rows.
%
%  Y = BUFFER(X,N,P) specifies an integer value P which controls the amount
%  of overlap or underlap in the buffered data frames.
%  - If P>0, there will be P samples of data from the end of one frame
%    (column) that will be repeated at the start of the next data frame.
%  - If P<0, the buffering operation will skip P samples of data after each
%    frame, effectively skipping over data in X, and thus reducing the
%    buffer "frame rate".
%  - If empty or omitted, P is assumed to be zero (no overlap or underlap).
%
%  Y = BUFFER(X,N,P,OPT) specifies an optional parameter OPT used in the
%  case of overlap or underlap buffering.
%  - If P>0 (overlap), OPT specifies an initial condition vector with
%    length equal to the overlap, P.  If empty or omitted, the initial
%    condition is assumed to be zeros(P,1).  Alternatively, OPT can be
%    set to 'nodelay' which removes initial conditions from the output and
%    starts buffering at the first data sample in X.
%  - If P<0 (underlap), OPT specifies a scalar offset indicating the number
%    of initial data samples in X to skip.  The offset must be in the range
%    [0, -P], where P<0.  If empty or omitted, the offset is zero.
%
%  Y = BUFFER(X, ...) returns a matrix of data frames, one frame per
%  column.  If X does not contain a multiple of (N-P) samples, zeros
%  will be appended to the last frame as needed.
%
%  [Y,Z] = BUFFER(X, ...) forces Y to contain only complete frames of data,
%  and returns any remaining samples (i.e., a partial frame) in vector Z.
%  The orientation of vector Z will be the same as X.  Z will be an empty
%  vector if the length of X is a multiple of (N-P).
%
%  [Y,Z,OPT] = BUFFER(X, ...) returns an options argument OPT which can
%  be used in the next call to BUFFER, and is primarily intended for use
%  in continuous buffering applications.
%
%  EXAMPLES:
%
%  % Example 1: Buffering with overlap.
%               % Buffer an entire signal vector into 8-sample frames,
%               % each frame overlapping the previous one by 4 samples.
%               x = 1:18;            % Example input data to be buffered
%               y = buffer(x, 8, 4); % Create overlapping buffer matrix
%
%  % Example 2: Buffering with underlap.
%               % Buffer an entire signal vector into 8-sample frames,
%               % skipping 4 samples between frames.  Return any partial
%               % buffer separately.
%               x = 1:40;                 % Example input data to be buffered
%               [y,z] = buffer(x, 8, -4); % Return last partial frame in z
%
%  % Example 3: Continuous buffering.
%               % Buffer a signal which itself is sequentially obtained one frame
%               % at a time.  We are going to obtain consecutive frames of 50 samples
%               % each, and we wish buffer these into frames of 24 samples with 8
%               % samples of overlap.  This "rebuffering" operation frequently arises
%               % when frame-based data is obtained from another source, such as a
%               % data acquisition peripheral.
%
%               % Notice in particular that, depending on the choice of parameters
%               % for size and overlap, successive buffers may contain differing
%               % numbers of frames (columns).  This is to be expected in general.
%               % Careful choice of rebuffering parameters may prevent this.
%
%               % Create a buffer of data to mimic an external data acquisition
%               x = buffer(1:1000, 50); % 50-sample nonoverlapping frames of data
%
%               % Loop over each frame of source data, to mimic the sequential
%               % arrival of each single frame of data:
%               z = [];  opt = [];
%               for i=1:size(x,2), % Loop over each source frame (column)
%               acq = x(:,i);      % Assume that this is what our data
%                                  % acquisition board returns
%
%               % y will contain a matrix of "rebuffered" data frames
%               % NOTE: For the first loop iteration, z and opt are empty
%               [y,z,opt] = buffer([z;acq], 24, 8, opt);
%               disp(y);        % Do something with the buffer of data
%               pause
%               end

% Author: D. Orofino
% Copyright 1988-2008 The MathWorks, Inc.

% The following comment, MATLAB compiler pragma, is necessary to avoid
% compiling this file instead of linking against the MEX-file.  Don't
% remove.
%# mex

error(message('signal:buffer:NotSupported'));

% [EOF] buffer.m

