function Hd = parallel(varargin)
%PARALLEL Create a parallel system of discrete-time filter objects.
%   Hd = PARALLEL(Hd1, Hd2, etc) constructs a parallel system of the filter
%   objects Hd1, Hd2, etc.  The block diagram looks like:
%
%           |->  Hd1 ->|
%           |          |
%      x ---|->  Hd2 ->|--> y
%           |          |
%           |-> etc. ->|
%
%   The filters Hd1, Hd2, ... must be operating either in double-precision
%   floating-point or single-precision floating-point.
%
%   Hd1, Hd2, ... must be either single-rate filters or multirate filters
%   in which case the rate change of each stage in the parallel structure
%   must be the same. Note that multirate filters require the DSP System
%   Toolbox.
%
%   Hd1, Hd2, ... can also be parallel or cascade filters themselves.
%
%   % EXAMPLE:
%   k1 = [-0.0154    0.9846   -0.3048    0.5601]; 
%   Hd1 = dfilt.latticeallpass(k1);
%   k2 = [-0.1294    0.8341   -0.4165];
%   Hd2 = dfilt.latticeallpass(k2);
%   Hpar = parallel(Hd1 ,Hd2);
%   x = randn(100,1); % Create a random input signal
%   y = filter(Hpar,x);
%   realizemdl(Hpar)    % Requires Simulink
%
%   See also DFILT/STRUCTURES 
  
%   Copyright 1988-2010 The MathWorks, Inc.

if nargin == 0,
    varargin = {dfilt.dffir(1),dfilt.dffir(1)};
end

Hd = dfilt.parallel;

Hd.FilterStructure = 'Parallel';

% Check that all are dfilts before starting to set parameters.
for k=1:length(varargin)
  if isnumeric(varargin{k})
    varargin{k} = dfilt.scalar(varargin{k});
  end
  if ~(isa(varargin{k}(end),'dfilt.abstractfilter') || isa(varargin{k}(end),'dfilt.multistage'))
      error(message('signal:dfilt:parallel:parallel:DFILTErr'));
  end
end

for k=1:length(varargin)
  Hd.Stage = [Hd.Stage; varargin{k}(:)];
end

checkvalidparallel(Hd);

