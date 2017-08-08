function [b, errstr, errid] = validate(this, nsections, diffdelay)
%VALIDATE 

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 3
    
    % If we are not given the differential delay, assume that the first
    % comb of the first channel is the correct size.
    diffdelay = length(this.Comb(1,1).Value);
    if nargin < 2
        
        % if we are not given the number of sections, assume that the
        % length of the first channel's integrator is correct.
        nsections = size(this.Integrator, 1);
    end
end

% Check that both the Integrator and Comb have the correct # of channels
if size(this.Integrator, 2) ~= size(this.Comb, 2)
    error(message('signal:filtstates:cic:validate:nchannelsMismatch', 'Comb', 'Integrator'))
end

% If this channel's Integrator or Comb length is not equal to
% NSections, error out.
if size(this.Integrator, 1) ~= nsections || ...
        size(this.Comb, 1) ~= nsections      
    error(message('signal:filtstates:cic:validate:lengthMismatch', 'Integrator', 'Comb'))  
end
    
% Check that each section's comb has a number of values equal to the
% differential delay.
for indx = 1:size(this.Comb, 1)
    for jndx = 1:size(this.Comb, 2)
        if length(this.Comb(indx, jndx).Value) ~= diffdelay
            error(message('signal:filtstates:cic:validate:lengthMismatchRows', 'Comb'))  
        end
    end
end

b = true;
errstr = '';
errid = '';

% [EOF]
