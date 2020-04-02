function RT = bva_rt(varargin)

%bva_rt - extract reaction times from EEG trigger 
%
% function RT = bva_rt(trigger,stimCode,response);
%
% Extract the reaction times from a trigger sequence made by 
% bva_readmarker.m. For extraction and epochising see
% bva_epoch.m
%
% Input:
%	trigger = triggers (cf. bva_readmarker.m.)
%	stimCode = stimulus code of interest
%	response = stimulus code of required response
% 	sampRate = sampling rate (otherwise time would be discrete)
%
% Output:
% 	rt = reaction times for correct trials
%
% requires:
%
% see also:  ERRP/io, bva_readtrigger.m
%

% Copyright (C) 2008-2012 Stefan Schinkel, HU Berlin
% http://people.physik.hu-berlin.de/~schinkel/
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% $Log$

% debug settings
debug = 0;
if debug;warning('on','all');else warning('off','all');end

%% check number of input arguments
error(nargchk(3,4,nargin))

%% check number of out arguments
error(nargoutchk(0,1,nargout))

varargin{5} = [];

%% check && assign input
if ~isempty(varargin{1}), trigger = varargin{1}; else error('ERRP:io:bva_rt:NoTrigger','No trigger given');end
if ~isempty(varargin{2}), stimCode = varargin{2}; else error('ERRP:io:bva_rt:NoCode','No Stimulus Code given');end
if ~isempty(varargin{3}), response = varargin{3}; else error('ERRP:io:bva_rt:NoReponse','No Response Code given');;end
if ~isempty(varargin{4}), sampRate = varargin{4}; else sampRate = []; disp('RT will be discrete!!!');end

%% match trigger
indStim = find(trigger(1,:) == stimCode);

%% check if trigger was found, or error
if isempty(indStim),
	error('ERRP:io:bva_epoch:TriggerNotFound',' The requested trigger was not Found');
end

%everything but the wanted trigger is considered wrong
indWrong =  find(trigger(1,indStim+1) ~= response);

% exclude wrong answers
indStim(indWrong) =[];

% index of correct trials in stream
indTime = trigger(2,indStim);

%reaction time
RT = trigger(2,indStim+1) - trigger(2,indStim);

%% account for sampling frequency, if given
if ~isempty(sampRate),
	RT = RT * 1/sampRate*1000;
end
