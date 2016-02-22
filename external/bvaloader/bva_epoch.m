function varargout = bva_epoch(varargin)

%bva_epoch - epoch continuous Brain Vision Data 
%
% function [eeg RT] = bva_epoch(data,trigger,stimCode,t0,t1,sampRate [,response])
%
% Epoch continous BVA data set 
%
% Input:
%	data = continous data  (channel x time)
%	trigger = triggers (cf. bva_readmaker.m)
%	stimCode = relevant trigger code
%	t0 = time pre-stimulus
%	t1 = time post-stimulus
%	sampRate = sampling rate of data
%	response = correct reponse - optional,  required for RTs
%
% Output:
%	eeg = eeg-set (chan x time x trial)
% 	rt = reaction times for trials (only computed if response given)
%
% requires:
%
% see also:  ERRP/io, bva_readmarker.m
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
error(nargchk(0,7,nargin))
%% check number of out arguments
error(nargoutchk(0,2,nargout))

%% check && assign input
varargin{8} = [];


if ~isempty(varargin{1}), data = varargin{1}; else error('ERRP:io:bva_epoch:NoData','No data given');end
if ~isempty(varargin{2}), trigger = varargin{2}; else error('ERRP:io:bva_epoch:NoTrigger','No trigger given');end
if ~isempty(varargin{3}), stimCode = varargin{3}; else error('ERRP:io:bva_epoch:NoCode','No Stimulus Code given');end
if ~isempty(varargin{4}), t0 = varargin{4}; else error('ERRP:io:bva_epoch:NoT0','No presimulus interval given');end
if ~isempty(varargin{5}), t1 = varargin{5}; else error('ERRP:io:bva_epoch:NoT1','No postsimulus interval given');end
if ~isempty(varargin{6}), sampRate = varargin{6}; else error('ERRP:io:bva_epoch:Nofs','No sampling rate given');end



%response is optional
if ~isempty(varargin{7}), response = varargin{7}; else response = []; RT = [];end


%% nessary params
channels = size(data,1);
dT = 1/sampRate*1000;

%% match trigger
indStim = find(trigger(1,:) == stimCode);

%% check if trigger was found, or error
if isempty(indStim),
	error('ERRP:io:bva_epoch:TriggerNotFound',' The requested trigger was not Found');
end

% find trials with incorrect response, if wanted
if ~isempty(response)
	
	%everything but the wanted trigger is considered wrong
	indWrong =  find(trigger(1,indStim+1) ~= response);

	% exclude wrong answers
	indStim(indWrong) =[];
	
	% index of correct trials in stream
	indTime = trigger(2,indStim);

	%reaction time
	RT = trigger(2,indStim+1) - trigger(2,indStim);
	
	%account for samling rate
	RT = RT * 1/sampRate*1000 ;
else

	%than we only need the time
	indTime = trigger(2,indStim);
end

%% the number of trials
trials = numel(indStim);

%% compute offset
offPre = t0/dT;
offPost = t1/dT;

%% loop through channels
for i=1:channels
	for j=1:trials
		eeg(i,:,j) = data(i,indTime(j)-abs(offPre):indTime(j)+offPost);
	end
end

varargout{1} = eeg;
varargout{2} = RT;
