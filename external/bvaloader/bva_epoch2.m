function varargout = bva_epoch2(varargin)

%bva_epoch2 - response-locked epoching of  continuous Brain Vision Data 
%
% function [eeg RT] = bva_epoch2(data,trigger,response,stimCode,t0,t1,sampRate,)
%
% Epoch continous BVA data set 
%
% Input:
%	data = continous data  (channel x time)
%	trigger = triggers (cf. bva_readmaker.m)
%	respCode = relevant response code
%	stimCode = relevant trigger code
%	t0 = time pre-stimulus
%	t1 = time post-stimulus
%	sampRate = sampling rate of data
%
% Output:
%	eeg = eeg-set (chan x time x trial)
%	RTs = reaction times
% requires:
%
% see also:  ERRP/io, bva_readmarker.m, bva_epoch.m 
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
error(nargoutchk(0,3,nargout))

%% check && assign input
varargin{8} = [];

if ~isempty(varargin{1}), data = varargin{1}; else error('ERRP:io:bva_epoch2:NoData','No data given');end
if ~isempty(varargin{2}), trigger = varargin{2}; else error('ERRP:io:bva_epoch2:NoTrigger','No trigger given');end
if ~isempty(varargin{3}), resCode = varargin{3}; else error('ERRP:io:bva_epoch2:NoResponse','No response code given');end
if ~isempty(varargin{4}), stimCode = varargin{4}; else error('ERRP:io:bva_epoch2:NoCode','No Stimulus Code given');end
if ~isempty(varargin{5}), t0 = varargin{5}; else error('ERRP:io:bva_epoch2:NoT0','No presimulus interval given');end
if ~isempty(varargin{6}), t1 = varargin{6}; else error('ERRP:io:bva_epoch2:NoT1','No postsimulus interval given');end
if ~isempty(varargin{7}), sampRate = varargin{7}; else error('ERRP:io:bva_epoch2:Nofs','No sampling rate given');end


%% nessary params
channels = size(data,1);
dT = 1/sampRate*1000;

% match response trigger
idResponse = find(trigger(1,:) == resCode);

%% check if trigger was found, or error
if isempty(idResponse),
	error('ERRP:io:bva_epoch2:TriggerNotFound',' The requested trigger was not Found');
end


% select subset for the chosen condition
idCondition = find(trigger(1,idResponse-1) == stimCode);

%calculate reaction times
reactionTime = trigger(2,idResponse(idCondition)) - trigger(2,idResponse(idCondition)-1);


%% the number of trials
trials = numel(idCondition);

%select the times
idTime = trigger(2,idResponse(idCondition));
%% compute offset
offPre = t0/dT;
offPost = t1/dT;

%% loop through channels
for i=1:channels
	for j=1:trials
		eeg(i,:,j) = data(i,idTime(j)-abs(offPre):idTime(j)+offPost);
	end
end

varargout{1} = eeg;
varargout{2} = reactionTime;
varargout{3} = idResponse(idCondition)-1;
