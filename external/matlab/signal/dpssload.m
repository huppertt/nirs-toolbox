function [E,V] = dpssload(N,NW)
%DPSSLOAD  Load discrete prolate spheroidal sequences from database.
%   [E,V] = DPSSLOAD(N,NW) are the  DPSSs E and their concentrations V, with 
%   length N and time-halfbandwidth product NW, as stored in the DPSS MAT-file 
%   database, 'dpss.mat'.  
%
%   % Example:
%   %   Create and load a Slepian sequence database in current directory.
% 
%   seq_length=[512,1024];
%   time_halfBW=[2.5,2.5];
%   [dps_seq1,lambda1]=dpss(seq_length(1),time_halfBW(1));
%   [dps_seq2,lambda2]=dpss(seq_length(2),time_halfBW(2));
% 
%   % Create databased dpss.mat in current working directory
%   dpsssave(time_halfBW(1),dps_seq1,lambda1);
%   dpsssave(time_halfBW(2),dps_seq2,lambda2);
% 
%   clear dps_seq1 lambda1      % clear workspace
%   [e,v] = dpssload(512,2.5);  % load sequences - length 512, timeBW 2.5
%   dpssclear(512,2.5);     % Remove from database
%   dpssclear(1024,2.5);    % Remove from database
%
%   See also DPSS, DPSSSAVE, DPSSDIR, DPSSCLEAR.

%   Author: T. Krauss
%   Copyright 1988-2004 The MathWorks, Inc.

index = dpssdir(N,NW);
if isempty(index)
    error(message('signal:dpssload:SignalErr'))
else
    key = index.wlist.key;
    data = load('dpss');
    E = data.(sprintf('E%g', key));
    V = data.(sprintf('V%g', key));
end

