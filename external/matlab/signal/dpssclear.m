function dpssclear(N,NW)
%DPSSCLEAR  Remove discrete prolate spheroidal sequences from database.
%   DPSSCLEAR(N,NW) removes the DPSSs with length N and time-halfbandwidth 
%   product NW, from the DPSS MAT-file database, 'dpss.mat'.  
%
%   % Example:
%   %   Create and manage Slepian sequences database in current directory.
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
%   dpssdir             % dpss database directory, before deletion
%   dpssclear(512,2.5)  % remove length-512,halftimeBW = 2.5 from database
%   dpssdir             % dpss database directory, after deletion
%   dpssclear(1024,2.5) % remove length-1024,halftimeBW = 2.5 from database
%
%   See also DPSS, DPSSSAVE, DPSSLOAD, DPSSDIR.

%   Author: T. Krauss
%   Copyright 1988-2004 The MathWorks, Inc.

narginchk(2,2) 
index = dpssdir;

if ~isempty(index)
    w = which('dpss.mat');

    i = find([index.N] == N);
    if isempty(i)
        error(message('signal:dpssclear:SignalErrL'))
    end
    j = find([index(i).wlist.NW] == NW);
    if isempty(j)
        error(message('signal:dpssclear:SignalErrLNW', 'NW', sprintf( '%g', NW )))
    end

    key = index(i).wlist(j).key;
    index(i).wlist(j) = [];
    if isempty(index(i).wlist)
        index(i) = []; %#ok<NASGU>
    end

    str = sprintf('E%g = []; V%g = [];',key,key);
    eval(str)
    save(w, sprintf('E%g', key), sprintf('V%g', key), 'index', '-append');
end
