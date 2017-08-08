function d = dpssdir(N,NW)
%DPSSDIR  Discrete prolate spheroidal sequence database directory.
%   DPSSDIR lists the directory of saved DPSSs in the file dpss.mat.
%   DPSSDIR(N) lists the DPSSs saved with length N.
%   DPSSDIR(NW,'NW') lists the DPSSs saved with time-halfbandwidth product NW.
%   DPSSDIR(N,NW) lists the DPSSs saved with length N and time-halfbandwidth 
%   product NW.
%
%   INDEX = DPSSDIR is a structure array describing the DPSS database.
%   Pass N and NW options as for the no output case to get a filtered INDEX.
%
%   % Example:
%   %   List the directory containing a Slepian sequence database.
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
%   dpssdir             % dpss database directory
%   dpssclear(512,2.5);     % Remove from database
%   dpssclear(1024,2.5);    % Remove from database
%
%   See also DPSS, DPSSSAVE, DPSSLOAD, DPSSCLEAR.

%   Author: T. Krauss
%   Copyright 1988-2009 The MathWorks, Inc.

N_FIXED = 0;
NW_FIXED = 0;

if nargin == 1
    N_FIXED = 1;
end
if nargin == 2
    if ischar(NW)
        NW_FIXED = 1;
        NW = N;
    else
        N_FIXED = 1;
        NW_FIXED = 1;
    end
end

index = [];
w = which('dpss.mat','-all');

doubled = 0;
if iscell(w)
    for i=2:length(w)
        doubled = ~strcmp(w{1},w{i});
        if doubled, break, end
    end
end

if doubled && length(w)>1
    warning(message('signal:dpssdir:Ignore', w{ 1 }))
end

if isempty(w)      % new dpss database
    if nargout == 0
        error(message('signal:dpssdir:FileNotFound'));
    end
else     % add this to existing dpss
    w = w{1};
    load(w, 'index');
    
    if nargout == 0
        disp(sprintf('%s',getString(message('signal:dpssdir:File',w))))
        disp(['    N     NW    ' getString(message('signal:dpssdir:VariableNames'))])
        disp('   ---   ----  ----------------')
    end

    [Nsort,ind] = sort([index.N]);
    index = index(ind);

    if N_FIXED
        ind = find(Nsort==N);
        index = index(ind);
    end

    for i = 1:length(index)
      [wlist,ind] = sort([index(i).wlist.NW]);
      index(i).wlist = index(i).wlist(ind);
      if NW_FIXED
          ind = find(wlist==NW);
          index(i).wlist = index(i).wlist(ind);
      end
 
      if nargout == 0
          for j = 1:length(index(i).wlist)

            key = index(i).wlist(j).key;
        
            args = {index(i).wlist(j).NW, key, key};
            if j == 1
                args = {index(i).N, args{:}};
                str = sprintf('%7.0f %5.2f  E%g, V%g', args{:}); %#ok
            else
                str = sprintf('        %5.2f  E%g, V%g', args{:}); %#ok
            end
            disp(str)

          end
      end
      
    end
    for i = length(index):-1:1
        if isempty(index(i).wlist)
            index(i) = [];
        end
    end
    if isempty(index) && nargout == 0
        error(message('signal:dpssdir:SignalErr'));
    end
end
    

if nargout > 0
    d = index;
end
