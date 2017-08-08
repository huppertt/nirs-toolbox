function [pks,frqs] = findpeaks(this,varargin)
%FINDPEAKS Find local peaks in data
%   PKS = FINDPEAKS(H) finds local peaks in data contained in the DSPDATA
%   object H.
%
%   [PKS,FRQS]= FINDPEAKS(H) also returns the frequencies FRQS at which the
%   PKS occur.
%
%   [...] = FINDPEAKS(H,'MINPEAKHEIGHT',MPH) finds only those peaks that
%   are greater than MINPEAKHEIGHT MPH. Specifying MPH may help in reducing
%   the processing time. MPH is a real valued scalar. The default value of
%   MPH is -Inf.  
%
%   [...] = FINDPEAKS(H,'MINPEAKDISTANCE',MPD) finds peaks that are at
%   least separated by MINPEAKDISTANCE MPD. MPD is a real valued positive
%   scalar specified in frequency units. This parameter may be specified to
%   ignore smaller peaks that may occur in close proximity to a large local
%   peak. For example, if a large local peak occurs at frequency Fp, then
%   all smaller peaks in the range (Fp-MPD, Fp+MPD) are ignored. If not
%   specified, MPD is assigned a value equal to the minimum distance
%   between two consecutive frequency points in the spectrum estimate.
%
%   [...] = FINDPEAKS(H,'THRESHOLD',TH) finds peaks that are at least
%   greater than their neighbhors by the THRESHOLD TH. TH is real valued
%   scalar greater than or equal to zero. The default value of TH is zero.
%
%   [...] = FINDPEAKS(X,'NPEAKS',NP) specifies the maximum number of peaks
%   to be found. NP is an integer greater than zero. If not specified, all
%   peaks are returned.
%
%   [...] = FINDPEAKS(H,'SORTSTR',STR) specifies the direction of sorting
%   of peaks. STR can take values of 'ascend','descend' or 'none'. If not
%   specified, STR takes the value of 'none' and the peaks are returned in
%   the order of their occurrence.
%
%   EXAMPLE
%      f1  = 0.5; f2 = 0.52; f3 = 0.9; f4 = 0.1; n = 0:255; 
%      x   = 10*cos(pi*f1*n)'+ 6*cos(pi*f2*n)'+ 0.5*cos(pi*f3*n)'+ ...
%            + 0.5*cos(pi*f4*n)';
%      H   = spectrum.periodogram;
%      h   = msspectrum(H,x); 
%      pks = findpeaks(h);
%
%      % To ignore peaks below 0.1 
%      [pks, frqs] = findpeaks(h,'MinPeakHeight',0.1);                                                       
%
%   See also DSPDATA/SFDR, FINDPEAKS

%   Copyright 2007-2010 The MathWorks, Inc.

error(nargchk(1,11,nargin,'struct'));
hopts = uddpvparse('dspopts.findpeaks',{'findpeaksopts',this},varargin{:});

TH  = hopts.Threshold;
PD  = hopts.MinPeakDistance;
PH  = hopts.MinPeakHeight;
NP  = hopts.NPeaks;
STR = hopts.SortStr;

Data = this.Data;
F = this.Frequencies;
K = length(F);

% Find equivalent of MinPeakDistance PD in terms of number of data points.
% This conversion is required since the findpeaks function, which is called
% upon later, requires MinPeakDistance in terms of number of data points
% (integer).
S = round(K*PD/(max(F)-min(F)));

% Setting S = 1 (minimum value allowed) if S < 1
if(isempty(S)) || (S < 1)
    S = 1;
end

if(K<=S)
    error(message('signal:dspdata:abstractps:findpeaks:largePeakDistance', 'MinPeakDistance', num2str( (max( F ) - min( F )) )));
else
    [pks,locs] = findpeaks(Data,'MinPeakHeight',PH,'MinPeakDistance',S,...
                'Threshold',TH,'NPeaks',NP,'SortStr',STR);
    frqs = F(locs); 
end

% [EOF]
