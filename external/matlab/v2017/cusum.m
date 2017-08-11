function [iUpper, iLower, upperSum, lowerSum] = cusum(x, varargin)
%CUSUM  Detect small changes in mean via cumulative sums
% [IUPPER, ILOWER] = CUSUM(X) returns the first index of the upper and
% lower cumulative sums of vector X that have moved beyond 5 standard
% deviations above and below the target mean, respectively.  IUPPER and
% ILOWER will be empty if all indices are within the tolerance.  The
% minimum detectable mean shift to detect is set to one standard deviation.
% The target mean and standard deviations are estimated from the first 25
% samples in the input signal.
%
% [IUPPER, ILOWER] = CUSUM(X, CLIMIT) sets the control limits that the
% upper and lower cumulative sums are allowed to drift from the mean.
% Specify the control limit, CLIMIT, in units of standard deviations from
% the mean. The default value of CLIMIT is 5.
%
% [IUPPER, ILOWER] = CUSUM(X, CLIMIT, MSHIFT) sets the minimum mean shift
% to detect.  Specify the mean shift, MSHIFT, in units of standard
% deviations away from the mean.  The default value of MSHIFT is 1.
%
% [IUPPER, ILOWER] = CUSUM(X, CLIMIT, MSHIFT, TMEAN) specifies the target
% mean of the overall signal from which to make the baseline measurement.
% If unspecified, TMEAN is computed as the mean of the first 25 samples of
% X.
%
% [IUPPER, ILOWER] = CUSUM(X, CLIMIT, MSHIFT, TMEAN, TDEV) specifies the
% target standard deviation from which to compute the upper and lower
% control limits.  If unspecified, TDEV is computed as the first 25 samples
% of X.
%
% [IUPPER, ILOWER] = CUSUM(...,'all') will return all indices of the upper
% and lower cumulative sums that are beyond the control limit.
%
% [IUPPER, ILOWER, UPPERSUM, LOWERSUM] = CUSUM(...) additionally returns
% the upper and lower cumulative sums.
%
% CUSUM(...) without output arguments plots the upper and lower
% cumulative sums normalized to one standard deviation above and below the
% target mean, respectively.
%
%   % Example:
%   %    Detect change in mean of a signal that has a mean shift in 
%   %    its data at the 50th sample point the signal.
%   load('meanshift','x')
%   plot(x)
%   title('Original Signal')
%   figure
%   cusum(x)
%
%   See also MEAN, FINDCHANGEPTS.

%   Copyright 2015 The MathWorks, Inc.

narginchk(1,6);

% validate input vector
validateattributes(x,{'numeric'},{'vector','real','finite'},'cusum','X',1);

% get presence 'all' input option
[allLimits,varargin] = getAllOption(varargin{:});
chkunusedopt(varargin);

% fetch and validate numeric arguments
[climit,mshift,tmean,tdev] = getNumericArguments(x,varargin{:});

% compute the upper and lower cumulative sums
[uppersum, lowersum] = cusums(x, tmean, mshift, tdev);

% fetch the indices where the sums exceed the borders
[iupper,ilower] = violations(uppersum,lowersum,climit,tdev,allLimits);

if nargout==0
  % plot control chart
  plotchart(uppersum,lowersum,tmean,tdev,climit,iupper,ilower);
else
  % copy variables to output
  iUpper = iupper;
  iLower = ilower;
  upperSum = uppersum;
  lowerSum = lowersum;
end

% -------------------------------------------------------------------------
function [allOption,varargin] = getAllOption(varargin)
allOption = false;

i = 1;
while i <= numel(varargin)
  if ischar(varargin{i}) && strncmpi(varargin{i},'all',length(varargin{i}))
    allOption = true;
    varargin(i)=[];
  else
    i = i+1;
  end
end

% -------------------------------------------------------------------------
function [climit,mshift,tmean,tdev] = getNumericArguments(x,varargin)

nvarargin = numel(varargin);
if nvarargin<4
  tdev = max(eps,std(x(1:min(25,length(x)))));
else
  tdev = varargin{4};
end

if nvarargin<3
  tmean = mean(x(1:min(25,length(x))));
else
  tmean = varargin{3};
end

if nvarargin<2
  mshift = 1;
else
  mshift = varargin{2};
end

if nvarargin<1
  climit = 5;
else
  climit = varargin{1};
end

validateattributes(climit,{'numeric'},{'scalar','real','positive','finite'},'cusum','CLIMIT',2);
validateattributes(mshift,{'numeric'},{'scalar','real','positive','finite'},'cusum','MSHIFT',3);
validateattributes(tmean,{'numeric'},{'scalar','real','finite'},'cusum','TMEAN',4);
validateattributes(tdev,{'numeric'},{'scalar','real','positive','finite'},'cusum','TDEV',5);

% -------------------------------------------------------------------------
function [uppersum, lowersum] = cusums(x, tmean, mshift, tdev)

uppersum = zeros(size(x));
lowersum = zeros(size(x));

uppersum(1) = 0;
lowersum(1) = 0;

for i=2:length(x)
  uppersum(i) = max(0, uppersum(i-1) + x(i) - tmean - mshift*tdev/2);
  lowersum(i) = min(0, lowersum(i-1) + x(i) - tmean + mshift*tdev/2);
end

% -------------------------------------------------------------------------
function [iupper,ilower] = violations(uppersum,lowersum,climit,tdev,allLimits)
if allLimits
  iupper = find(uppersum >  climit*tdev);
  ilower = find(lowersum < -climit*tdev);
else
  iupper = find(uppersum >  climit*tdev, 1, 'first');
  ilower = find(lowersum < -climit*tdev, 1, 'first');
end

% -------------------------------------------------------------------------
function  plotchart(uppersum,lowersum,tmean,tdev,climit,iupper,ilower)
newplot
n = length(uppersum);
hLines = plot(1:numel(uppersum),uppersum(:)/tdev,1:numel(lowersum),lowersum(:)/tdev);
line([1 n],climit*[1 1],'LineStyle',':','Color',hLines(1).Color);
line([1 n],-climit*[1 1],'LineStyle',':','Color',hLines(2).Color);

if ~isempty(iupper)
  line(iupper,uppersum(iupper)/tdev,'Color','r','Marker','o','LineStyle','none');
end

if ~isempty(ilower)
  line(ilower,lowersum(ilower)/tdev,'Color','r','Marker','o','LineStyle','none');
end

hAxes = hLines(1).Parent;
hAxes.YLim = [min(-climit-1,hAxes.YLim(1)) max(climit+1,hAxes.YLim(2))];
ylabel('StandardErrors');
xlabel('Samples');

