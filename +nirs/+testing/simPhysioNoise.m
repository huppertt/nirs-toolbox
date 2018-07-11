function [data, truth] = simPhysioNoise( data , cardiac_amp , resp_amp , mayer_amp )

if nargin<1 || ~exist('data','var') || isempty(data)
   [data, truth]=nirs.testing.simData; 
end
if nargin<2 || isempty(cariac_amp)
    cardiac_amp = .4;
end
if nargin<3 || isempty(resp_amp)
    resp_amp = .3;
end
if nargin<3 || isempty(mayer_amp)
    mayer_amp = .2;
end

nsamp = length(data.time);
time = data.time;
[nsamp,nchan] = size(data.data);

% Generate cardiac oscillations
cardiac_freq = 1 + .1*randn;
cardiac_phase = cumsum( .01*2*pi*randn(nsamp,1) , 1 );
cardiac_amp = cardiac_amp * std(data.data,0,1);
cardiac_data = cardiac_amp .* sin( 2*pi*cardiac_freq*time + cardiac_phase );

% Generate respiratory oscillations
resp_freq = .25 + .025*randn;
resp_phase = cumsum( .01*2*pi*randn(nsamp,1) , 1 );
resp_amp = resp_amp * std(data.data,0,1);
resp_data = resp_amp .* sin( 2*pi*resp_freq*time + resp_phase );

% Generate Mayer waves
mayer_freq = .1 + .01*randn;
mayer_phase = cumsum( .01*2*pi*randn(nsamp,1) , 1 );
mayer_amp = mayer_amp * std(data.data,0,1);
mayer_data = mayer_amp .* sin( 2*pi*mayer_freq*time + mayer_phase );

% Add physiological noises
data.data = data.data + cardiac_data + resp_data + mayer_data;

if(nargout>1)
    if(~exist('truth','var'))
        truth=[];
    end

end