function [data, truth] = simMotionArtifact( data , spikes_per_minute , shifts_per_minute )

if nargin<1 || ~exist('data','var') || isempty(data)
   [data, truth]=nirs.testing.simData; 
end

if nargin<2 || ~exist('spikes_per_minute','var') || isempty(spikes_per_minute), spikes_per_minute = 2; end
if nargin<3 || ~exist('shifts_per_minute','var') || isempty(shifts_per_minute), shifts_per_minute = .5; end

num_spikes = round( spikes_per_minute * (data.time(end)-data.time(1))/60 );
num_shifts = round( shifts_per_minute * (data.time(end)-data.time(1))/60 );

nsamp = length(data.time);

spike_inds = randi([2 nsamp-1],[1 num_spikes]);
shift_inds = randi([2 nsamp-1],[1 num_shifts]);

spike_amp_Z = 6 * randn([1 num_spikes]);
shift_amp_Z = 25 * randn([1 num_shifts]);

means = mean(data.data);
data.data = bsxfun(@plus, -log(data.data), log(means));
stds = std(data.data);

for i = 1:num_spikes
    
    amp = abs(spike_amp_Z(i)) .* stds;
    width = 1.9*abs(rand) + .1;
    
    t_peak = data.time(spike_inds(i));
    t_start = t_peak - width;
    t_end = t_peak + width;
    t_inds = find( (data.time > t_start) & (data.time < t_end) );
    
    for j = 1:ceil(length(t_inds)/2)
        
        t = data.time(t_inds(j)) - t_start;
        data.data(t_inds(j),:) = data.data(t_inds(j),:) + amp .^ (t/width);
        data.data(t_inds(end-j+1),:) = data.data(t_inds(end-j+1),:) + amp .^ (t/width);
    end
    
end

for i = 1:num_shifts
    
    shift_amt = shift_amp_Z(i) .* stds;
    data.data(shift_inds(i):end,:)  = bsxfun( @plus , data.data(shift_inds(i):end,:) , shift_amt );
    
end

data.data = exp( -bsxfun(@minus, data.data, log(means)) );

if(nargout>1)
    if(~exist('truth','var'))
        truth=[];
    end

end