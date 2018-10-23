function [data, truth] = simMotionArtifact( data , spikes_per_minute , shifts_per_minute )

if nargin<1 || ~exist('data','var') || isempty(data)
   [data, truth]=nirs.testing.simData; 
end

if nargin<2 || ~exist('spikes_per_minute','var') || isempty(spikes_per_minute), spikes_per_minute = 2; end
if nargin<3 || ~exist('shifts_per_minute','var') || isempty(shifts_per_minute), shifts_per_minute = .5; end

if length(data)>1
    for i = 1:length(data)
        data(i) = nirs.testing.simMotionArtifact( data(i) , spikes_per_minute , shifts_per_minute );
    end
    return
end

num_spikes = round( spikes_per_minute * (data.time(end)-data.time(1))/60 );
num_shifts = round( shifts_per_minute * (data.time(end)-data.time(1))/60 );

nsamp = length(data.time);

spike_inds = randi([2 nsamp-1],[1 num_spikes]);
shift_inds = randi([2 nsamp-1],[1 num_shifts]);

spike_amp_Z = 25 * randn([1 num_spikes]);
shift_amp_Z = 25 * randn([1 num_shifts]);

mu = mean(data.data);
stds = std(data.data);

for i = 1:num_spikes
    
    width = 9.9*rand + .1; % Spike duration of 0.1-10 seconds
    
    t_peak = data.time(spike_inds(i));
    t_start = t_peak - width/2;
    t_inds = find( (data.time > t_start) & (data.time <= t_peak) );
    spike_time = [data.time(t_inds); data.time(t_inds(end-1:-1:1))] - t_start;
    
    amp = abs(spike_amp_Z(i));
    amp = amp + .25*amp*randn([1 length(stds)]); % Jitter relative spike magnitude across channels
    amp = amp .* stds;
    amp = amp + bsxfun(@times, 2*rand(length(spike_time),length(stds))-1, .5*amp); % Add temporal jitter
    
    tau = width/2;

    spike_data = amp .^ (spike_time./tau);

    t_inds = t_inds(1):min(t_inds(1)+size(spike_data,1)-1,size(data.data,1));
    
    data.data(t_inds,:) = data.data(t_inds,:) + spike_data(1:length(t_inds),:);
        
end

for i = 1:num_shifts
    
    shift_amt = shift_amp_Z(i) .* stds;
    data.data(shift_inds(i):end,:)  = bsxfun( @plus , data.data(shift_inds(i):end,:) , shift_amt );
    
end

% Restore original mean intensity
data.data = bsxfun(@plus,bsxfun(@minus,data.data,mean(data.data)),mu);

% Prevent negative intensities
while any(data.data(:)<=0)
    has_neg = any(data.data<0);
    data.data(:,has_neg) = data.data(:,has_neg) + rand(1,sum(has_neg)) .* std(data.data(:,has_neg));
end

if(nargout>1)
    if(~exist('truth','var'))
        truth=[];
    end

end