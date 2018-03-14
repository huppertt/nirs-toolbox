function [stim_vectors,stim_names] = aux2stim_kmeans(aux,threshold,verbose)
% This function generate stimulus events from auxillary data
% [stim_vectors,stim_names] = nirs.util.aux2stim_kmeans( aux, [threshold], [verbose] )
if nargin<3 || isempty(verbose), verbose = false; end
if nargin<2 || isempty(threshold), threshold = .8; end

if isa(aux,'nirs.core.GenericData')
    aux = aux.data;
end

sil = zeros(size(aux,2),1);
num_onsets = nan(size(aux,2),1);
ss = zeros(size(aux));
for i = 1:size(aux,2)
    
    % Perform clustering
    z = zscore(aux(:,i));
    idx = kmeans(z,2);
    
    % Compute the best-case silhouette value
    within_1 = mean(pdist(z(idx==1)));
    within_2 = mean(pdist(z(idx==2)));
    within = min(within_1,within_2);
    between = mean(reshape(pdist2(z(idx==1),z(idx==2)),[],1));
    sil(i) = (between - within) / max(between,within);
    
    % Convert to stim marks, flipping if necessary
    tmp_s = logical( idx - 1 );
    b = [ones(size(z)) tmp_s] \ z;
    if b(2)<0
        tmp_s = ~tmp_s;
    end
    
    num_onsets(i) = sum(diff(tmp_s)==1);
    ss(:,i) = tmp_s;
end

% Determine which to keep
if ischar(threshold) && any(strcmpi({'best','max'},threshold))
    threshold = max(sil);
end
isusable = (sil>=threshold);
good_channels = find(isusable);

% Create stim vectors
stim_vectors = [];
stim_names = {};
for i = 1:length(good_channels)
    stim_vectors = [stim_vectors ss(:,good_channels(i))];
    stim_names = [stim_names {sprintf('aux_channel%i',good_channels(i))}];
end

if verbose
    tbl = table((1:size(aux,2))',sil,num_onsets,isusable,'VariableNames',{'Channel','Silhouette','NumOnsets','Keep'});
    disp(tbl);   
end

end
