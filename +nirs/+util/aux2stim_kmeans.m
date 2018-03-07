function [stim_vectors,stim_names] = aux2stim_kmeans(aux,threshold,verbose)
% This function generate stimulus events from auxillary data
% [stim_vectors,stim_names] = nirs.util.aux2stim_kmeans( aux, [threshold], [verbose] )
if nargin<3, verbose = false; end
if nargin<2, threshold = .001; end

if isa(aux,'nirs.core.GenericData')
    aux = aux.data;
end

if verbose
    fprintf('<strong>%7s  |  %16s  |  %8s</strong>\n','Channel','VarianceRatio','HasStim');
end

stim_vectors = [];
stim_names = {};
for i = 1:size(aux,2)
    
    tmp_s = logical( kmeans(aux(:,i),2) - 1 );
    [b,stats] = robustfit( tmp_s , aux(:,i) );
    varratio = var(stats.resid) / var(aux(:,i));
    if verbose
        fprintf('%7s  |  %16s  |  %8s\n',num2str(i),sprintf('%g',varratio),num2str(varratio<threshold));
    end
    if varratio>=threshold
        continue;
    end
    if b(2)<0
        tmp_s = ~tmp_s;
    end
    stim_vectors = [stim_vectors tmp_s];
    stim_names = [stim_names {sprintf('aux_channel%i',i)}];
end

end
