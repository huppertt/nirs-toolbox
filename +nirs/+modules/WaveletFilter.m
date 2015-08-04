classdef WaveletFilter < nirs.modules.AbstractModule
%% WaveletFilter - Filter to remove outliers (motion) and low freq characteristics.
%
% Options:
%     sthresh         - std dev to define outliers
%     wbasis          - wavelet basis function
%     removeScaling   - flag to remove low freq scaling coefficients
        
    properties
        sthresh         = 5;        % std dev to define outliers
        wbasis          = 'sym8';   % wavelet basis function
        removeScaling   = true;     % flag to remove low freq scaling coefficients
    end
    
    methods

        function obj = WaveletFilter( prevJob )
           obj.name = 'Remove Trend & Motion w/ Wavelets';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            for i = 1:length(data)
                for j = 1:size(data(i).data, 2)
                    y = data(i).data(:,j);

                    % max level
                    n = wmaxlev(size(y), obj.wbasis);

                    % decomposition
                    [w, l] = wavedec(y, n, obj.wbasis);

                    % cumulative indices
                    L = cumsum([0; l(1:n+1)]);
                    
                    % remove lowest freq components
                    if obj.removeScaling
                        w(1:L(2)) = 0;
                    end

                    % thresholding
                    for k = 2:n+1
                       idx = L(k)+1:L(k+1);
                       
                       % selection
                       lst = abs(w(idx)) ./ (mad(w(idx),0)/0.6745);
                       lst = lst > obj.sthresh;
                       
                       w(idx(lst)) = 0;
                    end

                    % estimated signal
                    y = waverec(w, l, 'sym8');

                    data(i).data(:,j) = y;
                end
            end
            
        end
    end
end