function [s, f] = sdf(y, Pmax, Fs, f)
    
%     Estimate the spectral density function (SDF) assuming and AR model.
%     
%     Args:
%         y  : time-series data
%         
%     Keyword Args:
%         Fs : sampling frequency
%         f  : vector of frequencies to estimate sdf for
%         
%     Returns:
%         s : spectral density
%         f : frequencies
%         
%     Example:
%         y = ar_sim([0, 0.1, 0, 0, 0, 0.7], 1000)
%         s, f = sdf(y, 25)
%         
%         sfft = np.abs(np.fft.fft(y))
%         sfft = sfft*sfft / 1000
%         sfft = sfft[0:501]
%         
%         plt.plot(f, np.concatenate((s[:,np.newaxis], sfft[:,np.newaxis]),1))
    
    n = size(y,1);
    
    if nargin < 3
        Fs = 2.0; % will give normalized freq
    end
    
    if nargin < 4
        f = (0:n)' * Fs / n;
        f = f(1:floor(n/2));
    end
            
    % fit an ar model
    [a, r] = nirs.math.ar_fit(y, Pmax);
    
    % separate ar coefs
    a = a(2:end);
    
    % just see the equation here: 
    % https://en.wikipedia.org/wiki/Spectral_density_estimation#Parametric_estimation
    k = (0:length(a)-1)';
    s = bsxfun(@times, k, f');
    s = -2*1i*pi*s / Fs;
    s = exp(s);
    s = bsxfun(@times, a, s);
    s = abs(1 - sum(s,1))' .^ 2;
    s = (var(r)/Fs) ./ s;
    
    % correct the dc term
    s(1) = sum(y) ^ 2 / n / Fs;