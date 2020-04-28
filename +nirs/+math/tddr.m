<<<<<<< HEAD
function signal_corrected = tddr( signal , sample_rate )
=======
function signal_corrected = tddr( signal , sample_rate,splitPosNeg )
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
% Perform Temporal Derivative Distribution Repair (TDDR) motion correction
%
% Usage:
%   signals_corrected = nirs.math.tddr( signals , sample_rate );
%
% Inputs:
%   signals: A [sample x channel] matrix of uncorrected optical density data
%   sample_rate: A scalar reflecting the rate of acquisition in Hz
%
% Outputs:
%   signals_corrected: A [sample x channel] matrix of corrected optical density data
%
%   Fishburn F.A., Ludlum R.S., Vaidya C.J., & Medvedev A.V. (2019). 
%   Temporal Derivative Distribution Repair (TDDR): A motion correction 
%   method for fNIRS. NeuroImage, 184, 171-179.
%   https://doi.org/10.1016/j.neuroimage.2018.09.025

<<<<<<< HEAD
=======

if(nargin<3)
    splitPosNeg=false;
end

>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
%% Iterate over each channel
nch = size(signal,2);
if nch>1
    signal_corrected = zeros(size(signal));
    for ch = 1:nch
<<<<<<< HEAD
        signal_corrected(:,ch) = nirs.math.tddr( signal(:,ch) , sample_rate );
=======
        signal_corrected(:,ch) = nirs.math.tddr( signal(:,ch) , sample_rate,splitPosNeg );
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
    end
    return
end

DC=median(signal);
signal=signal-DC;

%% Preprocess: Separate high and low frequencies
filter_cutoff = .5;
filter_order = 3;
Fc = filter_cutoff * 2/sample_rate;
if Fc<1
    [fb,fa] = butter(filter_order,Fc);
    signal_low = filtfilt(fb,fa,signal);
else
    signal_low = signal;
end
signal_high = signal - signal_low;

%% Initialize
tune = 4.685;
D = sqrt(eps(class(signal)));
mu = inf;
iter = 0;

%% Step 1. Compute temporal derivative of the signal
deriv = diff(signal_low);

%% Step 2. Initialize observation weights
w = ones(size(deriv));

%% Step 3. Iterative estimation of robust weights
while iter < 50
    
    iter = iter + 1;
    mu0 = mu;
    
    % Step 3a. Estimate weighted mean
    mu = sum( w .* deriv ) / sum( w );
    
<<<<<<< HEAD
    % Step 3b. Calculate absolute residuals of estimate
    dev = abs(deriv - mu);
=======
    if(splitPosNeg)
    lst=find((deriv - mu)>0);
    % Step 3b. Calculate absolute residuals of estimate
    dev = abs(deriv(lst) - mu);
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f

    % Step 3c. Robust estimate of standard deviation of the residuals
    sigma = 1.4826 * median(dev);

    % Step 3d. Scale deviations by standard deviation and tuning parameter
    r = dev / (sigma * tune);
    
    % Step 3e. Calculate new weights accoring to Tukey's biweight function
<<<<<<< HEAD
    w = ((1 - r.^2) .* (r < 1)) .^ 2;
=======
    w(lst) = ((1 - r.^2) .* (r < 1)) .^ 2;
    
    lst=find((deriv - mu)<=0);
    % Step 3b. Calculate absolute residuals of estimate
    dev = abs(deriv(lst) - mu);

    % Step 3c. Robust estimate of standard deviation of the residuals
    sigma = 1.4826 * median(dev);

    % Step 3d. Scale deviations by standard deviation and tuning parameter
    r = dev / (sigma * tune);
    
    % Step 3e. Calculate new weights accoring to Tukey's biweight function
    w(lst) = ((1 - r.^2) .* (r < 1)) .^ 2;
    
    else
        % Step 3b. Calculate absolute residuals of estimate
        dev = abs(deriv - mu);
        
        % Step 3c. Robust estimate of standard deviation of the residuals
        sigma = 1.4826 * median(dev);
        
        % Step 3d. Scale deviations by standard deviation and tuning parameter
        r = dev / (sigma * tune);
        
        % Step 3e. Calculate new weights accoring to Tukey's biweight function
    w = ((1 - r.^2) .* (r < 1)) .^ 2;
    end
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f

    % Step 3f. Terminate if new estimate is within machine-precision of old estimate
    if abs(mu-mu0) < D*max(abs(mu),abs(mu0))
        break;
    end

end

%% Step 4. Apply robust weights to centered derivative
new_deriv = w .* (deriv-mu);

%% Step 5. Integrate corrected derivative
signal_low_corrected = cumsum([0; new_deriv]);

%% Postprocess: Center the corrected signal
signal_low_corrected = signal_low_corrected - mean(signal_low_corrected);

%% Postprocess: Merge back with uncorrected high frequency component
signal_corrected = signal_low_corrected + signal_high;

signal_corrected=signal_corrected+DC;

end
