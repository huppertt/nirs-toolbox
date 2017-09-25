function corrected = tddc( signal )

[L,nch] = size(signal);
if nch>1
    corrected = zeros(size(signal));
    for ch = 1:nch
        corrected(:,ch) = nirs.math.tddc( signal(:,ch) );
    end
    return
end

%% Compute temporal derivative of the signal
tmp = diff(signal);

%% Determine rank of each value
[~,ind] = sort(tmp);
[~,rank] = sort(ind);

%% Get distribution parameters
mu = median(tmp);
sigma = 1.4826 .* mad(tmp,1);
pd = makedist('normal','mu',mu,'sigma',sigma);

%% Generate corrected values based on expectation from distribution
new_deriv = icdf( pd , rank ./ (L+1) );

%% Integrate corrected derivative to create corrected signal
corrected = cumsum([0; new_deriv]);

end