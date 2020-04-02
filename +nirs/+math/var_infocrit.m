function crit = var_infocrit( LogL, num_obs, num_var, num_param, criterion )
% Calculate information criterion from Log-Likelihood (BIC, AIC, AICc, CAIC, MAX)
% crit = mv_infocrit(LogL, num_obs, num_var, num_param, criterion)
%
% References -- 
% Eqn 7 and 8 of:
% Hurvich, Clifford M., and Chih?Ling Tsai. "A corrected Akaike information 
% criterion for vector autoregressive model selection." Journal of time series 
% analysis 14.3 (1993): 271-279.
% Eqn 2.2 of:
% Yanagihara, Hirokazu, Hirofumi Wakaki, and Yasunori Fujikoshi. "A consistency 
% property of the AIC for multivariate linear models when the dimension and the 
% sample size are large." Electronic Journal of Statistics 9.1 (2015): 869-897.
% 

if nargin<4
    criterion = 'BIC';
end

LogL = LogL + num_obs .* num_var;
unknown_params = (num_param .* num_var.^2 + num_var .* (num_var + 1)/2);
switch upper(criterion)
    case 'BIC'
        crit = LogL + unknown_params.* log(num_obs);
    case 'AIC'
        crit = LogL + 2.*unknown_params;
    case 'AICC'
        crit = LogL + 2.*num_obs.*unknown_params./(num_obs - num_var - num_param .* num_var - 1);
        crit((num_obs - num_var - num_param .* num_var - 1) <= 0) = nan;
    case 'CAIC'
        crit = LogL + unknown_params .* (1+log(num_obs));
    otherwise
        error('Unknown model selection criterion: %s',criterion);
end
end