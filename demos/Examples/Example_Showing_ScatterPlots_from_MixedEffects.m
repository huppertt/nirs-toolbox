% This example shows how to use the diagnostic info in the group-level
% Mixed effects models
% Note- this example is incomplete and will not run as is.  It is ment to
% copy/paste into your own analysis script

% Assume that we have data loaded and have run the subject level Stats
% models




j = nirs.modules.AR_IRLS();
SubjStats = j.run( hb );



%% group level
j = nirs.modules.MixedEffects( );
j.formula = 'beta ~ -1 + cond + age:cond + (1|subject) + (1|gender)';
j.dummyCoding = 'full';
j.include_diagnostics=true;  % This flag will cause the code to compute the 
% additional LME model info that we can use for diagntstics  
GroupStats = j.run( SubjStats );

% We can see that the variables table has an additional column holding the
% linear model object

%     source    detector    type            cond                   model      
%     ______    ________    ____    _____________________    _________________
% 
%     1         1           808     'Digit_Symbol'           [1x1 LinearModel]
%     1         1           808     'Stroop1'                [1x1 LinearModel]
%     1         1           808     'Stroop2'                [1x1 LinearModel]
%     1         1           808     'Stroop3'                [1x1 LinearModel]
%     1         1           808     'Shift_Attention'        [1x1 LinearModel]
%     1         1           808     'Stroop3-Stroop1'        [1x1 LinearModel]
%     1         1           808     'Stroop2-Stroop1'        [1x1 LinearModel]
%     1         1           808     'Digit_Symbol:age'       [1x1 LinearModel]
%     1         1           808     'Stroop1:age'            [1x1 LinearModel]
%     1         1           808     'Stroop2:age'            [1x1 LinearModel]

% This is a class LinearModel and has the following methods/properties

%  LinearModel Fitted multiple linear regression model.
%     LM = FITLM(...) fits a linear model to data. The fitted model LM is a
%     LinearModel that can predict a response as a linear function of
%     predictor variables and terms created from predictor variables.
%  
%     LinearModel methods:
%         addTerms - Add terms to linear model
%         removeTerms - Remove terms from linear model
%         step - Selectively add or remove terms from linear model
%         anova - Analysis of variance
%         coefCI - Coefficient confidence intervals
%         coefTest - Linear hypothesis test on coefficients
%         predict - Compute predicted values given predictor values
%         feval - Evaluate model as a function
%         random - Generate random response values given predictor values
%         dwtest - Durbin-Watson test for autocorrelation in residuals
%         plot - Summary plot of regression model
%         plotAdded - Plot of marginal effect of a single term
%         plotAdjustedResponse - Plot of response and one predictor
%         plotDiagnostics - Plot of regression diagnostics
%         plotEffects - Plot of main effects of predictors
%         plotInteraction - Plot of interaction effects of two predictors
%         plotResiduals - Plot of residuals
%         plotSlice - Plot of slices through fitted regression surface
%  
%     LinearModel properties:
%         Coefficients - Coefficients and related statistics
%         Rsquared - R-squared and adjusted R-squared
%         ModelCriterion - AIC and other model criteria
%         Fitted - Vector of fitted (predicted) values
%         Residuals - Table containing residuals of various types
%         ResponseName - Name of response
%         PredictorNames - Names of predictors
%         NumPredictors - Number of predictors
%         Variables - Table of variables used in fit
%         NumVariables - Number of variables used in fit
%         VariableNames - Names of variables used in fit
%         VariableInfo - Information about variables used in the fit
%         NumObservations - Number of observations in the fit
%         ObservationNames - Names of observations in the fit
%         ObservationInfo - Information about observations used in the fit
%         Diagnostics - Regression diagnostics
%         MSE - Mean squared error (estimate of residual variance)
%         RMSE - Root mean squared error (estimate of residual standard deviation)
%         DFE - Degrees of freedom for residuals
%         SSE - Error sum of squares
%         SST - Total sum of squares
%         SSR - Regression sum of squares
%         Steps - Stepwise fit results
%         Robust - Robust fit results
%         Formula - Representation of the model used in this fit
%         LogLikelihood - Log of likelihood function at coefficient estimates
%         CoefficientCovariance - Covariance matrix for coefficient estimates
%         CoefficientNames - Coefficient names
%         NumCoefficients - Number of coefficients
%         NumEstimatedCoefficients - Number of estimated coefficients



% We can make scatter plots for each entry of the model using
G.variables.model{1}.plotAdjustedResponse('age');  % This will plot model entry 1 
% Src-1,Det-1,type-808, cond-'Digit_Symbol' verses age
 
disp(G.variables.model{1})
% Linear regression model:
%     beta ~ cond + cond:age
% 
% Estimated Coefficients:
%                               Estimate          SE         tStat     pValue 
%                              ___________    __________    _______    _______
% 
%     cond_Digit_Symbol          0.0027823     0.0015762     1.7652    0.10099
%     cond_Digit_Symbol:age    -0.00056654    0.00033735    -1.6794    0.11694
% 
% Number of observations: 15, Error degrees of freedom: 13
% Root Mean Squared Error: 3.53
% R-squared: 0.178,  Adjusted R-Squared 0.115
% F-statistic vs. constant model: 2.82, p-value = 0.117
% G.variables.model{1}.plotAdjustedResponse('age')

% We can also plot outlier subjects 
G.variables.model{1}.plotDiagnostics

% We can do ANOVA
anova(G.variables.model{1})
%                 SumSq     DF    MeanSq      F       pValue 
%                 ______    __    ______    ______    _______
% 
%     cond        29.765     1    29.765    2.3899    0.14611
%     cond:age    35.125     1    35.125    2.8203    0.11694
%     Error       161.91    13    12.454   


% The ROI average function  will also keep this info

Region{1} = table(NaN,NaN,'VariableNames',{'source','detector'});
names{1}='WholeProbe';

ROItable = nirs.util.roiAverage(GroupStats,Region,names)
%           ROI                 Contrast              Beta            SE        DF         T            p               model              q     
%     ________________    _____________________    ___________    __________    ___    _________    __________    _________________    __________
% 
%     'WholeProbe:808'    'Digit_Symbol'             0.0028592     0.0008011    148       3.5692    0.00048322    [1x1 LinearModel]     0.0011275
%     'WholeProbe:808'    'Digit_Symbol:age'       -9.5464e-06    0.00018148    148    -0.052603       0.95812    [1x1 LinearModel]             1
%     'WholeProbe:808'    'Shift_Attention'          0.0039025    0.00086087    148       4.5332    1.1904e-05    [1x1 LinearModel]    5.5553e-05
%     'WholeProbe:808'    'Shift_Attention:age'    -0.00074776    0.00025726    148      -2.9066     0.0042154    [1x1 LinearModel]      0.007377
%     'WholeProbe:808'    'Stroop1'                  0.0045767    0.00076639    168       5.9717    1.3596e-08    [1x1 LinearModel]     9.517e-08
%     'WholeProbe:808'    'Stroop1:age'             -0.0004216    0.00015593    168      -2.7037     0.0075626    [1x1 LinearModel]      0.011764
%     'WholeProbe:808'    'Stroop2'                  0.0043506    0.00067907    168       6.4066     1.439e-09    [1x1 LinearModel]    2.0147e-08
%     'WholeProbe:808'    'Stroop2-Stroop1'        -0.00019453    0.00067794    168     -0.28694       0.77451    [1x1 LinearModel]       0.98574
%     'WholeProbe:808'    'Stroop2-Stroop1:age'    -1.0912e-05    0.00015553    168    -0.070161       0.94415    [1x1 LinearModel]             1
%     'WholeProbe:808'    'Stroop2:age'            -0.00062437    0.00015892    168      -3.9289    0.00012445    [1x1 LinearModel]    0.00034846
%     'WholeProbe:808'    'Stroop3'                  0.0038515    0.00087451    168       4.4042    1.8849e-05    [1x1 LinearModel]    6.5972e-05
%     'WholeProbe:808'    'Stroop3-Stroop1'          0.0002907    0.00075505    168        0.385       0.70072    [1x1 LinearModel]       0.98101
%     'WholeProbe:808'    'Stroop3-Stroop1:age'     1.2324e-06    0.00016504    168    0.0074672       0.99405    [1x1 LinearModel]       0.99405
%     'WholeProbe:808'    'Stroop3:age'            -0.00054567     0.0001814    168      -3.0082     0.0030326    [1x1 LinearModel]     0.0060652% 

% The "model" field in the table can then be used to show scatterplots, etc
% for the ROI values
disp(ROItable.model{1})
% Linear regression model:
%     beta ~ cond + cond:age
% 
% Estimated Coefficients:
%                               Estimate          SE          tStat        pValue  
%                              ___________    __________    _________    __________
% 
%     cond_Digit_Symbol          0.0028592     0.0008011       3.5692    0.00048322
%     cond_Digit_Symbol:age    -9.5464e-06    0.00018148    -0.052603       0.95812
% 
% 
% Number of observations: 150, Error degrees of freedom: 148
% Root Mean Squared Error: 1.03
% R-squared: 1.87e-05,  Adjusted R-Squared -0.00674
% F-statistic vs. constant model: 0.00277, p-value = 0.958

