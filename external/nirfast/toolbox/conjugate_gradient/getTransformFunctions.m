function [forwardFunc reverseFunc gradientFunc] = getTransformFunctions(parameterName)
%GETTRANSFORMFUNCTIONS Returns default transformation functions for
%parameter.
%   Returns log transformations for parameters in the interval
%   (0,infinity], reverse Sigmoid transformations for parameters in the
%   interval (0,1), and dummy transformations (no transformation) for
%   anything unknown.
if strcmp(parameterName,'mua')
    space = 'log';
elseif strcmp(parameterName,'mus')
    space = 'log';
elseif strcmp(parameterName,'HbO')
    space = 'log';
elseif strcmp(parameterName,'deoxyHb')
    space = 'log';
elseif strcmp(parameterName,'Water')
    space = 'reverseSigmoid'; %Log transformation appears to reconstruct better, but handles constraints suboptimally.
elseif strcmp(parameterName,'Fat')
    space = 'reverseSigmoid'; %Log transformation appears to reconstruct better, but handles constraints suboptimally.
elseif strcmp(parameterName,'S-Amplitude')
    space = 'log';
elseif strcmp(parameterName,'S-Power')
    space = 'log';
else
    space = 'untransformed';
end
[forwardFunc reverseFunc gradientFunc] = getTransformsForSpace(space);
end

function [forwardFunc reverseFunc gradientFunc] = getTransformsForSpace(space)
if strcmp(space,'log')
    %Use 'log' for parameters belonging to the interval (0,infinity].
    forwardFunc = @log;
    reverseFunc = @exp;
    gradientFunc = @logGradientTransform;
elseif strcmp(space,'reverseSigmoid')
    %Use 'reverseSigmoid' for parameters belonging to the interval (0,1).
    forwardFunc = @invSigmoid;
    reverseFunc = @sigmoid;
    gradientFunc = @reverseSigmoidGradientTransform;
elseif strcmp(space,'untransformed');
    %Use 'untransformed' for unconstrained parameters.
    forwardFunc = @dummyFunction;
    reverseFunc = @dummyFunction;
    gradientFunc = @dummyFunction;
else
    forwardFunc = @dummyFunction;
    reverseFunc = @dummyFunction;
    gradientFunc = @dummyFunction;
end
end

