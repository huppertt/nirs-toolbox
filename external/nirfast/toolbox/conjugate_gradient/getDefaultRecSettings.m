function [ recParam ] = getDefaultRecSettings( mesh, data, frequency )
%GETDEFAULTRECSETTINGS Creates default reconstruction parameter struct.
%   See comments.
%Maximum number of iterations for the algorithm to run.
recParam.numberIterations = 100;
%Whether or not to refine the initial line search estimate. Warning: this
%can dramatically increase run time.
recParam.refineQuadraticGuess = false;
%Whether or not to record the solution at each iteration of the algorithm.
%Warning: this can consume a significant quantity of memory.
recParam.recordX = false;
%Pixel basis to use in the reconstruction. Remove this field if pixel basis
%is unnecessary.
if mesh.dimension == 2
    recParam.pixelBasis = [30 30];
else
    recParam.pixelBasis = [30 30 30];
end
%Whether or not to use a hard prior (the mesh region labels) in the
%reconstruction.
recParam.useHardPrior = false;
%Percentage improvement in reconstruction at which to terminate the
%algorithm.
recParam.percentageImprovementCutoff = 1e-4;
%Number of iterations below percentage cutoff necessary to terminate
%reconstruction.
recParam.numberBelowCutoff = 2;
%The parameters to reconstruct. By default this is all the parameters in
%frequency domain, and mua or chromophores in continuous wave domain.
if strcmp(mesh.type,'stnd')
    recParam.reconstructionParameters = {'mua';'mus'};
    if frequency == 0
        recParam.reconstructionParameters =...
            recParam.reconstructionParameters(strcmp(recParam.reconstructionParameters,'mus')==0);
    end
elseif strcmp(mesh.type,'spec') || strcmp(mesh.type,'specPenn')
    recParam.reconstructionParameters = mesh.chromscattlist;
    if frequency == 0
        recParam.reconstructionParameters =...
            recParam.reconstructionParameters((strcmp(recParam.reconstructionParameters,'S-Amplitude') +...
            strcmp(recParam.reconstructionParameters,'S-Power')) == 0);
    end
end
%The coordinate spaces in which to perform the reconstruction.
for ii = 1:length(recParam.reconstructionParameters)
    [recParam.forwardTransformFunc{ii,1} recParam.reverseTransformFunc{ii,1} recParam.gradientTransformFunc{ii,1}] =...
        getTransformFunctions(recParam.reconstructionParameters{ii});
end
%Use wavelength dependent preconditioners by default.
recParam.useWavelengthDependentR = true;
%Recalculate preconditioner at least every 5 iterations.
recParam.preconditionerRecalculationInterval = 5;

end






