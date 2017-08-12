function [phi logAmpFrac phaseFrac] = calculateObjectiveFunction...
    (mesh, measurements, frequency, currentData,recParam,regParam)
%CALCULATEOBJECTIVEFUNCTION Calculates the value of the objective function
%for the forward data and measurements.
%   Runs forward model if forward data is not provided.
%If forward model data is not provided calculate it, and if regularisation
%is not provided assume none.
if ~isfield(mesh,'cgFormatLink')
    mesh = convertLinkFormat(mesh);
end
if nargin < 5
    currentData = calculateForwardData(mesh,frequency);
    regParam.lambda = 0;
elseif nargin < 6
    if ~ isfield(currentData,'paa')
        regParam = recParam;
        recParam = currentData;
        currentData = calculateForwardData(mesh,frequency);
    else
        regParam.lambda = 0;
    end
end
if ~isfield(regParam,'P')
    regParam.P = 1;
end
%Calculate residual vectors.
ampInd = 1:2:size(measurements.paa,2);
phasInd = 2:2:size(measurements.paa,2);
measAmp = measurements.paa(:,ampInd);
datAmp = currentData.paa(:,ampInd);
logAmpDiff = log(measAmp ./ datAmp);
measPhas = measurements.paa(:,phasInd);
datPhas = currentData.paa(:,phasInd);
phaseDiff = pi*(measPhas - datPhas)/180;
%Account for numerical error in phase calculation in continuous wave mode.
if frequency == 0
    phaseDiff(:) = 0;
end
%Make phase difference periodic (ie. phase difference lies in the interval
%[-pi,pi].
phaseDiff(phaseDiff > pi) = phaseDiff(phaseDiff > pi) - 2*pi;
phaseDiff(phaseDiff < -pi) = phaseDiff(phaseDiff < -pi)+ 2*pi;

%Calculate the fractional residual vector.
logAmpFrac = norm(flatten(logAmpDiff(mesh.cgFormatLink(:,3:end) == 1)))/...
    norm(flatten(log(measAmp(mesh.cgFormatLink(:,3:end) == 1))));
phaseFrac = norm(flatten(phaseDiff(mesh.cgFormatLink(:,3:end) == 1)))/...
    norm(flatten(pi*measPhas(mesh.cgFormatLink(:,3:end) == 1)/180));
%Calculate the normalised sum of squares error.
if strcmp(mesh.type,'stnd')
    b = [logAmpDiff(mesh.cgFormatLink(:,3) == 1); phaseDiff(mesh.cgFormatLink(:,3) == 1)];
    if isscalar(regParam.P)
        pMat = regParam.P;
    else
        pMat = sparse(regParam.P(:,1),regParam.P(:,2),regParam.P(:,3));
    end
    phi = b'*pMat*b/2;
elseif strcmp(mesh.type,'spec') || strcmp(mesh.type,'specPenn')
    phi = 0;
    for ii = 1:length(mesh.wv)
        b = [logAmpDiff((mesh.cgFormatLink(:,2+ii) == 1),ii);phaseDiff((mesh.cgFormatLink(:,2+ii) == 1),ii)];
        if isscalar(regParam.P)
            pMat = regParam.P;
        else
            currentP = regParam.P{ii};
            pMat = sparse(currentP(:,1),currentP(:,2),currentP(:,3));
        end
        phi = phi + b'*pMat*b/2;
    end
end


%%

% Add regularisation term to sum of squares error.
if regParam.lambda > 0
    x = constructParameterVector(mesh,recParam);
    
    b = x - reverseTransformParameterVector(regParam.x0,mesh,recParam);
    phi = phi + regParam.lambda*b'*regParam.qMat*b;
end

end

function [flattenedMatrix] = flatten(matrix)
flattenedMatrix = matrix(:);
end