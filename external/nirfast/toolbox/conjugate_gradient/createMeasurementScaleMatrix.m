function [ P currentData ] = createMeasurementScaleMatrix(mesh, measurements, frequency, currentData)
%CREATEMEASUREMENTSCALEMATRIX Creates precision matrix for measurements to
%normalise contribution of amplitude and phase.
%   Returns precision matrix P and forward model data. Precision matrix is
%   returned in sparse format.

if ~isfield(mesh,'cgFormatLink')
    mesh = convertLinkFormat(mesh);
end

%Calculate contribution from amplitude and phase.
if nargin < 4
    currentData = calculateForwardData(mesh,frequency);
end
    
    %%
    
if strcmp(mesh.type,'stnd')
    ampInd = 1;
    phasInd = 2;
    measAmp = measurements.paa((mesh.cgFormatLink(:,3) == 1),ampInd);
    datAmp = currentData.paa((mesh.cgFormatLink(:,3) == 1),ampInd);
    measPhas = measurements.paa((mesh.cgFormatLink(:,3) == 1),phasInd);
    datPhas = currentData.paa((mesh.cgFormatLink(:,3) == 1),phasInd);
    logAmpDiff = log(measAmp ./ datAmp);
    phaseDiff = pi*(measPhas - datPhas)/180;
    %Make phase diff periodic.
    phaseDiff(phaseDiff > pi) = phaseDiff(phaseDiff > pi) - 2*pi;
    phaseDiff(phaseDiff < -pi) = phaseDiff(phaseDiff < -pi)+ 2*pi;
    if frequency > 0
        b = [ones(size(logAmpDiff))*norm(logAmpDiff); ones(size(phaseDiff))*norm(phaseDiff)];
    else
        b = [ones(size(logAmpDiff))*norm(logAmpDiff); ones(size(phaseDiff))];
    end
    P = [(1:length(b))',(1:length(b))',1./(b.^2)];
elseif strcmp(mesh.type,'spec')
    P = cell(length(mesh.wv),1);
    for ii = 1:length(mesh.wv)
        ampInd = 2*ii - 1;
        phasInd = 2*ii;
        measAmp = measurements.paa((mesh.cgFormatLink(:,2+ii) == 1),ampInd);
        datAmp = currentData.paa((mesh.cgFormatLink(:,2+ii) == 1),ampInd);
        measPhas = measurements.paa((mesh.cgFormatLink(:,2+ii) == 1),phasInd);
        datPhas = currentData.paa((mesh.cgFormatLink(:,2+ii) == 1),phasInd);
        logAmpDiff = log(measAmp ./ datAmp);
        phaseDiff = pi*(measPhas - datPhas)/180;
        %Make phase diff periodic.
        phaseDiff(phaseDiff > pi) = phaseDiff(phaseDiff > pi) - 2*pi;
        phaseDiff(phaseDiff < -pi) = phaseDiff(phaseDiff < -pi)+ 2*pi;
        if frequency > 0
            b = [ones(size(logAmpDiff))*norm(logAmpDiff); ones(size(phaseDiff))*norm(phaseDiff)];
        else
            b = [ones(size(logAmpDiff))*norm(logAmpDiff); ones(size(phaseDiff))];
        end
        P{ii} = [(1:length(b))',(1:length(b))',1./(b.^2)];
    end
end
end

