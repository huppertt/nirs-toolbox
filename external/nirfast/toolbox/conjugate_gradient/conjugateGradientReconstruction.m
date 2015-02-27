function [mesh xBestVal objVal xSeries totalTime gradientTime...
    lineSearchTime] = conjugateGradientReconstruction...
    (mesh, measurements, frequency, recParam, regParam, output_fn)
%CONJUGATEGRADIENTRECONSTRUCTION Reconstructs parameters using
%gradient-based method.
%   Reconstructs parameter maps from continuous wave or frequency domain
%   measurements for standard or spectral problems.

%% Initialise variables.

%Convert link file if necessary.
mesh.link = measurements.link;
mesh = convertLinkFormat(mesh);

%If no reconstruction settings are provided, get default settings.
if nargin < 4
    recParam = getDefaultRecSettings(mesh,measurements,frequency);
end

%If no regularisation settings are provided, get default settings.
if nargin < 5
    [regParam recParam mesh] = getDefaultRegSettings(mesh,measurements,frequency,recParam);
end

% Initialise x (parameter) value matrix to 0 if not recording algorithm path.
if ~recParam.recordX
    xSeries = 0;
end

%Create recon pixel mesh if appropriate.
if isfield(recParam,'pixelBasis')
    if isstruct(recParam.pixelBasis)
        reconMesh = recParam.pixelBasis;
    else
        [mesh reconMesh] = createPixelMesh(mesh,recParam.pixelBasis);
    end
end

%Precalculate preconditioner.
if ~isfield(mesh,'R')
    mesh.R = generateRMatrix(mesh,frequency,recParam);
end

%Keep track of last time preconditioner was calculated.
iterationsSincePreconditionerCalculation = 0;

%Initialise time recording variables.
totalTimeStart = tic;
lineSearchTime = 0;
gradientTime = 0;

%% Initialise parameter vector.

%Create parameter (optical or chromophore properties) vector.
x = constructParameterVector(mesh,recParam);
%Transform parameter vector into appropriate coordinate space for
%reconstruction.
x = forwardTransformParameterVector(x,mesh,recParam);

% Store initial parameter vector for use in regularisation (this is used in
% (x - x0)'Q(x-x0) term).
regParam.x0 = x;

%% Calculate gradient and calculate solution fitness.

timeStart = tic;
% try
    if exist('reconMesh','var')
        [z currentData] = calculateGradient(mesh, measurements,frequency,recParam,regParam,reconMesh);
    else
        [z currentData] = calculateGradient(mesh, measurements,frequency,recParam,regParam);
    end
% catch exObj
%     fprintf([exObj.message '\n']);
%     error('Terminating because of inability to calculate gradient.');
% end
gradientTime = gradientTime + toc(timeStart);

fprintf('Iteration 1\n');
fprintf('Initial z norm = %e\n',...
    norm(z));

%Create objective function value storage.
objVal = zeros(recParam.numberIterations+1,1);
%Calculate solution fitness value, and fractional sizes of log amplitude
%and phase errors.
[objVal(1) logAmpFrac phaseFrac] =...
    calculateObjectiveFunction(mesh,measurements,frequency,currentData,recParam,regParam);

fprintf('Initial percentage projection error = %e\n',...
    objVal(1));
fprintf('Log amplitude error is %e percent, and phase error is %e percent.\n',...
    logAmpFrac*100,phaseFrac*100);
xBest = x;
xBestVal = objVal(1);

if recParam.recordX
    xSeries = zeros(recParam.numberIterations,length(x));
    xSeries(1,:) = reverseTransformParameterVector(x,mesh,recParam);
end

%% Initialise conjugate gradient variables.

d = -z;
p = d;
n = 1;
numberBelowCutoff = 0;
% Initial alpha of 0.1 is not optimal. Error catching should deal with the
% case where the movement is so large that it makes the MASS matrix badly
% conditioned, but this may be changed in future.
alpha = 1e-1;
%% Loop until conditions satisfied.
while(n <= recParam.numberIterations)
    timeStart = tic;
    %% Perform line search.
    lineSearchCompleted = false;
    attempt = 1;
    maxAttempt = 6;
    % If step size is so large it makes the MASS matrix badly
    % conditioned, reduce it until it is small enough not to. Scaling by
    % ratio of gradients may also work, but has been observed to lead to
    % inferior final result.
    while (~lineSearchCompleted) && (attempt <= maxAttempt)
        try
            alpha = quadraticLineSearch(mesh,measurements,frequency,x,d,alpha,objVal(n),recParam,regParam);
            lineSearchCompleted = true;
        catch exObj
            alpha = alpha/2^5;
            fprintf([exObj.message '\n']);
            fprintf('Reducing alpha because error occurred in line search. New alpha is %e.\n',alpha);
            attempt = attempt + 1;
            %If line search step fails and preconditioner is old,
            %recalculate it.
            if iterationsSincePreconditionerCalculation > 0
                fprintf('Recalculating preconditioner to attempt to reduce numerical instability.\n');
                mesh = rmfield(mesh,'R');
                mesh.R = generateRMatrix(mesh,frequency,recParam);
                iterationsSincePreconditionerCalculation = 0;
            end
        end
    end
    if attempt > maxAttempt
        disp('Terminating because of numerical error');
        break;
    end
    if isnan(alpha)
        disp('Terminating because alpha is NaN');
        break;
    end
    if isinf(alpha)
        disp('Terminating because alpha is inf');
        break;
    end
    lineSearchTime = lineSearchTime + toc(timeStart);
    clear xNew
    
    %% Update parameters.
    xNew = x + alpha.*d;
    
    %Ensure parameters are physically valid. If correct parameter
    %transformations are used, this should not occur.
    xNew = reverseTransformParameterVector(xNew,mesh,recParam);
    xNew = sanitiseXNew(mesh,xNew,recParam);
    
    mesh = changeMeshParameters(mesh,xNew,recParam);
    xNew = forwardTransformParameterVector(xNew,mesh,recParam);
    
    iterationsSincePreconditionerCalculation = iterationsSincePreconditionerCalculation + 1;
    
    %Update preconditioner if it is sufficiently old.
    
    if isfield(recParam,'preconditionerRecalculationInterval') && (iterationsSincePreconditionerCalculation >= recParam.preconditionerRecalculationInterval)
        fprintf('Recalculating preconditioner\n');
        mesh = rmfield(mesh,'R');
        mesh.R = generateRMatrix(mesh,frequency,recParam);
        iterationsSincePreconditionerCalculation = 0;
    end
    
    
    %% Update regularisation weight if appropriate.
    if regParam.lambda > 0
        regParam.lambda = regParam.lambda/regParam.denom;
    end
    
    %% If maximum number of iterations has been reached, calculate solution
    %% fitness and choose best solution to output.
    if n == recParam.numberIterations
        [objVal(n+1) logAmpFrac phaseFrac] =...
            calculateObjectiveFunction(mesh,measurements,frequency,recParam,regParam);
        if(objVal(n+1)<xBestVal)
            xBest = xNew;
            xBestVal = objVal(n+1);
        end
        x = xNew;
        if recParam.recordX
            xSeries(n+1,:) = reverseTransformParameterVector(x,mesh,recParam);
        end
        break;
    end
    
    %% Calculate new gradient and solution fitness.
    
    if exist('zNew','var')
        clear zNew
    end
    clear currentData
    timeStart = tic;
    
    %Calculate new gradient and fitness.
    try
        if exist('reconMesh','var')
            [zNew currentData] = calculateGradient(mesh, measurements,frequency,recParam,regParam,reconMesh);
        else
            [zNew currentData] = calculateGradient(mesh, measurements,frequency,recParam,regParam);
        end
    catch exObj
        %If gradient calculation fails and preconditioner is old,
        %recalculate it.
        if iterationsSincePreconditionerCalculation > 0
            try
                fprintf('Recalculating preconditioner to attempt to reduce numerical instability.\n');
                mesh = rmfield(mesh,'R');
                mesh.R = generateRMatrix(mesh,frequency,recParam);
                iterationsSincePreconditionerCalculation = 0;
                if exist('reconMesh','var')
                    [zNew currentData] = calculateGradient(mesh, measurements,frequency,recParam,regParam,reconMesh);
                else
                    [zNew currentData] = calculateGradient(mesh, measurements,frequency,recParam,regParam);
                end
            catch exObj
                fprintf([exObj.message '\n']);
                fprintf('Terminating because of numerical instability.');
                break;
            end
        else
            fprintf([exObj.message '\n']);
            fprintf('Terminating because of numerical instability.');
            break;
        end
    end
    gradientTime = gradientTime + toc(timeStart);
    fprintf('Iteration %i\n',n+1);
    fprintf('New z norm = %e\n',...
        norm(zNew));
    
    [objVal(n+1) logAmpFrac phaseFrac] =...
        calculateObjectiveFunction(mesh,measurements,frequency,currentData,recParam,regParam);%Calculate regularised objective value
    
    fprintf('phi = %e. %e percent improvement.\n',...
        objVal(n+1),100*(1-objVal(n+1)/objVal(n)));
    fprintf('Log amplitude error is %e percent, and phase error is %e percent.\n',...
        logAmpFrac*100,phaseFrac*100);
    if(objVal(n+1)<xBestVal)
        xBest = xNew;
        xBestVal = objVal(n+1);
    end
    
    %% Terminate algorithm if termination conditions are satisfied.
    
    if ((objVal(n+1)-objVal(n))/objVal(n)) > -(recParam.percentageImprovementCutoff/100)
        numberBelowCutoff = numberBelowCutoff + 1;
        if numberBelowCutoff >= recParam.numberBelowCutoff
            disp('break due to small improvement');
            break;
        end
    else
        numberBelowCutoff = 0;
    end
    if (objVal(n+1) > objVal(n)) && (iterationsSincePreconditionerCalculation > 0)
        fprintf('Recalculating preconditioner to attempt to reduce numerical instability.\n');
        mesh = rmfield(mesh,'R');
        mesh.R = generateRMatrix(mesh,frequency,recParam);
        iterationsSincePreconditionerCalculation = 0;
    end
    if isfield(recParam,'minLogAmpFrac') && isfield(recParam,'minPhaseFrac')
        if (logAmpFrac < recParam.minLogAmpFrac) && (phaseFrac < recParam.minPhaseFrac)
            disp('break due to fractional measurement error threshold reached');
            break;
        end
    elseif isfield(recParam,'minLogAmpFrac')
        if logAmpFrac < recParam.minLogAmpFrac
            disp('break due to fractional log amplitude error threshold reached');
            break;
        end
    elseif isfield(recParam,'minPhaseFrac')
        if phaseFrac < recParam.minPhaseFrac
            disp('break due to fractional phase error threshold reached');
            break;
        end
    end
    
    %% Update conjugate gradient variables and repeat.
    pNew = -zNew;
    beta = max(0,pNew'*(pNew-p)/(p'*p));
    fprintf('beta = %e\n',beta);
    dNew = pNew + beta*d;
    
    x = xNew;
    if recParam.recordX
        xSeries(n+1,:) = reverseTransformParameterVector(x,mesh,recParam);
    end
    p = pNew;
    d = dNew;
    
    %% write to solution file
    mesh_out = changeMeshParameters(mesh,reverseTransformParameterVector(x,mesh,recParam),recParam);
    if strcmp(mesh_out.type,'stnd')
          if n == 1
            fid = fopen([output_fn '_mua.sol'],'w');
          else
            fid = fopen([output_fn '_mua.sol'],'a');
          end
          fprintf(fid,'solution %g ',n);
          fprintf(fid,'-size=%g ',length(mesh_out.nodes));
          fprintf(fid,'-components=1 ');
          fprintf(fid,'-type=nodal\n');
          fprintf(fid,'%f ',mesh_out.mua);
          fprintf(fid,'\n');
          fclose(fid);

          if n == 1
            fid = fopen([output_fn '_mus.sol'],'w');
          else
            fid = fopen([output_fn '_mus.sol'],'a');
          end
          fprintf(fid,'solution %g ',n);
          fprintf(fid,'-size=%g ',length(mesh_out.nodes));
          fprintf(fid,'-components=1 ');
          fprintf(fid,'-type=nodal\n');
          fprintf(fid,'%f ',mesh_out.mus);
          fprintf(fid,'\n');
          fclose(fid);
    else
        all_sol = char(mesh_out.chromscattlist);
        [n_allsol,junk] = size(all_sol);
        chrom_scatt = mesh_out.conc;
        chrom_scatt(:,end+1) = mesh_out.sa;
        chrom_scatt(:,end+1) = mesh_out.sp;
        for j = 1:n_allsol
            tmp_string = [strcat(output_fn,'_',all_sol(j,:),'.sol')];
            if (n==1)
                fid = fopen(tmp_string,'w');
            else
                fid = fopen(tmp_string,'a');
            end
            fprintf(fid,'solution %d ',n);
            fprintf(fid,'-size=%g ',length(mesh_out.nodes));
            fprintf(fid,'-components=1 ');
            fprintf(fid,'-type=nodal\n');
            fprintf(fid,'%f ',chrom_scatt(:,j));
            fprintf(fid,'\n');
            fclose(fid);
            clear tmp_string
        end
    end
    
    n = n + 1;
    
end
%% Update mesh to best solution.
mesh = changeMeshParameters(mesh,reverseTransformParameterVector(xBest,mesh,recParam),recParam);
totalTime = toc(totalTimeStart);
disp(['Gradient calculations took ' sec2string(gradientTime) '.']);
disp(['Line search calculations took ' sec2string(lineSearchTime) '.']);
disp(['The reconstruction process took ' sec2string(totalTime) '.']);
end

function [z currentData] = calculateGradient(mesh, measurements, frequency, recParam, regParam, reconMesh)

%Calculate unregularised, untransformed gradient.
[z currentData] = adjointGradient(mesh, measurements,frequency,recParam,regParam);

%Regularise gradient.
if regParam.lambda > 0
    x = constructParameterVector(mesh,recParam);
    b = x - reverseTransformParameterVector(regParam.x0,mesh,recParam);
    zReg = regParam.lambda*(regParam.qMat.^0.5)*b/2;
    z = z + zReg;
end

%Transform gradient.
z = transformGradient(z,mesh,recParam);

% Smooth gradient using recon mesh if appropriate.
if nargin >= 6
    z = pixelBasisSmooth(mesh,z,reconMesh);
end

% Homogenise gradient using hard prior if appropriate.
if recParam.useHardPrior
    z = applyHardPrior(mesh,z);
end
end


