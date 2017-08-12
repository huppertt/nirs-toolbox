function [ alpha ] = quadraticLineSearch( mesh,...
    measurements, frequency, x, d, s1, originalFitness, recParam,regParam)
%QUADRATICLINESEARCH Calculates approximate line minimum using quadratic
%fitting.
%   Adapted from http://www-personal.umich.edu/~murty/611/611slides9.pdf

%Find three points s0, sm, s1 such that s0 < sm < s1 and f(s0) > s(sm) and
%f(s1) > f(sm).

%% If no maximum number of forward model evaluations is given, use a
%% default of 50.
if isfield(recParam,'maxFEM')
    maxFEM = recParam.maxFEM;
else
    maxFEM = 50;
end
forwardSolves = 0;
% s0, sm and s1 are step sizes such that s0 < sm < s1.
s0 = 0;
f_s0 = originalFitness;
% f_s0, f_sm, and f_s1 are the fitness values at s0, sm, and s1
% respectively.
fprintf('Beginning line search\n');
f_s1 = calculateObjMov(mesh, measurements, frequency,x,d,s1,recParam,regParam);
forwardSolves = forwardSolves + 1;
%If f_s1 is less than f_s0 widen the search until f_s1 is greater than
%f_sm.
if f_s1 < f_s0
    sm = s1;
    s1 = 2*s1;
    f_sm = f_s1;
    f_s1 = calculateObjMov(mesh, measurements, frequency,x, d, s1,recParam,regParam);
    forwardSolves = forwardSolves + 1;
    while f_s1 < f_sm
        s0 = sm;
        sm = s1;
        s1 = 2*s1;
        f_s0 = f_sm;
        f_sm = f_s1;
        f_s1 = calculateObjMov(mesh, measurements, frequency,x, d, s1,recParam,regParam);
        forwardSolves = forwardSolves + 1;
        %Terminate early if maximum FEM evaluations reached.
        if forwardSolves > maxFEM
            break
        end
    end
else
    %If f_s1 is greater than or equal to f_s0 shrink the search until f_sm
    %is greater than f_s0.
    sm = s1/2;
    f_sm = calculateObjMov(mesh, measurements, frequency,x, d, sm,recParam,regParam);
    forwardSolves = forwardSolves + 1;
    while f_sm > f_s0
        s1 = sm;
        sm = sm/2;
        f_s1 = f_sm;
        f_sm = calculateObjMov(mesh, measurements, frequency,x, d, sm,recParam,regParam);
        forwardSolves = forwardSolves + 1;
        %Terminate early if maximum FEM evaluations reached.
        if forwardSolves > maxFEM
            break
        end
    end
end
%Fit quadratic and calculate minimum.
alpha = calculateMinimum(s0,sm,s1,f_s0,f_sm,f_s1);
% fprintf('s0 is %f, sm is %f, s1 is %f, and alpha is %f\n',s0,sm,s1,alpha);
% fprintf('Initial fitness was %f, f_s0 is %f, f_sm is %f, f_s1 is %f, and f_alpha is %f\n',originalFitness,f_s0,f_sm,f_s1,calculateObjMov(mesh, measurements, frequency,x, d, alpha,recParam,regParam));


%Improve accuracy of guess by re-estimating using interval obtained from
%three out of four points s0, sm, s1, and alpha.
breakCriteria = max(1e-9*alpha,realmin*1e2);
while ((s1 - s0) > breakCriteria) && isfield(recParam,'refineQuadraticGuess') && recParam.refineQuadraticGuess
    if forwardSolves > maxFEM
        break
    end
    f_alpha = calculateObjMov(mesh, measurements, frequency,x, d, alpha ,recParam,regParam);
    forwardSolves = forwardSolves + 1;
    if alpha < sm
        if f_alpha < f_sm
            s1 = sm;
            f_s1 = f_sm;
            sm = alpha;
            f_sm = f_alpha;
        elseif f_alpha > f_sm
            s0 = alpha;
            f_s0 = f_alpha;            
        else
            break;
        end
    elseif alpha > sm
        if f_alpha < f_sm
            s0 = sm;
            f_s0 = f_sm;
            sm = alpha;
            f_sm = f_alpha;
        elseif f_alpha > f_sm  
            s1 = alpha;
            f_s1 = f_alpha;
        else
            break;
        end
    else
        break;
    end
    oldAlpha = alpha;
    alpha = calculateMinimum(s0,sm,s1,f_s0,f_sm,f_s1);
    if isinf(alpha) || isnan(alpha) || (alpha < s0) || (alpha > s1)
        alpha = oldAlpha;
        break;
    end
end
fprintf('%i fem calculations used in line search.\n',forwardSolves);
end

function value = calculateObjMov(mesh, measurements, frequency,x, d , s,recParam,regParam)
%CALCULATEOBJMOV Alters mesh properties as required and calculates
%objective function value.
%
xNew = x + s*d;
xNew = reverseTransformParameterVector(xNew,mesh,recParam);
xNew = sanitiseXNew(mesh,xNew,recParam);
testMesh = changeMeshParameters(mesh,xNew,recParam);
value = calculateObjectiveFunction(testMesh,measurements,frequency,recParam,regParam);
clear testMesh;
end

function alpha = calculateMinimum(s0,sm,s1,f_s0,f_sm,f_s1)
%CALCULATEMINIMUM Fits quadratic to three points and calculates minimum.
%   
a = (s1*f_sm+s0*(f_s1-f_sm)+(sm-s1)*f_s0-f_s1*sm);%/denom;
b = -(s1^2*f_sm+s0^2*(f_s1-f_sm)+(sm^2-s1^2)*f_s0-f_s1*sm^2);%/denom;
alpha = -b/(2*a);
end
