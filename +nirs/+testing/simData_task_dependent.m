function [data, truth] = simData_task_dependent(probe, noise, stim, beta, channels, basis )
%SIMDATA Simulates NIRS data by adding a task to baseline noise.

time=(0:1/10:300)';

if nargin < 3 || isempty(stim)
    stim = nirs.testing.randStimDesign(time, 2, 7, 1);
end

if(isa(stim,'function_handle'))
    stim = stim(time);
end


if nargin < 2 || isempty(noise)  
    noise = [1 2*(ones(stim.count,1))];
end

if nargin<1
    probe=[];
end

noise = nirs.testing.simARNoise_task_dependent(probe,time,[],[],stim,noise);
    
  

if nargin < 4 || isempty(beta)
    beta = [];
end
if nargin < 5 
    channels=[];
end
if nargin < 6
    basis=[];
end


[data, truth] = nirs.testing.simData( noise, stim, beta, channels, basis );
