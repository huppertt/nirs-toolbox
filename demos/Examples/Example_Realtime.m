%This example will show how to use the real-time plugins in the NIRS-toolbox

% The code uses a realtime version of the core.Data class
raw=nirs.realtime.core.Data;

%This class has the following methods (in addition to everything the static data class has
% 
% reset  - clears any data and stimulus events           
% adddata - adds data and time to class
% addevent  - add an impulse stim mark      
% addeventStart - start a block (duration) event
% addeventEnd - close an open (started) event

% The data class has a field called
%  updatefunction: {} - this specifies the pipeline to run in real time.
% the data class is responsible for its own real-time processing

% a listener passes data from a source to stores it in data
MyListener = nirs.realtime.listeners.simulator;
% Th simulator listener takes data from a source (e.g.
% @nirs.testing.SimData) and moves it to the data_output

MyListener.datasource = @nirs.testing.simData;  
% the data source can either be a function handle like above or an existing
% data class (e.g load existing data and then rerun a pipeline)

% the third part of the realtime code is to preform analysis using job
% pipelines created using nirs.realtime.modules.*
% this is similar to static module pipelines
myRTjob = nirs.realtime.modules.MotionCorrect;
myRTjob = nirs.realtime.modules.BandPass(myRTjob);

% this is assigned to the realtime data class's update function
raw.updatefunction = myRTjob;

% finally the data class is assigned to the output of the listener
MyListener.data_output=raw;

% then the listener is started and stopped like this;
MyListener.start;
pause(10);
MyListener.stop;

% if we look at raw (on the command line), we see that the data magically
% appeared here.  We use handles (which are like Matlab's version of
% passing by reference) so that any object that is holding onto a reference
% of raw is also updated.

% if we do the same thing but then plot the data
MyListener.start;
raw.gui;  % draws the interactive plot.  raw.draw would do the same

% note- you draw first and then start the listener, but the raw.probe field 
% is often created in the initialization of the listener (depending on the data source) 
% so starting then drawing is more reliable 

pause(10);
MyListener.stop;





