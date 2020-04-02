function [ result ] = sec2string( time )
%SEC2STRING Summary of this function goes here
%   Detailed explanation goes here
hours = floor(time/3600);
minutes = floor(mod(time/60,60));
seconds = floor(mod(time,60));
result = [num2str(hours) ' hours and ' num2str(minutes) ' minutes and '...
    num2str(seconds) ' seconds'];

end

