function out = calibrate_data( data, calib, model )
%CALIBRATE_DATA Summary of this function goes here
%   Detailed explanation goes here

    c = mean( log( calib.data ),1 ); % calibration measurement
    m = mean( log( model.data ),1 ); % model measurement
    
    x = m - c; % calibration term
    
    d = log( data.data ); % log data
    
    d = bsxfun(@plus,d,x); % calibrated log data
    
    out = data; % output
    out.data = exp( d ); % calibrate data
    
end

