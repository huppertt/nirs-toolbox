function [ xNew ] = invSigmoid( x )
%INVSIGMOID Summary of this function goes here
%   Detailed explanation goes here
xNew = -log(1./x - 1);

end

