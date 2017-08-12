function [ xNew ] = sigmoid( x )
%SIGMOID Summary of this function goes here
%   Detailed explanation goes here
xNew = 1./(1 + exp(-x));

end

