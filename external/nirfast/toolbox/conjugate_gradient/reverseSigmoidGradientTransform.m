function [ zTrans ] = reverseSigmoidGradientTransform( z, x, xTrans )
%REVERSESIGMOIDTRANSFORM Summary of this function goes here
%   Detailed explanation goes here
zTrans = z.*(exp(xTrans)./((exp(xTrans) + 1).^2));

end

