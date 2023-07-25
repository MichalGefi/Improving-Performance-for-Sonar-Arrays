function [xs] = LR_place_estimator(data,W,NB,lambda,ys)
% using trained weights w and applying LR
xs = (W.')*base_function(data,NB,lambda,ys);
end

