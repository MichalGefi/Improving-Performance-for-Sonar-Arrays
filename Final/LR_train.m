function [W] = LR_train(T,X,NB,lambda,ys)
%%%inputs
% T[NE,NT] - a set of NT arrays with NE element locations each
% X[NB*(NE+1),NT] - a set of NT arrays with NE element each mesured at NB different angles

%%% output
% W[NB*NE,NE] - trained weight matrix
PHI = base_function(X,NB,lambda,ys); % PHI[NB*NE,NT]
W = pinv(PHI.') * T.';

%NE = size(T,1);
%W2 = repmat(eye(NE),[NB,1])./NB; %wiegth for AVg method

end

