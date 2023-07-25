function [xs] = AVG_place_estimator(data,NB,lambda,ys)
% takes X[NB*(NE+1),1] and returns xs[NE,1]
% requiesrs global parameter ys
k=2*pi/lambda;
NE = (size(data,1)-NB)/NB;

beta = zeros(NB,1);
meas = zeros(NE*NB,1);
for j = 1:NB
    beta(j) = data(1+(j-1)*(NE+1));
    meas(1+(j-1)*NE:j*NE) = data(2+(j-1)*(NE+1):j*(NE+1));
    meas(1+(j-1)*NE:j*NE) = meas(1+(j-1)*NE:j*NE).*exp(1i*k.*cosd(beta(j)).*ys);
end

rep_beta = repelem(beta,NE);
PHI = -1*angle(meas)./(k.*sind(rep_beta));

xs = zeros(NE,1);
for j = 1:NB
    xs = xs + PHI(1+(j-1)*NE:j*NE);
end

xs = xs./NB;
end

