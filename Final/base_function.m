function [PHI] = base_function(x,NB,lambda,ys)
k = 2*pi/lambda; % of wave
NE = (size(x,1)-NB)/NB;
if isvector(x) %base function used for LR_place_estimator (phi(x))
    % takes X[NB*(NE+1),1] and returns PHI[NB*NE,1]
    beta = zeros(NB,1);
    meas = zeros(NE*NB,1);
    for j = 1:NB
        beta(j) = x(1+(j-1)*(NE+1));
        meas(1+(j-1)*NE:j*NE) = x(2+(j-1)*(NE+1):j*(NE+1));
        meas(1+(j-1)*NE:j*NE) = meas(1+(j-1)*NE:j*NE).*exp(1i*k.*cosd(beta(j)).*ys);
    end
    rep_beta = repelem(beta,NE);

else %base function used for LR_train (PHI(X))
    % takes X[NB*(NE+1),NT] and returns PHI[NB*NE,NT]
    NT = size(x,2);
    beta = zeros(NB,NT);
    meas = zeros(NE*NB,NT);
    for n = 1:NT
        for m = 1:NB
            beta(m,n) = x(1+(m-1)*(NE+1),n);
            meas(1+(m-1)*NE:m*NE,n) = x(2+(m-1)*(NE+1):m*(NE+1),n);
            meas(1+(m-1)*NE:m*NE,n) = meas(1+(m-1)*NE:m*NE,n).*exp(1i*k.*cosd(beta(m,n)).*ys);
        end
    end
    rep_beta = zeros(NE*NB,NT);
    for n = 1:NT
        rep_beta(:,n) = repelem(beta(:,n),NE);
    end

end
    PHI = -1*angle(meas)./(k.*sind(rep_beta));
end