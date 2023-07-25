function [out_est,out_angle] = phased_array(lambda,xs,D,meas,N)
k=2*pi/lambda;
NE = length(xs);
ys = (-D*((NE-1)/2):D:D*((NE-1)/2)).';
boundery = 1;
theta = linspace(boundery,180-boundery,N); % posible estimated angles

% estimated measurments of the plane wave field
phi_theta = zeros(N,NE); % Matrix of estimated measurments for each posible angle
for j = 1:length(theta)
    phi_theta(j,:) = exp(-1i*(k.*(xs*sind(theta(j))+ys*cosd(theta(j)))));
end

% applying coefficients to each element (all 1 for phased array)
coef = ones(N,NE);
for j = 1:length(theta)
    phi_theta(j,:) = coef(j,:).*phi_theta(j,:);%Horizontal vector of elemnts field estimated measurements
end

% applying estimator
out_est = zeros(N,1);
for j = 1:length(theta)
    out_est(j,:) = abs((conj(phi_theta(j,:))*meas)./NE);
end
out_est = out_est./trapz(theta,out_est); % normalization


% Output %
max_out_vals = [max(out_est(1:round(length(out_est)/2))),max(out_est(round(length(out_est)/2):end))];
indexs_of_max_valas = [find(out_est(1:round(length(out_est)/2))>=max_out_vals(1)), (length(out_est)/2 + find(out_est(round(length(out_est)/2):end)>=max_out_vals(2)))] ;
out_angle = theta(indexs_of_max_valas(1)); % selecting the correct angle