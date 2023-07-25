function [T,X] = create_data(NE,lambda,D,NT,beta,dev_sigma,dev_max,curve_sigma,curve_max,meas_noise)

%%% inputs
% NE- number of elemnts
% lambda- lambda of wave
% D- distance between elements (y axis)
% NT- number of training sets to generate
% beta- array of angle for incoming waves used for calibration
% sigma- std of deveation of locations
% curve_max- max coef of array curve

%%% outputs
% T- generated true locations (xs_dev)
% X- wave mesurments at true element locations

k=2*pi/lambda; % of wave
NB = length(beta);


%%% Element locations setup %%%
ys = (-D*((NE-1)/2):D:D*((NE-1)/2)).';

xs_true = zeros(NE, NT);
% Adding loaction curve and deviation
curve = truncnormrnd([1,NT], 0, curve_sigma, -curve_max, curve_max);
dev = truncnormrnd([NE,NT],0, dev_sigma, -dev_max, dev_max);
xs_true = xs_true + dev ;
for m = 1:NT
    xs_true(:,m) = xs_true(:,m) + ((ys).^2) .*curve(m)./(NE*D);% add curve
    xs_true(:,m) = xs_true(:,m)-mean(xs_true(:,m));
end
    
T = xs_true;

data = zeros(NB*(NE+1),NT);
for n=1:NT
    for m = 1:NB
        meas = exp(-1i*(k.*(xs_true(:,n)*sind(beta(m))+ys*cosd(beta(m))))); % measurments
        % meas = exp(-1i*k*sind(beta(m)).*xs_true(:,n)); % measurments
        meas = meas + normrnd(0,meas_noise/sqrt(2),[NE,1])+ 1i*normrnd(0,meas_noise/sqrt(2),[NE,1]); % Added mesurment noise
        data(1+(m-1)*(1+NE):m*(1+NE),n) = [beta(m);meas];
    end
end

X = data;

end