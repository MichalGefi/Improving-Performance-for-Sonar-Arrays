function [ideal_angle,NOM_angle,AVG_angle,LR_angle] = err_curve_run(param,train)
%% parameters %%
%runtime parameters
NT = 10^5; % number of arrays used for training
NB = 5; % number of angles
N = 10^3; % phased array estimator

%Array parameters
NE = 40; % number of elements in the array
lambda = 1; % wavelength
D = 1; % distance in y axes between elements
boundery = 2; % in degrees
meas_noise = 0.1; % stdv of mesurment device noise
dev_sigma = 0.1; % stdv of deviations for training
dev_max = 0.25; % maximal deviation
curve_sigma = param(end)*0.7; % stdv of curve for training
curve_max = param(end)*1.5; % maximal curve


%% Untested input %%
%parameters for untesteded array
dev_sigma_test = 0.05; % stdv of deviations for test
curve_test = param; % curve for test

%% Setup %%
% Elements
ys = (-D*((NE-1)/2):D:D*((NE-1)/2))';
xs = zeros(NE, 1);
xs_true = zeros(NE, 1) + truncnormrnd([NE,1],0, dev_sigma_test, -dev_max, dev_max)+(ys.^2)*curve_test./(NE*D);
xs_true = xs_true-mean(xs_true);

% measurments
k=2*pi/lambda; % of wave

beta = unifrnd(boundery, 90,[1,NB]);
% beta = linspace(boundery,90-boundery,NB);
data = zeros(NB*(NE+1),1);
for j = 1:NB %running on every angle
    meas = exp(-1i*(k.*(xs_true*sind(beta(j))+ys*cosd(beta(j)))));
    % meas = exp(-1i*k*sind(beta(j)).*xs_true);
    meas = meas + normrnd(0,meas_noise/sqrt(2),[NE,1])+ 1i*normrnd(0,meas_noise/sqrt(2),[NE,1]); % Added mesurment noise
    data(1+(j-1)*(1+NE):j*(1+NE)) = [beta(j);meas];
end

%% Linear Regression %%
persistent W

if train
    %Creating data for LR training
    beta2 = unifrnd(boundery, 90,[1,NB]);
    % beta2 = linspace(boundery,90-boundery,NB);
    [T,X] = create_data(NE,lambda,D,NT, beta2 ,dev_sigma,dev_max,curve_sigma,curve_max,meas_noise);
    
    %Traning LR
    W = LR_train(T,X,NB,lambda,ys);
end

%% Applying location estimatmators %%
xs_LR_est = LR_place_estimator(data,W,NB,lambda,ys);
xs_AVG_est = AVG_place_estimator(data,NB,lambda,ys);

%% Applaying Phased array angle estimator %%
Phased_array = true; % select true is you want to estimate angle with Phased Array
if Phased_array
    % data for Phased Array
    % alpha = unifrnd(boundery, 90);
    alpha = 40;
    meas = exp(-1i*(k.*(xs_true*sind(alpha)+ys*cosd(alpha))));
    meas = meas + normrnd(0,meas_noise/sqrt(2),[NE,1])+ 1i*normrnd(0,meas_noise/sqrt(2),[NE,1]); % Added mesurment noise
    
    % Applaying Phased array
    [~,ideal_angle] = phased_array(lambda,xs_true,D,meas,N); % using ideal locations
    [~,NOM_angle] = phased_array(lambda,xs,D,meas,N); % using nominal locations
    [~,AVG_angle] = phased_array(lambda,xs_AVG_est,D,meas,N); % using AVG locations
    [~,LR_angle] = phased_array(lambda,xs_LR_est,D,meas,N); % using LR locations

end

end

