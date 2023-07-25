%% README %%
% This script showcases a demonstration of the Phased Array results
% appling improvments using the elements location estimator (AVG method)
% Due to the stochastic nature in order to achive reliable and consistant results
% parametrs must be setup correctly and random variables within limits accoring to them.
% Since real data and SONAR capabilities are classified, all parameters are variables.

tic
%% Wave Pararmeters %%
lambda = 1;
phase = 0; % not tested
alpha = 25;
boundery = 2; %degrees

%% Array Pararmeters %%
axis = 1; % of the array (0->x, 1->y)
D = 1*lambda;% the distance between elements
NE = 30; % Number of elements

%% Array Pararmeters %%
sigma = 0.12*lambda; % The std of the diviation in the locations
curve = 0.19*(lambda^-1); % curve parameter
meas_noise = 0.3; % measurment_noise

%% Runtime Pararmeters %%
N2 = 10^5; % Angle numerical accuracy parameter
N3 = 10; % Number of plane waves sent to the system in known angles

%% Setup %%
k=2*pi/lambda;
theta = linspace(boundery,90-boundery,N2); % Estimated posible angles
phi_theta = zeros(N2,NE); % Matrix of estimated measurments for each posible angle
phi_theta_est = zeros(N2,NE); % Matrix of estimated measurments for each posible angle using estimated locations

%% Element locations setup %%
ys = -D*((NE-1)/2):D:D*((NE-1)/2);
ys_dev = ys;
xs = zeros(1, NE);
xs_dev = zeros(1, NE) + normrnd(0,sigma,[1,NE]) + (ys.^2)*curve./(NE*D);
xs_dev = xs_dev-mean(xs_dev); 

%% element location estimator %%
%Running over different angles of the wave to deal with noise
estimated_xs = AVG_Est(k, boundery, N3, NE, phase, xs_dev, ys_dev, meas_noise);

%% phased array %%

% vector of elemnts plane wave field measurments with and without deveation
meas = exp(-1i*(k.*(xs*sind(alpha)+ys*cosd(alpha))+deg2rad(phase))); %nominal
meas_dev = exp(-1i*(k.*(xs_dev*sind(alpha)+ys_dev*cosd(alpha))+deg2rad(phase))); %with deviations and curve
meas_dev = meas_dev+meas_noise.*(rand(1,length(meas_dev))-0.5); % added noise

% estimated measurments of the plane wave field
for j = 1:length(theta)
    phi_theta(j,:) = exp(-1i*(k.*(xs*sind(theta(j))+ys*cosd(theta(j)))));
    phi_theta_est(j,:) = exp(-1i*(k.*(estimated_xs*sind(theta(j))+ys*cosd(theta(j)))));

end

% applying coefficients to each element (all 1 for phased array)
coef = ones(N2,NE);
for j = 1:length(theta)
    phi_theta(j,:) = coef(j,:).*phi_theta(j,:);%Horizontal vector of elemnts field guessed measurements
    phi_theta_est(j,:) = coef(j,:).*phi_theta_est(j,:);%Horizontal vector of elemnts field estimated measurements
end

% applying estimator
out_array = zeros(N2,1);
out_array_dev = zeros(N2,1);
out_array_est = zeros(N2,1);

for j = 1:length(theta)
    out_array(j,:) = abs((conj(phi_theta(j,:))*meas.')./NE);
    out_array_dev(j,:) = abs((conj(phi_theta(j,:))*meas_dev.')./NE);
    out_array_est(j,:) = abs((conj(phi_theta_est(j,:))*meas_dev.')./NE);
end

% normalization
out_array = out_array./trapz(theta,out_array);
out_array_dev = out_array_dev./trapz(theta,out_array_dev);
out_array_est = out_array_est./trapz(theta,out_array_est);

%% Output %%
% nominal
max_out_vals = [max(out_array(1:round(length(out_array)/2))),max(out_array(round(length(out_array)/2):end))];
indexs_of_max_valas = [find(out_array(1:round(length(out_array)/2))>=max_out_vals(1)), (length(out_array)/2 + find(out_array(round(length(out_array)/2):end)>=max_out_vals(2)))] ;

% with deviations curve and noise
max_out_vals_dev = [max(out_array_dev(1:round(length(out_array_dev)/2))),max(out_array_dev(round(length(out_array_dev)/2):end))];
indexs_of_max_valas_dev = [find(out_array_dev(1:round(length(out_array_dev)/2))>=max_out_vals_dev(1)), (length(out_array_dev)/2 + find(out_array_dev(round(length(out_array_dev)/2):end)>=max_out_vals_dev(2)))] ;


% with deviations curve and noise using estimated locations
max_out_vals_est = [max(out_array_est(1:round(length(out_array_est)/2))),max(out_array_est(round(length(out_array_est)/2):end))];
indexs_of_max_valas_est = [find(out_array_est(1:round(length(out_array_est)/2))>=max_out_vals_est(1)), (length(out_array_est)/2 + find(out_array_est(round(length(out_array_est)/2):end)>=max_out_vals_est(2)))] ;

% selecting the correct angle
if alpha < 90
    out = theta(indexs_of_max_valas(1));
    out_dev = theta(indexs_of_max_valas_dev(1));
    out_est = theta(indexs_of_max_valas_est(1));
else
    out = theta(indexs_of_max_valas(2));
    out_dev = theta(indexs_of_max_valas_dev(2));
    out_est = theta(indexs_of_max_valas_est(2));    
end

% output angle
out_angle = phased_array(k, N2, NE, xs, ys, meas, boundery, coef, alpha);
out_angle_dev = phased_array(k, N2, NE, xs, ys, meas_dev, boundery, coef, alpha);
out_angle_est = phased_array(k, N2, NE, estimated_xs, ys, meas_dev, boundery, coef, alpha);

% error calculations
err_prc=100*abs((out-alpha)/alpha);
err_prc_dev=100*abs((out_dev-alpha)/alpha);
err_prc_est=100*abs((out_est-alpha)/alpha);

T = toc; disp(['Execution took ' sprintf('%.2f', 1000*T) 'ms.']);

%% Plots %%
% Array elements locations
figure();
hold on;
scatter(xs,ys,10,"filled","b")
scatter(xs_dev,ys_dev,10,"filled","r")
scatter(estimated_xs,ys_dev,10,'x',"g")
hold off;
title("Array elements locations");
formatSpec = "NE=%d, D=%0.1f, curve=%0.2d";
subtitle(sprintf(formatSpec,NE,D,curve));
legend('Nominal','With deviation and curve', 'Estimated')
xlabel('X [wavelenght]') 
ylabel('Y [wavelenght]')
lim_x = 2*max(max(abs(xs_dev)),max(abs(estimated_xs)));
xlim([-lim_x lim_x])
ylim([-1.3*NE*D/2 1.3*NE*D/2])


% angle estimation
figure();
hold on
plot(theta, out_array,"b");
plot(theta, out_array_dev, "r");
plot(theta, out_array_est, "g" )
hold off
title_str = sprintf(['Phased Array: True angle=%.2f' char(176)],alpha);
title(title_str);
dataline = "NE=%d, D=%0.1f[Waveleght], Noise=%0.2f;";
outline = "Nominal estimated angle %0.3f"+char(176)+" Error=%0.3f";
outline_dev = "Estimated angle with error %0.3f"+char(176)+" Error=%0.3f";
outline_est = "Estimated angle with algorithm %0.3f"+char(176)+" Error=%0.3f";
subtitle({sprintf(dataline,NE,D,meas_noise),sprintf(outline,out,err_prc)+"%",sprintf(outline_dev,out_dev,err_prc_dev)+"%",sprintf(outline_est,out_est,err_prc_est)+"%"});
legend('Nominal','Deviation curve and error', 'With estimated locations');
xlabel("angle[deg]");
ylabel('Estimator'); 
