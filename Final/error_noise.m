%% README
% This script runs the project using the "err_noise_run" function.
% Using this code one can review improvment in SONAR accuracy
% as the mesurments noise increases.
% Due to the stochastic nature, in order to showcase a trend a median accuracy is used
% A median is used since outlier are relatively common.

%% Parameters
N = 500; % number of curves mesured
K = 20; % number of mesurments averged over
section = 50; % the amount of time the LR will be trained 
boundery = 2;
alpha = 40;
param = linspace(10^-45, 1,N);

%% Main section
wb1 = waitbar(0,'Executing...');
tic
err_mat = zeros([3,K,N]);
for j = 1:N
    if (1 == mod(j,round(N/section))) && j~=N
        [~,~,~,~] = err_curve_run(param(j+N/section-1),true);
    end
    for k = 1:K
    [~,NOM,AVG,LR] = err_curve_run(param(j),false);
    err_mat(1,k,j)= abs(alpha-NOM)./alpha;
    err_mat(2,k,j)= abs(alpha-AVG)./alpha;
    err_mat(3,k,j)= abs(alpha-LR)./alpha;
    waitbar((k/K+j-1)/N,wb1);
    end
end

err = zeros([3,N]);
for j = 1:N
    err(1,j) = median(err_mat(1,:,j));
    err(2,j) = median(err_mat(2,:,j));
    err(3,j) = median(err_mat(3,:,j));
end

delete(wb1)
T = toc; disp(['Execution took ' sprintf('%.2f', T/60) ' minutes.']);

%% Plots
l=150; % moving avarge (used to showcase general trend better)
figure();
hold on
    plot(param,smooth(err(1,:),l),"b");
    plot(param,smooth(err(2,:),l),"g");
    plot(param,smooth(err(3,:),l),"magenta");
hold off
title("Phased Array Error","fontsize", 28);
%formatSpec = "NE=%d, D=%0.1f, r=%d";
%subtitle(sprintf(formatSpec,NE,D,r));
legend('Nominal','AVG estimated','LR estimated',"fontsize", 20, 'Location', 'NorthWest');
xlabel('Noise level',"fontsize", 20); 
ylabel('Error (Log-Scale)',"fontsize", 20);
set(gca, 'YScale', 'log');

