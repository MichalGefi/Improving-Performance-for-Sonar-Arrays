%% README
% This script runs the project using the "err_dev_run" function.
% Using this code one can review improvment in SONAR accuracy
% as the stdv of the elements location deviation increases.
% Due to the stochastic nature, in order to showcase a trend a median accuracy is used
% A median is used since outlier are relatively common.

%% Parameters
N = 400; % number of stdv's mesured
K = 100; % number of mesurments averged over
boundery = 1;
alpha = 30;
param = linspace(10^-6, 1,N);

%% setup- used to train the LR weights
% A matrix "W" will be generated
wb0 = waitbar(0,'Initializing setup...');
tic
[~,~,~,~] = err_dev_run(param(1),true);
T = toc; disp(['Training LR took ' sprintf('%.2f', T) ' seconds.']);
delete(wb0)

%% Main section
wb1 = waitbar(0,'Executing...');
tic
err_mat = zeros([3,K,N]);
for k = 1:K
    for j = 1:N
    [~,NOM,AVG,LR] = err_dev_run(param(j),false);
    err_mat(1,k,j)= abs(alpha-NOM)./alpha;
    err_mat(2,k,j)= abs(alpha-AVG)./alpha;
    err_mat(3,k,j)= abs(alpha-LR)./alpha;
    waitbar((j/N+k-1)/K,wb1);
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
l=25; % moving avarge (used to showcase general trend better)
figure();
hold on
    plot(param,smooth(err(1,:),l),"b");
    plot(param,smooth(err(2,:),l),"g");
    plot(param,smooth(err(3,:),l),"magenta");
hold off
title("Phased Array Error","fontsize", 28);
%formatSpec = "NE=%d, D=%0.1f, r=%d";
%subtitle(sprintf(formatSpec,NE,D,r));
legend('Nominal','AVG estimated','LR estimated',"fontsize", 20, 'Location', 'SouthEast');
xlabel('σ of deviation [wavelength]',"fontsize", 20); 
ylabel('Error (Log-Scale)',"fontsize", 20);
ylim([0, 2*max(err(1,:))])
set(gca, 'YScale', 'log');
