%% README %%
% this script showcases a demonstration of the elements location estimator using an average over several wave mesurments.
% Due to the stochastic nature in order to achive reliable and consistant results
% parametrs must be setup correctly and random variables within limits accoring to them.
% Since real data and SONAR capabilities are classified, all parameters are variables.

tic
%% Wave Pararmeters %%
lambda = 1;
boundery = 2; % degrees
N3 = 20; % Number of plane waves sent to the system in known angles

%% Array Pararmeters %%
D = 0.5*lambda;% the distance between elements
NE = 50; % Number of elements

%% Project Array Pararmeters %%
dev = 1*lambda; % maximal devation distance
sigma1 = 0.07;
curve = 0.1; % 5*(10^-2)*(lambda^-1); %curve parameter
meas_noise = 0.3; %measurment_noise


%% Runtime Pararmeters %%
N2 = 10^4; % Angle numerical accuracy factor

%% Setup %%
k=2*pi/lambda;
xs = 0;
xs_dev = 0;
ys = 0;
ys_dev=0;
theta = linspace(boundery,180-boundery,N2); % Estimated posible angles
phi_theta = zeros(N2,NE); % Matrix of estimated measurments for each posible angle

%% Elements %%
if NE>1
        ys = -D*((NE-1)/2):D:D*((NE-1)/2);
        ys_dev = ys;
        xs = zeros(1, NE);
        xs_dev = zeros(1, NE)+dev.*(normrnd(0,sigma1,[1,NE]))+(ys.^2)*curve./(NE*D);
        xs_dev = xs_dev-mean(xs_dev);
end

% Running over different angles of the wave to deal with noise
alpha = linspace(boundery,180-boundery,N3);
estimated_xs = zeros(1, NE);
for j = 1:length(alpha)
% vector of elemnts plane wave field measurments with and without deveation
    meas = exp(-1i*(k.*(xs*sind(alpha(j))+ys*cosd(alpha(j))))); %nominal
    meas_dev = exp(-1i*(k.*(xs_dev*sind(alpha(j))+ys_dev*cosd(alpha(j))))); %with deviations and curve
    meas_dev = meas_dev+meas_noise.*(rand(1,length(meas_dev))-0.5); % added noise

% multipling the measured field by factors that cancle the depedence on
% ys_dev and the pahse of the known wave. We assume that there is no error
% in ys and that we know it percicly
    mult_factor = exp(1i*(k.*ys_dev*cosd(alpha(j))));
    temp = meas_dev.*mult_factor;

% extracting xs
    temp = angle(temp);
    estimated_xs = estimated_xs - temp./(k.*sind(alpha(j)));
end
estimated_xs = estimated_xs/N3;

T = toc; disp(['Execution took ' sprintf('%.2f', 1000*T) 'ms.']);

%% Plots %%
% Array elements locations
figure();
hold on;
scatter(ys,xs,10,"filled","b")
scatter(ys_dev,xs_dev,10,"filled","r")
scatter(ys_dev,estimated_xs,10,'x',"g")
% scatter(xs,ys,10,"filled","b")
% scatter(xs_dev,ys_dev,10,"filled","r")
% scatter(estimated_xs,ys_dev,10,'x',"g")
hold off;
title("Array elements locations");
formatSpec = "NE=%d, D=%0.1f, Devation=%0.2f curve=%0.2f";
subtitle(sprintf(formatSpec,NE,D,dev,curve));
legend('Nominal','With deviation and curve', 'Estimated')
xlabel('X [wavelenght]') 
ylabel('Y [wavelenght]')
lim_x = 1.3*max(max(abs(xs_dev)),max(abs(estimated_xs)));
ylim([-lim_x lim_x])
xlim([-1.3*NE*D/2 1.3*NE*D/2])


