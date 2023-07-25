tic
%% README %%
% This script showcases a demonstration of the Arrays Radiation pattern
% appling the elements location estimator for assesment (AVG method)
% Due to the stochastic nature in order to achive reliable and consistant results
% parametrs must be setup correctly and random variables within limits accoring to them.
% Since real data and SONAR capabilities are classified, all parameters are variables.

%% Wave Pararmeters %%
lambda = 1;
phase = 0; % not tested
boundery = 2; %degrees
N3 = 20; %Number of plane waves sent to the system in known angles

%% Array Pararmeters %%
D = 1*lambda;%the distance between elements
NE = 50; %Number of elements

%% Array Pararmeters %%
sigma = 0.12*lambda; % The std of the diviation in the locations
curve = 0.19*(lambda^-1); % curve parameter
meas_noise = 0.3; % measurment_noise

%% View distace %%
range_x = 250;
range_y = 2000;

%% Runtime Pararmeters %%
N=3; %Numerical accuracy factor
N_x=range_x*N*lambda;
N_y=range_y*N*lambda;

%% Setup %%
k=2*pi/lambda;
xs = 0;
xs_dev = 0;
ys = 0;
ys_dev=0;
Amp = ones(1,NE);
z = zeros(NE,N_y,N_x);
z_dev = zeros(NE,N_y,N_x);
z_est = zeros(NE,N_y,N_x);

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
    meas = exp(-1i*(k.*(xs*sind(alpha(j))+ys*cosd(alpha(j)))+deg2rad(phase))); %nominal
    meas_dev = exp(-1i*(k.*(xs_dev*sind(alpha(j))+ys_dev*cosd(alpha(j)))+deg2rad(phase))); %with deviations and curve
    meas_dev = meas_dev+meas_noise.*(rand(1,length(meas_dev))-0.5); % added noise

% multipling the measured field by factors that cancle the depedence on
% ys_dev and the pahse of the known wave. We assume that there is no error
% in ys and that we know it percicly
    mult_factor = exp(1i*(k.*ys_dev*cosd(alpha(j))+deg2rad(phase)));
    temp = meas_dev.*mult_factor;

% extracting xs
    temp = angle(temp);
    estimated_xs = estimated_xs - temp./(k.*sind(alpha(j)));
end
estimated_xs = estimated_xs/N3;

%% Field calculations %%
wb1 = waitbar(0,'Loading...');
[x,y] = meshgrid(linspace(-(range_x)*lambda,(range_x)*lambda, N_x), linspace(-(range_y)*lambda,(range_y)*lambda, N_y)); %Grid
for m = 1:NE
    z(m,:,:) = -1i/4*Amp(m)*conj(besselh(0,k*sqrt((x-xs(m)).^2+(y-ys(m)).^2)));
    z_dev(m,:,:) = -1i/4*Amp(m)*conj(besselh(0,k*sqrt((x-xs_dev(m)).^2+(y-ys_dev(m)).^2)));
    z_est(m,:,:) = -1i/4*Amp(m)*conj(besselh(0,k*sqrt((x-estimated_xs(m)).^2+(y-ys_dev(m)).^2)));
    waitbar(m/NE,wb1)
end

z = squeeze(sum(z));
z_dev = squeeze(sum(z_dev));
z_est = squeeze(sum(z_est));

T = toc; disp(['Execution took ' sprintf('%.2f', T) ' seconds.']);
delete(wb1)

%% Plots %%
% Array elements locations
figure();
hold on;
scatter(xs,ys,10,"filled","b")
scatter(xs_dev,ys_dev,10,"filled","r")
scatter(estimated_xs,ys_dev,10,'x',"g")
hold off;
title("Array elements locations");
formatSpec = "NE=%d, D=%0.1f, Devation=%0.2d curve=%0.2d";
subtitle(sprintf(formatSpec,NE,D,dev,curve));
legend('Nominal','With deviation and curve', 'Estimated')
xlabel('X [wavelenght]') 
ylabel('Y [wavelenght]')
xlim([-3*(dev+0.2) 3*(dev+0.2)])
ylim([-1.3*NE*D/2 1.3*NE*D/2])

% Radiation Pattern
figure();
hold on;
val=0.1^2;
contour(y,x,abs(z).^2,[val,val],'b');
contour(y,x,abs(z_dev).^2,[val,val],'r');
contour(y,x,abs(z_est).^2,[val,val],'g');
hold off;
title("Radiation Pattern");
formatSpec = "NE=%d, D=%0.1f, Devation=%0.2f;";
subtitle(sprintf(formatSpec,NE,D,dev));
legend('Nominal','With deviation', 'Estimated');
xlabel('X [wavelenght]'); 
ylabel('Y [wavelenght]');
