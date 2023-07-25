tic

%% Wave Pararmeters %%
lambda = 1;
phase = 0; % not tested
alpha = 64; % degrees
boundery = 2; % degrees

%% Array Pararmeters %%
axis = 0; % Of the array (0->x, 1->y) (unresolved issue for axis=1)
D = 0.3*lambda;% The distance between elements
NE = 55; % Number of elements

%% Project Array Pararmeters %%
dev = 5*(10^-2)*lambda; %maximal devation distance
sigma1 = 0.1;
curve = 5*(10^-2)*(lambda^-1); %curve parameter
meas_noise = 3*(10^-0); %measurment_noise


%% Runtime Pararmeters %%
N2 = 10^4; %Angle numerical accuracy factor

%% Setup %%
k=2*pi/lambda;
xs = 0;
xs_dev = 0;
ys = 0;
ys_dev=0;
theta = linspace(boundery,180-boundery,N2); %Estimated posible angles
phi_theta = zeros(N2,NE); %Matrix of estimated measurments for each posible angle

%% Elements %%
if NE>1
    if axis==0
        xs = -D*((NE-1)/2):D:D*((NE-1)/2);
        xs_dev = xs;
        ys = zeros(1, NE);
        ys_dev = zeros(1, NE)+dev.*(normrnd(0,sigma1,[1,NE]))+(xs.^2)*curve./(NE*D);
        ys_dev = ys_dev-mean(ys_dev);
    else
        ys = -D*((NE-1)/2):D:D*((NE-1)/2);
        ys_dev = ys;
        xs = zeros(1, NE);
        xs_dev = zeros(1, NE)+dev.*(normrnd(0,sigma1,[1,NE]))+(ys.^2)*curve./(NE*D);
        xs_dev = xs_dev-mean(xs_dev);
    end
end

%vector of elemnts plane wave field measurments with and without deveation
meas = exp(-1i*(k.*(xs*sind(alpha)+ys*cosd(alpha))+deg2rad(phase))); %nominal
meas_dev = exp(-1i*(k.*(xs_dev*sind(alpha)+ys_dev*cosd(alpha))+deg2rad(phase))); %with deviations and curve
meas_dev = meas_dev+meas_noise.*(rand(1,length(meas_dev))-0.5); % added noise

%estimated measurments of the plane wave field
for j = 1:length(theta)
    phi_theta(j,:) = exp(-1i*(k.*(xs*sind(theta(j))+ys*cosd(theta(j)))));
end

%applying coefficients to each element (all 1 for phased array)
coef = ones(N2,NE);
for j = 1:length(theta)
    coef(j,:) = coef(j,:).*phi_theta(j,:);%Horizontal vector of elemnts field guessed measurements
end

%applying estimator
out_array = zeros(N2,1);
out_array_dev = zeros(N2,1);
for j = 1:length(theta)
    out_array(j,:) = abs((conj(coef(j,:))*meas.')./NE);
    out_array_dev(j,:) = abs((conj(coef(j,:))*meas_dev.')./NE);
end

% normalization
out_array = out_array./trapz(theta,out_array);
out_array_dev = out_array_dev./trapz(theta,out_array_dev);



%% Output %%
%nominal
max_out_vals = maxk(out_array, 2);
indexs_of_max_valas = find(out_array>=max_out_vals(2));
out = theta(indexs_of_max_valas(1));
err_prc=100*abs((out-alpha)/alpha);

%with deviations curve and noise
max_out_vals_dev = maxk(out_array_dev, 2);
indexs_of_max_valas_dev = find(out_array_dev>=max_out_vals_dev(2));
out_dev = theta(indexs_of_max_valas_dev(1));
err_prc_dev=100*abs((out_dev-alpha)/alpha);


%% Plots %%
% Array elements locations
figure();
hold on;
scatter(xs,ys,10,"filled","b")
scatter(xs_dev,ys_dev,10,"filled","r")
hold off;
title("Array elements locations");
formatSpec = "NE=%d, D=%0.1f, Devation=%0.2d curve=%0.2d";
% subtitle(sprintf(formatSpec,NE,D,dev,curve));
legend('Nominal','With deviation and curve')
xlabel('X [wavelenght]') 
ylabel('Y [wavelenght]')

% axis=abs(1-axis); %used for plot visuals
if axis==0
    xlim([-1.3*NE*D/2 1.3*NE*D/2])
    ylim([-3*(dev+0.2) 3*(dev+0.2)])
else
    xlim([-3*(dev+0.2) 3*(dev+0.2)])
    ylim([-1.3*NE*D/2 1.3*NE*D/2])
end

% angle estimation
figure();
hold on
plot(theta,out_array,"b");
plot(theta,out_array_dev, "r");
hold off
title_str = sprintf(['Phased Array: True angle=%.2f' char(176)],alpha);
title(title_str);
dataline = "NE=%d, D=%0.1f[Waveleght], Noise=%0.2f;";
outline = "Nominal estimated angle %0.3f"+char(176)+" Error=%0.3f";
outline_dev = "Actual estimated angle %0.3f"+char(176)+" Error=%0.3f";
subtitle({sprintf(dataline,NE,D,meas_noise),sprintf(outline,out,err_prc)+"%",sprintf(outline_dev,out_dev,err_prc_dev)+"%"});
legend('Nominal','Deviation curve and error');
xlabel("angle[deg]");
ylabel('Estimator'); 

toc 