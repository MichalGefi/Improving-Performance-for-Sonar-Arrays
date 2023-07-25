tic

%% Wave Pararmeters %%
lambda = 1;
phase = 0; %not tested
alpha = 30; %degrees
boundery = 2; %degrees

%% Array Pararmeters %%
D=0.5*lambda;%the distance between elements
NE = 45; %Number of elements
dev = 0.3*lambda; %maximal devation distance
axis = 0; % of the array (0->x, 1->y) (unresolved issue for axis=1)

%% Runtime Pararmeters %%
N2 = 10^5; %Angle numerical accuracy factor

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
        ys_dev = zeros(1, NE)+2*dev.*(rand(1,NE)-0.5);
    else
        ys = -D*((NE-1)/2):D:D*((NE-1)/2);
        ys_dev = ys;
        xs = zeros(1, NE);
        xs_dev = zeros(1, NE)+2*dev.*(rand(1,NE)-0.5);
    end
end

%vector of elemnts plane wave field measurments with and without deveation
meas = exp(-1i*(k.*(xs*sind(alpha)+ys*cosd(alpha))+deg2rad(phase)));
meas_dev = exp(-1i*(k.*(xs_dev*sind(alpha)+ys_dev*cosd(alpha))+deg2rad(phase)));

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
%without deviation
max_out_vals = maxk(out_array, 2);
indexs_of_max_valas = find(out_array>=max_out_vals(2));
out = theta(indexs_of_max_valas(1));
err_prc=100*abs((out-alpha)/alpha);
%with deviation
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
formatSpec = "NE=%d, D=%0.1f, Devation=%0.2f[Waveleght]";
subtitle(sprintf(formatSpec,NE,D,dev));
legend('Nominal','With deviation')
xlabel('X [wavelenght]') 
ylabel('Y [wavelenght]')
% axis=abs(1-axis); %used for plot visuals
if axis==0
    xlim([-1.3*NE*D/2 1.3*NE*D/2])
    ylim([-5*dev 5*dev])
else
    xlim([-5*dev 5*dev])
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
dataline = "NE=%d, D=%0.1f[Waveleght], Devation=%0.2f[Waveleght];";
outline = "Estimated angle without deviation %0.3f"+char(176)+" Error=%0.3f";
outline_dev = "Estimated angle with deviation %0.3f"+char(176)+" Error=%0.3f";
subtitle({sprintf(dataline,NE,D,dev),sprintf(outline,out,err_prc)+"%",sprintf(outline_dev,out_dev,err_prc_dev)+"%"});
legend('Nominal','With deviation');
xlabel("angle[deg]");
ylabel('Y'); 

toc 