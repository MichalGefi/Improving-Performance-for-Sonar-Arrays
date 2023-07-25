tic

%% Wave Pararmeters %%
lambda = 1;
alpha = 42;
boundery = 2;

%% Array Pararmeters %%
D=0.3*lambda; %the distance between elements
NE = 50; %Number of elements

%% Runtime Pararmeters %%
N = 10^5;

%% Setup %%
k=2*pi/lambda;
theta = linspace(boundery,180-boundery,N); %Estimated posible angles

%% Elements %%
if NE>1 
    xs = -D*((NE-1)/2):D:D*((NE-1)/2); 
    ys = zeros(1, NE);
else
    xs = 0;
    ys = 0;
end

meas = exp(-1i*(k.*(xs*sind(alpha)+ys*cosd(alpha)))); %Horizontal vector of elemnts plane wave field measurments
phi_theta = zeros(N,NE); %Matrix of estimated measurments for each posible angle

for j = 1:length(theta)
    phi_theta(j,:) = exp(-1i*(k.*(xs*sind(theta(j))+ys*cosd(theta(j)))));
end

coef = ones(N,NE);
for j = 1:length(theta)
    coef(j,:) = coef(j,:).*phi_theta(j,:);%Horizontal vector of elemnts field guessed measurements
end

out_array = zeros(N,1);
for j = 1:length(theta)
    out_array(j,:) = abs((conj(coef(j,:))*meas.')./NE);
end

out_array = out_array./trapz(theta,out_array); % normalization

max_out_vals = maxk(out_array, 2);
indexs_of_max_valas = find(out_array>=max_out_vals(2));
out = theta(indexs_of_max_valas(1));
out_angle = out;
err_prc=100*abs((out-alpha)/alpha);

%% Plots %%
% Array elements locations
figure();
hold on;
scatter(xs,ys,10,"filled","b")
hold off;
title("Array elements locations");
formatSpec = "NE=%d, D=%0.1f";
subtitle(sprintf(formatSpec,NE,D));
legend('Nominal')
xlabel('X [wavelenght]') 
ylabel('Y [wavelenght]')
ylim([-1.3*NE*D/2 1.3*NE*D/2])
xlim([-0.1, 0.1])

% angle estimation
figure();
hold on
plot(theta,out_array,"b");
hold off
title_str = sprintf(['Phased Array: True angle=%.2f' char(176)],alpha);
title(title_str);
dataline = "NE=%d, D=%0.1f[Wavelenght]";
outline = "Estimated angle %0.3f"+char(176)+" Error=%0.3f";
subtitle({sprintf(dataline,NE,D),sprintf(outline,out,err_prc)+"%"});
legend('Nominal');
xlabel("angle[deg]");
ylabel('Estimator'); 
toc 
