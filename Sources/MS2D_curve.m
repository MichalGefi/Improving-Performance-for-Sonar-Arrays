tic

%% Array Pararmeters %%
axis = 1; %of the array (0->x, 1->y)
lambda = 1; %wavelength
D = 0.3*lambda; %Distance between elements
NE = 30; %Number of elements

dev = 0*lambda; %maximal devation distance
sigma1 = 0 ;
curve = 5*lambda*(10^-3); %curve parameter

%% View distace %%
if axis==0
    range_x = 10;
    range_y = 260;
else
    range_x = 260;
    range_y = 10;
end

%% Runtime Pararmeters %%
N=30; %Numerical accuracy factor
N_x=range_x*N*lambda;
N_y=range_y*N*lambda;

%% Setup %%
k=2*pi/lambda;
Amp = ones(1,NE);
z = zeros(NE,N_y,N_x);
z_dev = zeros(NE,N_y,N_x);
xs = 0;
xs_dev = 0;
ys = 0;
ys_dev=0;

%% Elements %%
if NE>1
    if axis==0
        xs = -D*((NE-1)/2):D:D*((NE-1)/2);
        xs_dev = xs;
        ys = zeros(1, NE);
        ys_dev = zeros(1, NE)+dev.*(normrnd(0,sigma1,[1,NE]))+(xs.^2)*curve;
        ys_dev = ys_dev-mean(ys_dev);
    else
        ys = -D*((NE-1)/2):D:D*((NE-1)/2);
        ys_dev = ys;
        xs = zeros(1, NE);
        xs_dev = zeros(1, NE)+dev.*(normrnd(0,sigma1,[1,NE]))+(ys.^2)*curve;
        xs_dev = xs_dev-mean(xs_dev);
    end
end


%% Field calculations %%
[x,y] = meshgrid(linspace(-(range_x)*lambda,(range_x)*lambda, N_x), linspace(-(range_y)*lambda,(range_y)*lambda, N_y)); %Grid

for m = 1:NE
    z(m,:,:) = -1i/4*Amp(m)*conj(besselh(0,k*sqrt((x-xs(m)).^2+(y-ys(m)).^2)));
    z_dev(m,:,:) = -1i/4*Amp(m)*conj(besselh(0,k*sqrt((x-xs_dev(m)).^2+(y-ys_dev(m)).^2)));
end
Z = squeeze(sum(z));
Z_dev = squeeze(sum(z_dev));

%% Plots %%


% Array elements locations
figure();
hold on;
scatter(xs,ys,10,"filled","b")
scatter(xs_dev,ys_dev,10,"filled","r")
hold off;
title("Array elements locations");
formatSpec = "NE=%d, D=%0.1f, Devation=%0.2d curve=%0.2d";
subtitle(sprintf(formatSpec,NE,D,dev,curve));
legend('Nominal','With deviation and curve')
xlabel('X [wavelenght]') 
ylabel('Y [wavelenght]')
if axis==0
    xlim([-1.3*NE*D/2 1.3*NE*D/2])
    ylim([-3*(dev+0.2) 3*(dev+0.2)])
else
    xlim([-3*(dev+0.2) 3*(dev+0.2)])
    ylim([-1.3*NE*D/2 1.3*NE*D/2])
end

% Radiation Pattern
val=0.15^2;
figure();
hold on;
contour(x,y,abs(Z).^2,[val,val],'b');
contour(x,y,abs(Z_dev).^2,[val,val],'r');
hold off;
title("Radiation Pattern");
formatSpec = "NE=%d, D=%0.1f, Devation=%0.2d curve=%0.2d";
subtitle(sprintf(formatSpec,NE,D,dev, curve));
legend('Nominal','With deviation and curve');
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');

% Fields
figure();
sgtitle("Array radiation");
subplot(2,2,1);
image(real(Z),'CDataMapping', 'scaled');
title("Field");
subplot(2,2,2);
image(log10(1+(abs(Z).^2)),'CDataMapping', 'scaled');
title("Power (Log scale)");
subplot(2,2,3);
image(real(Z_dev),'CDataMapping', 'scaled');
title("Field with curve");
subplot(2,2,4);
image(log10(1+(abs(Z_dev).^2)),'CDataMapping', 'scaled');
title("Power (Log scale) with curve");

toc 


