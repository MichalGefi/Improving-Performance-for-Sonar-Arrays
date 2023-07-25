tic
%% Array Pararmeters %%
axis = 1; %of the array (0->x, 1->y)
lambda = 1; %wavelength
D = 0.3*lambda; %Distance between elements
NE = 5; %Number of elements

%% View distace %%
% edit this section do display the required range

if axis==0
    range_x = 5;
    range_y = 26;
else
    range_x = 26;
    range_y = 5;
end

%% Runtime Pararmeters %%
N = 50; %Numerical accuracy factor
N_x=range_x*N*lambda;
N_y=range_y*N*lambda;


%% Setup %%
k=2*pi/lambda;
Amp = ones(1,NE);
z = zeros(NE,N_y,N_x);
xs = 0;
ys = 0;
ax_x = -range_x*lambda:range_x*lambda;
ax_y = -range_y*lambda:range_y*lambda;

%% Elements %%
if NE>1
    if axis==0
        xs = -D*((NE-1)/2):D:D*((NE-1)/2);
        ys = zeros(1, NE);
    else
        ys = -D*((NE-1)/2):D:D*((NE-1)/2);
        xs = zeros(1, NE);
    end
end

%%% Field calculations %%%
[x,y] = meshgrid(linspace(-(range_x)*lambda,(range_x)*lambda, N_x), linspace(-(range_y)*lambda,(range_y)*lambda, N_y));

for m = 1:NE
    z(m,:,:) = -1i/4*Amp(m)*conj(besselh(0,k*sqrt((x-xs(m)).^2+(y-ys(m)).^2)));
end
Z = squeeze(sum(z));

%% Plots %%
f=figure();
sgtitle("Array radiation");

% Array elements locations
subplot(2,2,1);
scatter(xs,ys,10,"filled","b")
title("Array elements locations");
xlabel('X [wavelenght]') 
ylabel('Y [wavelenght]')
if axis==0
    xlim([-1.3*NE*D/2 1.3*NE*D/2])
    ylim([-0.5 0.5])
else
    xlim([-0.5 0.5])
    ylim([-1.3*NE*D/2 1.3*NE*D/2])
end

% Radiation Pattern
val= 6.5*10^-3;
subplot(2,2,2);
contour(x,y,abs(Z).^2,[val,val],'b');
title("Radiation Pattern Shape");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');

% Field
subplot(2,2,3);
image('XData',ax_x,'YData',ax_y,'CData',real(Z),'CDataMapping', 'scaled');
title("Field");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range_x*lambda/2 range_x*lambda/2]);
ylim([-range_y*lambda/2 range_y*lambda/2]);
subplot(2,2,4);
image('XData',ax_x,'YData',ax_y,'CData',log10(1+(abs(Z).^2)),'CDataMapping', 'scaled');
title("Power (Log Scale)");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range_x*lambda/2 range_x*lambda/2]);
ylim([-range_y*lambda/2 range_y*lambda/2]);

%% ADDITIONAL PLOTS %%
% xindex = round(0.8*N);
% Zx = Z(xindex,:);
% figure();
% subplot(2,2,1);
% plot(y,real(Zx))
% title("One line const X Real")
% subplot(2,2,2);
% plot(y,abs(Zx))
% title("One line const X Abs")
% 
% 
% % yval= 0.5*lambda;
% % val = yval/(10+M/2)*lambda;
% % yindex = floor((val+1)/2*N);
% yindex = round(0.7*N);
% Zy = Z(:,yindex);
% Zy = transpose(Zy);
% subplot(2,2,3);
% plot(x(1,:), real(Zy))
% title("One line const Y Real")
% subplot(2,2,4);
% plot(x(1,:),abs(Zy))
% title("One line const Y Abs")

toc 


