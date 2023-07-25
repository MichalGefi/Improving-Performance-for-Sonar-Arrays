%% Pararmeters %%
lambda = 1;
k=2*pi/lambda;
W=1*lambda; %the widthes of the slits
D=4*lambda;%the distance between the slits
Amp = 1;
xs = -5*lambda; %xs must be smaller than 0
ys = 5*lambda;
range_x=10;
range_y=10;

N=10^3;
%N must be odd
if(mod(N,2) == 0)
    N = N+1;
end

%% Main %%
tic

[x,y] = meshgrid(linspace(-range_x*lambda,range_x*lambda, N), linspace(-range_y*lambda,range_y*lambda, N));
phi_0 =  -1i/4*Amp*conj(besselh(0,k*sqrt((x-xs).^2+(y-ys).^2)));%point source
step = (x<0).*1;
phi_1 = phi_0.*step;%plane wave half space


window = (x == 0).*(((-W/2<(y-D/2))&((y-D/2)<W/2))|((-W/2<(y+D/2))&((y+D/2)<W/2)));

sources = phi_0.*window.*cos(atan2((y-ys),(x-xs)));

column_index = floor(N/2)+1;

current = zeros(N,N);
total = zeros(N,N);

for m = 1:N
    if (sources(m,column_index) ~= 0)
        current(:,:) = -1i/4*sources(m,column_index)*conj(besselh(0,k*sqrt((x-x(m,column_index)).^2+(y-y(m, column_index)).^2)));
        total = total + current;
    end
end

phi_2 = ~step.*total;

phi = phi_1 + phi_2;
toc

%% Plots %%
ax_x = -range_x*lambda:range_x*lambda;
ax_y = -range_y*lambda:range_y*lambda;
figure();
sgtitle("Spherical wave radiating through two slits");

% Field real value
subplot(1,2,1);
image('XData',ax_x,'YData',ax_y,'CData',real(phi),'CDataMapping', 'scaled');
title("Field");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range_x*lambda range_x*lambda]);
ylim([-range_y*lambda range_y*lambda]);

% Field absolute value
subplot(1,2,2);
image('XData',ax_x,'YData',ax_y,'CData',log10(1+(abs(phi).^2)),'CDataMapping', 'scaled');% +1 used for visual clarity
title("Power (Log Scale)");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range_x*lambda range_x*lambda]);
ylim([-range_y*lambda range_y*lambda]);