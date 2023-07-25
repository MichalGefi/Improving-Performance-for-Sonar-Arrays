%% Pararmeters %%
lambda = 1;
k=2*pi/lambda;
W=1*lambda; %the width of the slit
Amp = 2;
xs = -2*lambda; %xs must be smaller than 0
ys = 7*lambda;

N=10^3;
%N must be odd
if(mod(N,2) == 0)
    N = N+1;
end

%% Main %%
tic

[x,y] = meshgrid(linspace(-10*lambda,10*lambda, N), linspace(-10*lambda,10*lambda, N));
phi_0 = -1i/4*Amp*conj(besselh(0,k*sqrt((x-xs).^2+(y-ys).^2)));%point source
step = (x<0).*1;
phi_1 = phi_0.*step;%plane wave half space

window = (x == 0).*(abs(y)<W/2);

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
sgtitle("Spherical wave radiating through a slit");

% Field real value
subplot(1,2,1);
image('XData',ax_x,'YData',ax_y,'CData',2*real(phi),'CDataMapping', 'scaled');
title("Field");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range_x*lambda range_x*lambda]);
ylim([-range_y*lambda range_y*lambda]);

% Field absolute value
subplot(1,2,2);
image('XData',ax_x,'YData',ax_y,'CData',log10(1+abs(phi).^2),'CDataMapping', 'scaled');% +1 used for visual clarity
title("Power (Log Scale)");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range_x*lambda range_x*lambda]);
ylim([-range_y*lambda range_y*lambda]);