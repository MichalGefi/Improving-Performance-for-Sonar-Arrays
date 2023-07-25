tic
%% Pararmeters %%
lambda = 1;
k=2*pi/lambda;
N=1000*lambda; %The number of rows and columns in the matrix
d=3;
x1 = -d*lambda;
y1 = 0;
x2 = d*lambda;
y2 = 0;

range_x=10;
range_y=10;

xx=[x1,x2];
yy=[y1,y2];

%% Field calculations %%
[x,y] = meshgrid(linspace(-range_x*lambda,range_x*lambda, N), linspace(-range_y*lambda,range_y*lambda, N));
z1 = -1i/4*conj(besselh(0,k*sqrt((x-x1).^2+(y-y1).^2)));
z2 = -1i/4*conj(besselh(0,k*sqrt((x-x2).^2+(y-y2).^2)));
Z = z1+z2;

%% Plots %%
ax_x = -range_x*lambda:range_x*lambda;
ax_y = -range_y*lambda:range_y*lambda;

figure();
sgtitle("Two point sources radiation");

% Field
subplot(1,2,1);
image('XData',ax_x,'YData',ax_y,'CData',real(Z),'CDataMapping', 'scaled');
title("Field");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range_x*lambda range_x*lambda]);
ylim([-range_y*lambda range_y*lambda]);

%power
subplot(1,2,2);
image('XData',ax_x,'YData',ax_y,'CData',log10(1+(abs(Z).^2)),'CDataMapping', 'scaled');
title("Power (Log Scale)");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range_x*lambda range_x*lambda]);
ylim([-range_y*lambda range_y*lambda]);

toc