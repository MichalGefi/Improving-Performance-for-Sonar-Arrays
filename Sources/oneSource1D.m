lambda = 1;
k=2*pi/lambda;
N=1000; %The number of rows and columns in the matrix
range_x=10;

x = linspace(-range_x*lambda,range_x*lambda, N);
y = -1i/4*conj(besselh(0,k*abs(x)));

%% Plots %%
figure();
sgtitle("Point source radiation in 1D");

% Field real value
subplot(1,2,1);
plot(x,real(y));
title("Field");
xlabel('X [wavelenght]');
ylabel('phi');
xlim([-range_x*lambda range_x*lambda]);

% Field absolute value
subplot(1,2,2);
plot(x,log10(1+(abs(y).^2)));
title("Power (Log Scale)");
xlabel('X [wavelenght]');
ylabel('power');
xlim([-range_x*lambda range_x*lambda]);