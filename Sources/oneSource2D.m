tic
%% Pararmeters %%
lambda = 1;
k=2*pi/lambda;
N=3000; %Numerical accuracy factor
range = 10;

%% Field calculations %%
[x,y] = meshgrid(linspace(-(range)*lambda,(range)*lambda, N), linspace(-(range)*lambda,(range)*lambda, N)); %Grid
Z = -1i/4*conj(besselh(0,k*sqrt(x.^2+y.^2)));

ax = -range*lambda:range*lambda;

%% Plots %%
figure();
sgtitle("Point source radiation");

% Array elements locations
subplot(2,2,1);
scatter(0,0,30,"filled","b")
title("Array elements locations");
xlabel('X [wavelenght]') 
ylabel('Y [wavelenght]')
xlim([-1 1])
ylim([-1 1])

% Radiation Pattern
subplot(2,2,2);
val=0.03^2;
contour(x,y,abs(Z).^2,[val,val],'b');
title("Radiation Pattern");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');

% Field
subplot(2,2,3);
image('XData',ax,'YData',ax,'CData',real(Z),'CDataMapping', 'scaled');
title("Field");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range*lambda range*lambda]);
ylim([-range*lambda range*lambda]);

%power
subplot(2,2,4);
image('XData',ax,'YData',ax,'CData',log10(1+(abs(Z).^2)),'CDataMapping', 'scaled');
title("Power (Log Scale)");
xlabel('X [wavelenght]');
ylabel('Y [wavelenght]');
xlim([-range*lambda range*lambda]);
ylim([-range*lambda range*lambda]);

toc