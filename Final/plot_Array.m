function [] = plot_Array(ys,plot_x_lim,xs,xs_true,xs_AVG_est,xs_LR_est)
% Array elements locations
figure();
hold on
    scatter(ys,xs,50,"X","b")
    scatter(ys,xs_AVG_est,100,"X","g");
    scatter(ys,xs_LR_est,100,"X","magenta");
    scatter(ys,xs_true,15,"filled","r")
hold off
title("Array elements locations","fontsize", 28);
%formatSpec = "NE=%d, D=%0.1f, r=%d";
%subtitle(sprintf(formatSpec,NE,D,r));
legend('Nominal','AVG estimated','LR estimated','True',"fontsize", 20);
xlabel('X [wavelength]',"fontsize", 20); 
ylabel('Y [wavelength]',"fontsize", 20);
ylim([-plot_x_lim plot_x_lim])
end

