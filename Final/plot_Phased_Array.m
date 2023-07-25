function [] = plot_Phased_Array(N,alpha,ideal_out,ideal_angle,NOM_out,NOM_angle,AVG_out,AVG_angle,LR_out,LR_angle)
boundery = 1;
theta = linspace(boundery,90,N/2); % posible estimated angles
% ploting estimator
figure();
hold on;
    plot(theta,ideal_out(1:(N/2)),"r");
    plot(theta,NOM_out(1:(N/2)),"b");
    plot(theta,AVG_out(1:(N/2)),"g");
    plot(theta,LR_out(1:(N/2)),"magenta");
hold off;
legend('Ideal','Nominal','AVG estimated','LR estimated',"fontsize", 18);
title_str = sprintf('Phased Array  Angle=%g°' ,alpha);
title(title_str,"fontsize", 28);
l1 = "Estimated angles: Ideal=%.3f°, Nominal=%.3f°, Avg=%.3f°, LR=%.3f°";
subtitle({sprintf(l1,ideal_angle,NOM_angle,AVG_angle,LR_angle)},"fontsize", 20);
xlabel("Angle [deg]","fontsize", 20);
ylabel('Estimator',"fontsize", 20);
ylim([0 1.05*max(ideal_out)])
end

