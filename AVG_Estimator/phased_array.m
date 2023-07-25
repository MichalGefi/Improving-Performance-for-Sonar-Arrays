function [out_angle] = phased_array(k, N2, NE, xs, ys, meas, boundery, coef, alpha)
    theta = linspace(boundery,180-boundery,N2); %Estimated posible angles
    % estimated measurments of the plane wave field
    phi_theta = zeros(N2,NE); %Matrix of estimated measurments for each posible angle
    for j = 1:length(theta)
        phi_theta(j,:) = exp(-1i*(k.*(xs*sind(theta(j))+ys*cosd(theta(j)))));
    end

    % applying coefficients to each element (all 1 for phased array)
    for j = 1:length(theta)
        phi_theta(j,:) = coef(j,:).*phi_theta(j,:);%Horizontal vector of elemnts field estimated measurements
    end

    % applying estimator
    out_array = zeros(N2,1);

    for j = 1:length(theta)
        out_array(j,:) = abs((conj(phi_theta(j,:))*meas.')./NE);
    end
    
    out_array = out_array./trapz(theta,out_array); % normalization

    %% Output %%
    max_out_vals = [max(out_array(1:round(length(out_array)/2))),max(out_array(round(length(out_array)/2):end))];
    indexs_of_max_valas = [find(out_array(1:round(length(out_array)/2))>=max_out_vals(1)), (length(out_array)/2 + find(out_array(round(length(out_array)/2):end)>=max_out_vals(2)))] ;
 
    %selecting the correct angle
    if alpha < 90
        out_angle = theta(indexs_of_max_valas(1));
    else
        out_angle = theta(indexs_of_max_valas(2));  
    end
end

