function [estimated_xs] = AVG_Est(k, boundery, N3, NE, phase, xs_dev, ys_dev, meas_noise)
    %% element location estimator %%
    % Running over different angles of the wave to deal with noise
    beta = linspace(boundery,180-boundery,N3);
    estimated_xs = zeros(1, NE);
    for j = 1:length(beta)
    % vector of elemnts plane wave field measurments with and without deveation
        meas_dev = exp(-1i*(k.*(xs_dev*sind(beta(j))+ys_dev*cosd(beta(j)))+deg2rad(phase))); %with deviations and curve
        meas_dev = meas_dev+meas_noise.*(rand(1,length(meas_dev))-0.5); % added noise
    
    % multipling the measured field by factors that cancle the depedence on
    % ys_dev and the pahse of the known wave. We assume that there is no error
    % in ys and that we know it percicly
        mult_factor = exp(1i*(k.*ys_dev*cosd(beta(j))+deg2rad(phase)));
        temp = meas_dev.*mult_factor;
    
    % extracting xs
        temp = angle(temp);
        estimated_xs = estimated_xs - temp./(k.*sind(beta(j)));
    end
    estimated_xs = estimated_xs/N3;
end

