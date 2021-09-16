function [ratio] = ratio_meas_2(images_pulse_ifft_avg_tot, pic, N_mean, plotflag)

if size(pic,2) == 1
    temp = pic
    pic = temp*ones(1,size(images_pulse_ifft_avg_tot,1)-1);
end
ratio = zeros(size(images_pulse_ifft_avg_tot,1) - 1,1);

x = [1:1:size(images_pulse_ifft_avg_tot,2)];
for v = 1:size(images_pulse_ifft_avg_tot,1)-1
    N_start = 3300;%1500;
    N_stop = size(images_pulse_ifft_avg_tot,2)-100;
    % fit en log
    if pic(v)-300 > N_start
    c = polyfit(x([N_start:pic(v)-300 pic(v)+300:N_stop]),log(images_pulse_ifft_avg_tot(v,[N_start:pic(v)-300 pic(v)+300:N_stop])),1);
    else
    c = polyfit(x([N_start:N_stop]),log(images_pulse_ifft_avg_tot(v,[N_start:N_stop])),1);
    end
    y_est = polyval(c,x);
    % hold on; plot(x,log(images_pulse_ifft_avg_tot(v,:))); plot(x,y_est);
    %apply the fit in non-log
    E = exp(y_est);
    if plotflag
        clf
        hold on; plot(x,images_pulse_ifft_avg_tot(v,:)); plot(x,exp(y_est));
        plot(pic(v),mean(E(pic(v)-N_mean:pic(v)+N_mean)),'.','MarkerSize', 20)
        pause(1)
    end
    %measure the ratio
    ratio(v) = mean(images_pulse_ifft_avg_tot(v,pic(v)-N_mean:pic(v)+N_mean))/ (mean(E(pic(v)-N_mean:pic(v)+N_mean)));
end

end