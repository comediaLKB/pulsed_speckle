function [s_mp, rho_mp, s_norm] = MP_law(Image, flag)

%code that plot the MP probalibity distribution and singular value
%histogram.
% the difference with the wikipedia formula comes from the fact that
% wikipedia does add all zeros to the singular values if one dimention is
% bigger than the other one.

[U S V] = svd(Image);
s = [diag(S)];
s_norm  = s./sqrt(sum(s.^2)/size(S,1));

gamma = size(S,2) / size(S,1);

%MP law
s_min = abs(1-sqrt(1/gamma));
s_max = abs(1+sqrt(1/gamma));
s_mp = linspace(s_min,s_max,1000);
rho_mp = max(1,gamma)./(pi*s_mp) .* sqrt((s_mp.^2 - s_min.^2) .* (s_max.^2 - s_mp.^2));

if flag
    figure(1)
    set(gca,'fontsize',25)
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    hold on
    box on
    plot(s_mp, rho_mp, '.')
    histogram(s_norm,100, 'Normalization', 'pdf')
    xlabel('s')
    ylabel('P(s)')
    title('Marchenko-Pastur law')
    
    figure(2)
    set(gca,'fontsize',25)
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    hold on
    box on
    plot(s_mp.^2, rho_mp./s_mp/2, '.')
    histogram(s_norm.^2,100, 'Normalization', 'pdf')
    xlabel('$s^2$')
    ylabel('P($s^2$)')
    title('Marchenko-Pastur law')
end

end