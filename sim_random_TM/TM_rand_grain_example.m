%% 

clearvars

% find out which user
[~, user] = system('id -F');
user = user(1:(end-1)); % remove linebreak
if strcmp(user, 'Bernhard')
    git_repo_dir    = '/Users/bernhard/Documents/Projects/Pulses in random media/data_analysis/_Codes';
elseif strcmp(user, 'Louisiane') %??? need to enter username and path
    git_repo_dir    = '???';
end

% add git repo
run([git_repo_dir, '/add_folder_to_search_path.m'])

%%

n_slm = 32^2;
n_ccd = 32^2;

n_real = 5;

bin = 1;

grain_size_TM = 10;
form_flag = 2;

svd_idx = 1;

tic

kd_rand = zeros(sqrt(n_ccd),sqrt(n_ccd));
ca_rand = zeros(sqrt(n_ccd),sqrt(n_ccd));
out_focus = zeros(sqrt(n_ccd),sqrt(n_ccd));
kd_svd = zeros(length(svd_idx),sqrt(n_ccd),sqrt(n_ccd));
ca_svd = zeros(length(svd_idx),sqrt(n_ccd),sqrt(n_ccd));
kd_svd_po = zeros(length(svd_idx),sqrt(n_ccd),sqrt(n_ccd));
ca_svd_po = zeros(length(svd_idx),sqrt(n_ccd),sqrt(n_ccd));
field_dist_svd = [];
field_dist_rand = [];
avg_svd_int = zeros(sqrt(n_ccd),sqrt(n_ccd));
avg_rand_int = zeros(sqrt(n_ccd),sqrt(n_ccd));
fprintf('\n');
for idx_real = 1:n_real
   
    TM = randTM(n_ccd, n_slm, grain_size_TM, form_flag, 0);
    
    [TM_bin, n_ccd_bin] = bin_TM(TM,bin);
    
    [TM_filt, TM_pic_fft] = filter_TM(TM_bin, 2, 1);
    
    % random input
    out_rand = reshape(TM*rand_vector(n_slm), sqrt(n_ccd), sqrt(n_ccd)); 
    kd_rand = kd_rand + k_dist(out_rand)/n_real;
    ca_rand = ca_rand + corr_auto(out_rand)/n_real;

    % focus
    [focus_vec, ~] = focus_center(TM_filt);
    out_focus = out_focus + reshape(TM*focus_vec, sqrt(n_ccd), sqrt(n_ccd))/n_real;

    % svd
    for idx = 1:length(svd_idx)
        
        % full control
        svd_vec = svd_vector(TM_filt,svd_idx(idx));
        out_svd = reshape(TM*svd_vec, sqrt(n_ccd), sqrt(n_ccd));
        kd_svd(idx,:,:) = squeeze( kd_svd(idx,:,:) ) + k_dist(out_svd)/n_real;
        ca_svd(idx,:,:) = squeeze( ca_svd(idx,:,:) ) + corr_auto(out_svd)/n_real;
        
        % phase only
        svd_po_vec = exp(1i*angle(svd_vec));
        out_svd_po = reshape(TM*svd_po_vec, sqrt(n_ccd), sqrt(n_ccd));
        kd_svd_po(idx,:,:) = squeeze( kd_svd_po(idx,:,:) ) + k_dist(out_svd_po)/n_real;
        ca_svd_po(idx,:,:) = squeeze( ca_svd_po(idx,:,:) ) + corr_auto(out_svd_po)/n_real;
        
    end
    
    fprintf('real %d / %d finished!\n',idx_real,n_real);
    
end

toc

%save('data_randTM_grain_N10_g2_52vec','kd_svd','ca_svd','kd_svd_po','ca_svd_po','kd_rand','ca_rand','out_focus','svd_idx','grain_size_TM','n_slm','n_ccd')

%% plot overall TM k-dist

figure(89)
clf
set(gcf,'defaultAxesTickLabelInterpreter','latex',...
        'defaulttextinterpreter','latex',...
        'defaultLegendInterpreter','latex');
box on
imagesc(TM_pic_fft)
xlabel('$x$')
ylabel('$y$')
set(gca, 'FontSize' , 18)

%saveas(gcf,'./elongated_grains_ext/TM_out_fft_antifocus.png')

%% show example speckles

x = (1:sqrt(n_ccd)) - sqrt(n_ccd)/2 - 1;

figure(12)
clf
set(gcf,'defaultAxesTickLabelInterpreter','latex',...
        'defaulttextinterpreter','latex',...
        'defaultLegendInterpreter','latex');
set(gcf, 'units', 'centimeters', 'position', [3 3 30 12])

subplot(1,2,1)
imagesc(x,x,abs(out_rand))
%imagesc(x,x,avg_rand_int)
xlabel('$x$')
ylabel('$y$')
title('random')
set(gca, 'FontSize' , 18)

subplot(1,2,2)
imagesc(x,x,abs(out_svd_po))
%imagesc(x,x,avg_svd_int)
xlabel('$x$')
ylabel('$y$')
title('1st svd')
set(gca, 'FontSize' , 18)

%saveas(gcf,'svd_last_vs_random_specle_ex_flatOTF.png')

%% plot k-dist

k = (1:sqrt(n_ccd)) - sqrt(n_ccd)/2 - 1;

svd_plot_idx = 1;

figure(11)
clf
set(gcf,'defaultAxesTickLabelInterpreter','latex',...
        'defaulttextinterpreter','latex',...
        'defaultLegendInterpreter','latex');
set(gcf, 'units', 'centimeters', 'position', [3 3 40 11])

subplot(1,3,1)
imagesc(k,k,abs(kd_rand))
xlabel('$k_x$')
ylabel('$k_y$')
title('random')
set(gca, 'FontSize' , 18)

subplot(1,3,2)
%imagesc(k,k,abs(squeeze(kd_svd(svd_plot_idx,:,:))))
imagesc(k,k,abs(squeeze(kd_svd_po(svd_plot_idx,:,:))))
xlabel('$k_x$')
ylabel('$k_y$')
title(sprintf('%d / %d svd',svd_idx(svd_plot_idx),n_ccd))
set(gca, 'FontSize' , 18)

k_mean = (sqrt(n_ccd)/2):(sqrt(n_ccd)/2+2);
subplot(1,3,3)
box on
hold on
%plot(k,mean(abs(squeeze(kd_svd(svd_plot_idx,:,:)))) / sum(mean(abs(squeeze(kd_svd(svd_plot_idx,:,:))))), 'linewidth', 2)
plot(k,mean(abs(squeeze(kd_svd_po(svd_plot_idx,k_mean,:))),1) / sum(mean(abs(squeeze(kd_svd_po(svd_plot_idx,k_mean,:))),1)), 'linewidth', 2)
plot(k,mean(abs(squeeze(kd_svd_po(svd_plot_idx,:,k_mean))),2) / sum(mean(abs(squeeze(kd_svd_po(svd_plot_idx,:,k_mean))),2)), 'linewidth', 2)
plot(k,mean(abs(kd_rand)) / sum(mean(abs(kd_rand))), 'linewidth', 2)
xlabel('$k_x$')
ylabel('$A_k$')
xlim([k(1) k(end)])
title('central cut')
legend('x','y','location','best')
set(gca, 'FontSize' , 18)

%saveas(gcf,'./elongated_grains_ext/k-dist_svd1_simple.png')

%% plot auto-correlation

svd_plot_idx = 1;

x = (1:sqrt(n_ccd)) - sqrt(n_ccd)/2 - 1;
center = sqrt(n_ccd)/2 + 1;
center_foc = center-1;

figure(10)
clf
set(gcf,'defaultAxesTickLabelInterpreter','latex',...
        'defaulttextinterpreter','latex',...
        'defaultLegendInterpreter','latex');
set(gcf, 'units', 'centimeters', 'position', [3 3 40 11])

subplot(1,3,1)
imagesc(x,x,abs(ca_rand))
xlabel('$x$')
ylabel('$y$')
title('random')
set(gca, 'FontSize' , 18)

subplot(1,3,2)
%imagesc(x,x,abs(squeeze(ca_svd(svd_plot_idx,:,:))))
imagesc(x,x,abs(squeeze(ca_svd_po(svd_plot_idx,:,:))))
xlabel('$x$')
ylabel('$y$')
title('1st svd')
set(gca, 'FontSize' , 18)

subplot(1,3,3)
box on
hold on
%plot(x,abs(squeeze(ca_svd(svd_plot_idx,center,:))) / abs(squeeze(ca_svd(svd_plot_idx,center,center))), 'linewidth', 2)
plot(x,abs(squeeze(ca_svd_po(svd_plot_idx,:,center))) / abs(squeeze(ca_svd_po(svd_plot_idx,center,center))), 'linewidth', 2)
plot(x,abs(squeeze(ca_svd_po(svd_plot_idx,center,:))) / abs(squeeze(ca_svd_po(svd_plot_idx,center,center))), 'linewidth', 2)
plot(x,abs(ca_rand(center,:)) / abs(ca_rand(center,center)), 'linewidth', 2)
%plot(x,circshift(abs(out_focus(center_foc,:)),1) / abs(out_focus(center_foc,center_foc)), 'linewidth', 2)
xlabel('$x$')
ylabel('auto corr norm')
xlim([x(1) x(end)])
legend('x','y','location','best')
%title('1st svd')
set(gca, 'FontSize' , 18)

%saveas(gcf,'./elongated_grains_ext/autocorr_svd1_po_antifocus.png')

%% fit grain size

fwhm_ca_svd_x = zeros(1,length(svd_idx));
fwhm_ca_svd_po_x = zeros(1,length(svd_idx));
fwhm_ca_svd_y = zeros(1,length(svd_idx));
fwhm_ca_svd_po_y = zeros(1,length(svd_idx));
for idx = 1:length(svd_idx)
    
    norm_ca_svd_slice = abs(squeeze(ca_svd(idx,center,:))) / abs(squeeze(ca_svd(idx,center,center)));
    norm_ca_svd_po_slice = abs(squeeze(ca_svd_po(idx,center,:))) / abs(squeeze(ca_svd_po(idx,center,center)));

%     [~, min_idx_l] = min(abs(norm_ca_svd_slice(1:center) - 0.5));
%     
%     fwhm_ca_svd(idx) = 2*abs(x(min_idx_l));
    
    ftype = fittype('b + (1-b)*exp(-x.^2 / (2*a^2))','coefficients',{'a','b'});
    fo = fitoptions('Method', 'NonlinearLeastSquares','StartPoint', [grain_size_TM / 2*sqrt(2*log(2)), 0]);
    ca_fit = fit(x',norm_ca_svd_slice,ftype,fo);
    ca_po_fit = fit(x',norm_ca_svd_po_slice,ftype,fo);

    fwhm_ca_svd_x(idx) = 2*sqrt(2*log(2)) * ca_fit.a;
    fwhm_ca_svd_po_x(idx) = 2*sqrt(2*log(2)) * ca_po_fit.a;
    
    
    norm_ca_svd_slice = abs(squeeze(ca_svd(idx,:,center))) / abs(squeeze(ca_svd(idx,center,center)));
    norm_ca_svd_po_slice = abs(squeeze(ca_svd_po(idx,:,center))) / abs(squeeze(ca_svd_po(idx,center,center)));

%     [~, min_idx_l] = min(abs(norm_ca_svd_slice(1:center) - 0.5));
%     
%     fwhm_ca_svd(idx) = 2*abs(x(min_idx_l));
    
    ftype = fittype('b + (1-b)*exp(-x.^2 / (2*a^2))','coefficients',{'a','b'});
    fo = fitoptions('Method', 'NonlinearLeastSquares','StartPoint', [grain_size_TM / 2*sqrt(2*log(2)), 0]);
    ca_fit = fit(x',norm_ca_svd_slice',ftype,fo);
    ca_po_fit = fit(x',norm_ca_svd_po_slice',ftype,fo);

    fwhm_ca_svd_y(idx) = 2*sqrt(2*log(2)) * ca_fit.a;
    fwhm_ca_svd_po_y(idx) = 2*sqrt(2*log(2)) * ca_po_fit.a;

%     figure(5)
%     clf
%     box on
%     hold on
%     plot(x,norm_ca_svd_slice)
%     plot(ca_fit)
%     xlim([x(1) x(end)])
%     
%     pause(0.2)
end

norm_ca_rand_slice = abs(squeeze(ca_rand(center,:))) / abs(squeeze(ca_rand(center,center)));

% [~, min_idx_l] = min(abs(norm_ca_rand_slice(1:center) - 0.5));
% fwhm_ca_rand = 2*abs(x(min_idx_l));

ftype = fittype('b + (1-b)*exp(-x.^2 / (2*a^2))','coefficients',{'a','b'});
fo = fitoptions('Method', 'NonlinearLeastSquares','StartPoint', [grain_size_TM / 2*sqrt(2*log(2)), 0]);
ca_fit = fit(x',norm_ca_rand_slice',ftype,fo);
fwhm_ca_rand = 2*sqrt(2*log(2)) * ca_fit.a;


%% plot grain size

figure(99)
clf
set(gcf,'defaultAxesTickLabelInterpreter','latex',...
        'defaulttextinterpreter','latex',...
        'defaultLegendInterpreter','latex');
set(gcf, 'units', 'centimeters', 'position', [3 3 20 15])

box on
hold on
h1 = plot(svd_idx, fwhm_ca_svd_po_x,'o');
h2 = plot(svd_idx, fwhm_ca_svd_po_y,'o');
plot([svd_idx(1) svd_idx(end)], grain_size_TM*[1 1],'k-','linewidth',2)
%plot([svd_idx(1) svd_idx(end)], fwhm_ca_rand*[1 1],'--')
set(h1,'MarkerFaceColor',get(h1,'Color'))
set(h2,'MarkerFaceColor',get(h2,'Color'))
xlim([svd_idx(1) svd_idx(end)])
xlabel('svd vector number')
ylabel('grain size (px)')
%title('1st svd')
legend('phase only X','phase only Y','location','best')
set(gca, 'FontSize' , 18)

%saveas(gcf,'svd_grain_size.png')

%% k dist video

for idx = 3%length(svd_idx)
    
    figure(11)
    clf
    set(gcf,'defaultAxesTickLabelInterpreter','latex',...
        'defaulttextinterpreter','latex',...
        'defaultLegendInterpreter','latex');
    set(gcf, 'units', 'centimeters', 'position', [3 3 40 11])
    
    subplot(1,3,1)
    imagesc(k,k,abs(kd_rand))
    xlabel('$k_x$')
    ylabel('$k_y$')
    title('random')
    set(gca, 'FontSize' , 18)
    
    subplot(1,3,2)
    imagesc(k,k,abs(squeeze(kd_svd(idx,:,:))))
    xlabel('$k_x$')
    ylabel('$k_y$')
    title('1st svd')
    title(sprintf('svd %d / %d',svd_idx(idx),n_ccd))
    set(gca, 'FontSize' , 18)
    
    subplot(1,3,3)
    box on
    hold on
    plot(k,mean(abs(squeeze(kd_svd(idx,:,:)))) / sum(mean(abs(squeeze(kd_svd(idx,:,:))))), 'linewidth', 2)
    plot(k,mean(abs(kd_rand)) / sum(mean(abs(kd_rand))), 'linewidth', 2)
    xlabel('$k_x$')
    ylabel('$A_k$')
    xlim([k(1) k(end)])
    %title('1st svd')
    set(gca, 'FontSize' , 18)

    F(idx) = getframe(gcf);
    
    pause(0.33)
end

%write_movie(F, 'svd_k_dist', 3);

%% statistics

x = linspace(0,9,100);
dist_rylgh = exp(-x);

figure(112)
clf
set(gcf,'defaultAxesTickLabelInterpreter','latex',...
    'defaulttextinterpreter','latex',...
    'defaultLegendInterpreter','latex');
set(gcf, 'units', 'centimeters', 'position', [3 3 20 15])

box on
hold on
histogram(abs(field_dist_svd(:)).^2 / mean(abs(field_dist_svd(:)).^2),'Normalization','pdf', 'DisplayStyle', 'stairs', 'linewidth',2)
histogram(abs(field_dist_rand(:)).^2 / mean(abs(field_dist_rand(:)).^2),'Normalization','pdf', 'DisplayStyle', 'stairs', 'linewidth',2)
plot(x, dist_rylgh, 'k--','linewidth',2)
set(gca,'YScale','log')
xlim([0.01 9])
xlabel('$I / \langle I \rangle$')
ylabel('pdf')
legend('svd','rand','rayleigh','location','northeast')
title('100 realizations (1024$\times$1024)')
set(gca, 'FontSize' , 18)

%saveas(gcf,'randTMg_int_dist.png')
