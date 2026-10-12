
% =====parameters =====
bg_CH4 = 1922.6; d14bg = 350; d13bg = -47.6; dDbg = -99; ethanebg = 1.5435; % ppb
d14ff = -1000; d14apm = -799; d14wet = 37.7;
d13ff = -44.8;  d13apm = -75.8; d13wt = -67.8;
ethaneff = 0.0652;                    
%% Uniform source-signature ranges
% d14C [fossil, permafrost, wetland]
d14_min = [-1000, -944, -28.1];
d14_max = [-1000, -655, 103.5];

% d13C [fossil, permafrost, wetland]
d13_min = [-51, -81.6, -83];
d13_max = [-35, -61.2, -55];

% Ethane:CH4 fossil ratio
ethane_min = 0.009;
ethane_max = 0.20;

unc_d14apm_age = 125;
d14_age_halfwidth = sqrt(3) * unc_d14apm_age;
%noise
noise_ch4=3; noise_d14C=5; noise_d13C=0.2; noise_eth=0.07;
K = 500; rng(1);
den_floor = 10;           
max_resamples = 5;      
cond_thresh = 1e-10;    
lambda_ridge = 1e-6;       
d14vec_mean = [d14ff, d14apm, d14wet];
d13vec_mean = [d13ff, d13apm, d13wt];
ethvec_mean = [ethaneff, 0, 0];
d14_min_vec = d14_min;
d14_max_vec = d14_max;
d13_min_vec = d13_min;
d13_max_vec = d13_max;
A_mean = [ 1, 1, 1;
           d14vec_mean;
           d13vec_mean;
           1000*ethvec_mean ];         % ethane row scaled by 1000
w = [ 1/noise_ch4^2; 1/noise_d14C^2; 1/noise_d13C^2; 1/noise_eth^2 ];
Wsqrt = diag(sqrt(w));
solve_wls_stable = @(A,b) stable_wls(A,b,Wsqrt,cond_thresh,lambda_ridge);



pf_max   = 500;   pf_step = 10;      % Permafrost CH4 grid (ppb)
EX_max   = 700;   EX_step = 10;      % Total excess CH4 grid (ppb)
frac_step = 0.02; 

pf_vec = (0:pf_step:pf_max)';    
EX_vec = (0:EX_step:EX_max)';   

pf_max_plot = max(pf_vec);
EX_max_plot = max(EX_vec);
fracPF_lines = [0.2 0.4 0.6 0.8 1.0];  

[half_ppb_50, half_pct_50] = run_one_fraction_total(0.2, pf_vec, EX_vec, ...
    bg_CH4, d14bg, d13bg, ethanebg, ...
    d14vec_mean, d13vec_mean, ethaneff, ...
    d14_min_vec, d14_max_vec, d13_min_vec, d13_max_vec, ethane_min, ethane_max, ...
    noise_ch4, noise_d14C, noise_d13C, noise_eth, ...
    K, den_floor, max_resamples, Wsqrt, cond_thresh, lambda_ridge);



%Blue to Dark Orange, 7 steps                                           

colourpalette = [0.000   0.400   0.400
    0.220   0.664   0.700
    0.600   0.940   0.980
    0.900   1.000   1.000
    1.000   0.793   0.600
    1.000   0.560   0.200
    0.600   0.250   0.000];
%
sat_factor = 0.65;   % 1 = original, 0 = grayscale
cmap_hsv = rgb2hsv(colourpalette);

cmap_hsv(:,2) = cmap_hsv(:,2) * sat_factor;  % reduce saturation only

colourpalette = hsv2rgb(cmap_hsv);




figure('Color','w');
tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

nexttile;
Mppb = half_ppb_50;
Mppb(Mppb == 0) = NaN;
imagesc(pf_vec, EX_vec, Mppb','AlphaData', ~isnan(Mppb'));
axis xy tight
colormap(colourpalette);
%xlabel('Permafrost CH_4 (ppb)');
ylabel('Total excess CH_4 (ppb)');
%cb1 = colorbar('Location','northoutside'); ylabel(cb1,'Uncertainty (ppb)');
set(gca,'FontSize',20);
caxis([0, 100])
cb1 = colorbar;
cb1.Location = 'northoutside';        % puts it above the tile
cb1.Label.String = 'Absolute Uncertainty (ppb)';
cb1.Label.FontSize = 20;              % title size
cb1.FontSize = 20;                    % tick size
cb1.Ticks = [0 20 40 60 80];
set(gca, 'XTick', [0 50 100 150 200 250 300 350 400 450 500]);       set(gca, 'XTickLabel', {'', '', '200', '', '400', ''});   grid on;
set(gca, 'XTickLabel', []);
hold on;
for fp = fracPF_lines
    x_line = linspace(1, pf_max_plot, 400);
    y_line = x_line ./ fp;
    mask = (y_line >= 0) & (y_line <= EX_max_plot) & (y_line >= x_line);
    if ~any(mask), continue; end

    if abs(fp-1.0) < 1e-6
        plot(x_line(mask), y_line(mask), 'k-', 'LineWidth', 1.8);
        label_str = '100%';
    else
        plot(x_line(mask), y_line(mask), 'k--', 'LineWidth', 0.8);
        label_str = sprintf('%d%%', round(fp*100));
    end

    idx_label = round(0.7 * nnz(mask));
    if idx_label > 0
        tx = x_line(mask);
        ty = y_line(mask);
        text(tx(idx_label), ty(idx_label), label_str, ...
            'Color','k','FontSize',20, ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom');
    end
end

line1 = plot(nan,nan,'k--','LineWidth',0.8);
lg=legend(line1, {'% Permafrost CH_4'}, 'Location','southeast');
lg.Box = 'off';
lg.Color = 'none';
lg.FontSize = 20;
hold off;


nexttile;
hold on
Mpct = half_pct_50;
Mpct(Mpct == 0) = NaN;
imagesc(pf_vec, EX_vec, Mpct','AlphaData', ~isnan(Mppb'));
axis xy tight
colormap(colourpalette);
caxis([0 105]);
cb2 = colorbar;
cb2.Ticks = [0 15 30 45 60 75 90 100];
cb2.TickLabels = string(cb2.Ticks);
cb2.Location = 'northoutside';
cb2.Label.String = 'Percentage Uncertainty (%)';
cb2.Label.FontSize = 20;
cb2.FontSize = 20;
set(gca,'FontSize',20);
set(gca, 'YTickLabel', []);
set(gca, 'XTick', [0 50 100 150 200 250 300 350 400 450 500]);       
set(gca, 'XTickLabel', []);
for fp = fracPF_lines
    x_line = linspace(1, pf_max_plot, 400);
    y_line = x_line ./ fp;
    mask = (y_line >= 0) & (y_line <= EX_max_plot) & (y_line >= x_line);
    if ~any(mask), continue; end

    if abs(fp-1.0) < 1e-6
        plot(x_line(mask), y_line(mask), 'k-', 'LineWidth', 1.8);
        label_str = '100%';
    else
        plot(x_line(mask), y_line(mask), 'k--', 'LineWidth', 0.8);
        label_str = sprintf('%d%%', round(fp*100));
    end

    idx_label = round(0.7 * nnz(mask));
    if idx_label > 0
        tx = x_line(mask);
        ty = y_line(mask);
        text(tx(idx_label), ty(idx_label), label_str, ...
            'Color','k','FontSize',20, ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom');
    end
end
set(gca,'Layer','top')
grid on;
box on;
hold off;

FW_fixed = 100;                      
pf_vec = (0:pf_step:pf_max)';
fracTot_vec = (0:frac_step:1.0)';      
nPF = numel(pf_vec); nFracTot = numel(fracTot_vec);

PF_true = zeros(nPF,nFracTot);
PF_det_meas_ppb = nan(nPF,nFracTot);
PF_det_meas_pct = nan(nPF,nFracTot);
PF_rec_mean = nan(nPF,nFracTot);
PF_rec_p16  = nan(nPF,nFracTot);
PF_rec_p84  = nan(nPF,nFracTot);
PF_width1684_ppb = nan(nPF,nFracTot);
PF_width1684_pct = nan(nPF,nFracTot);
Bias_mean   = nan(nPF,nFracTot);

for ii = 1:nPF
    PF = pf_vec(ii);
    EX = PF + FW_fixed;             
    CH4_true = bg_CH4 + EX;

    for jj = 1:nFracTot
        fTot = fracTot_vec(jj);

        FF = fTot * EX;
        WT = FW_fixed - FF;

        if WT < -1e-12 || FF > FW_fixed + 1e-12
            continue;
        end
        WT = max(WT,0);

        PF_true(ii,jj) = PF;

        d14_obs = (d14bg*bg_CH4 + d14ff*FF + d14apm*PF + d14wet*WT)/CH4_true;
        d13_obs = (d13bg*bg_CH4 + d13ff*FF + d13apm*PF + d13wt*WT)/CH4_true;
        eth_obs = ethanebg + 1000*ethaneff*FF;


        CH4_meas = CH4_true + noise_ch4*randn;
        d14_meas = d14_obs + noise_d14C*randn;
        d13_meas = d13_obs + noise_d13C*randn;
        eth_meas = eth_obs + noise_eth*randn;
        b_meas=[CH4_meas-bg_CH4;
                d14_meas*CH4_meas - d14bg*bg_CH4;
                d13_meas*CH4_meas - d13bg*bg_CH4;
                eth_meas - ethanebg];
        x_meas = solve_wls_stable(A_mean,b_meas);
        x_meas = max(x_meas,0);
        diff_ppb = abs(x_meas(2)-PF);
        PF_det_meas_ppb(ii,jj)=diff_ppb;
        PF_det_meas_pct(ii,jj)=100*diff_ppb/max(PF,den_floor);

        x_pf = nan(1,K);
        for k=1:K
            CH4_k=CH4_true + noise_ch4*randn;
            d14_k=d14_obs + noise_d14C*randn;
            d13_k=d13_obs + noise_d13C*randn;
            eth_k=eth_obs + noise_eth*randn;
            b_k=[CH4_k-bg_CH4;
                 d14_k*CH4_k - d14bg*bg_CH4;
                 d13_k*CH4_k - d13bg*bg_CH4;
                 eth_k - ethanebg];

            d14v = d14_min_vec + (d14_max_vec-d14_min_vec).*rand(1,3);
            d13v = d13_min_vec + (d13_max_vec-d13_min_vec).*rand(1,3);
            ethv1 = ethane_min + (ethane_max-ethane_min)*rand;

            A_k=[1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
            x_draw=solve_wls_stable(A_k,b_k);

            if any(~isfinite(x_draw)) || any(abs(x_draw)>1e6)
                tried=1;
                while tried<max_resamples && (~all(isfinite(x_draw)) || any(abs(x_draw)>1e6))
                    d14v = d14_min_vec + (d14_max_vec-d14_min_vec).*rand(1,3);
                    d13v = d13_min_vec + (d13_max_vec-d13_min_vec).*rand(1,3);
                    ethv1 = ethane_min + (ethane_max-ethane_min)*rand;
                    A_k=[1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                    x_draw=solve_wls_stable(A_k,b_k);
                    tried=tried+1;
                end
            end

            x_pf(k)=max(x_draw(2),0);
        end

        PF_rec_mean(ii,jj)=mean(x_pf,'omitnan');
        PF_rec_p16(ii,jj)=prctile(x_pf,16);
        PF_rec_p84(ii,jj)=prctile(x_pf,84);
        PF_width1684_ppb(ii,jj)=PF_rec_p84(ii,jj)-PF_rec_p16(ii,jj);
        PF_width1684_pct(ii,jj)=100*PF_width1684_ppb(ii,jj)/max(PF,den_floor);
        Bias_mean(ii,jj)=PF_rec_mean(ii,jj)-PF;
    end
end

PF_det_meas_pct(PF_det_meas_pct>100)=100; PF_det_meas_pct(PF_det_meas_pct<0)=0;
PF_width1684_pct(PF_width1684_pct>100)=100; PF_width1684_pct(PF_width1684_pct<0)=0;

% Plot absolute uncertainty
nexttile
Mppb_fracTot=(PF_width1684_ppb/2)';
imagesc(pf_vec,100*fracTot_vec,Mppb_fracTot,'AlphaData',~isnan(Mppb_fracTot)); axis xy tight
colormap(colourpalette);
ylabel('Fossil Fraction (%)');
caxis([0 100]);
set(gca,'FontSize',20);
set(gca,'XTick',[0 50 100 150 200 250 300 350 400 450 500]);
set(gca,'XTickLabel',[]);
set(gca,'YTick',0:20:100);
set(gca,'YTickLabel',{'0','20','40','60','80','100'});
set(gca,'Color','w');

xline_pf=linspace(0,pf_max_plot,500);
ymax_tot=100*FW_fixed./(xline_pf+FW_fixed);
hold on
plot(xline_pf,ymax_tot,'k-','LineWidth',1.6);
hold off
grid on

% Plot percentage uncertainty
nexttile
Mpct_fracTot=(PF_width1684_pct/2)';
imagesc(pf_vec,100*fracTot_vec,Mpct_fracTot,'AlphaData',~isnan(Mpct_fracTot)); axis xy tight
colormap(colourpalette);
caxis([0 105]);
set(gca,'FontSize',20);
set(gca,'YTickLabel',[]);
set(gca,'XTick',[0 50 100 150 200 250 300 350 400 450 500]);
set(gca,'XTickLabel',[]);
set(gca,'YTick',0:20:100);
set(gca,'Color','w');
hold on
plot(xline_pf,ymax_tot,'k-','LineWidth',1.6);
hold off
grid on


%for d14

pf_max = 500; pf_step = 10;
FF_fixed = 20; WT_fixed = 80;

d14apm_min = -1000; d14apm_max = 0; d14apm_step = 30;  

den_floor = 10;
max_resamples = 5; cond_thresh = 1e-10; lambda_ridge = 1e-6;

w = [ 1/noise_ch4^2; 1/noise_d14C^2; 1/noise_d13C^2; 1/noise_eth^2 ];
Wsqrt = diag(sqrt(w));
solve_wls_stable = @(A,b) stable_wls(A,b,Wsqrt,cond_thresh,lambda_ridge);

pf_vec    = (0:pf_step:pf_max)';                      
d14apm_vec = (d14apm_min:d14apm_step:d14apm_max)';    
nPF = numel(pf_vec); nD14 = numel(d14apm_vec);

PF_true = zeros(nPF,nD14);
PF_det_meas_ppb = nan(nPF,nD14);  
PF_det_meas_pct = nan(nPF,nD14);  
PF_rec_mean = nan(nPF,nD14);
PF_rec_p16  = nan(nPF,nD14);
PF_rec_p84  = nan(nPF,nD14);
PF_width1684_ppb = nan(nPF,nD14);
PF_width1684_pct = nan(nPF,nD14);
Bias_mean   = nan(nPF,nD14);

for ii = 1:nPF
    PF = pf_vec(ii);
    for jj = 1:nD14
        d14apm_this = d14apm_vec(jj);
        FF = FF_fixed; WT = WT_fixed;
        PF_true(ii,jj) = PF;
        CH4_true = bg_CH4 + FF + PF + WT;

        d14_obs = (d14bg*bg_CH4 + d14ff*FF + d14apm_this*PF + d14wet*WT)/CH4_true;
        d13_obs = (d13bg*bg_CH4 + d13ff*FF + d13apm*PF      + d13wt*WT)/CH4_true;
        eth_obs = ethanebg + 1000*ethaneff*FF;

        d14vec_mean = [d14ff, d14apm_this, d14wet];
        d13vec_mean = [d13ff, d13apm,      d13wt];
        ethvec_mean = [ethaneff, 0, 0];
        A_mean = [1,1,1; d14vec_mean; d13vec_mean; 1000*ethvec_mean];


        CH4_meas = CH4_true + noise_ch4*randn;
        d14_meas = d14_obs + noise_d14C*randn;
        d13_meas = d13_obs + noise_d13C*randn;
        eth_meas = eth_obs + noise_eth*randn;
        b_meas = [ CH4_meas - bg_CH4;
                   d14_meas*CH4_meas - d14bg*bg_CH4;
                   d13_meas*CH4_meas - d13bg*bg_CH4;
                   eth_meas - ethanebg ];
        x_meas = solve_wls_stable(A_mean,b_meas);
        x_meas = max(x_meas,0);
        diff_ppb = abs(x_meas(2) - PF);
        PF_det_meas_ppb(ii,jj) = diff_ppb;
        PF_det_meas_pct(ii,jj) = 100 * diff_ppb / max(PF, den_floor);

        x_pf = nan(1,K);
        for k=1:K
            CH4_k=CH4_true + noise_ch4*randn;
            d14_k=d14_obs + noise_d14C*randn;
            d13_k=d13_obs + noise_d13C*randn;
            eth_k=eth_obs + noise_eth*randn;
            b_k = [ CH4_k - bg_CH4;
                    d14_k*CH4_k - d14bg*bg_CH4;
                    d13_k*CH4_k - d13bg*bg_CH4;
                    eth_k - ethanebg ];

            d14ff_k = -1000;
            d14apm_k = (d14apm_this - d14_age_halfwidth) + ...
                       (2*d14_age_halfwidth)*rand;
            d14wet_k = d14_min_vec(3) + ...
                       (d14_max_vec(3)-d14_min_vec(3))*rand;
            d14v = [d14ff_k, d14apm_k, d14wet_k];

            d13v = d13_min_vec + (d13_max_vec-d13_min_vec).*rand(1,3);
            ethv1 = ethane_min + (ethane_max-ethane_min)*rand;

            A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
            x_draw = solve_wls_stable(A_k,b_k);
            if any(~isfinite(x_draw)) || any(abs(x_draw)>1e6)
                tried=1;
                while tried<max_resamples && (~all(isfinite(x_draw)) || any(abs(x_draw)>1e6))
                    d14ff_k = -1000;
                    d14apm_k = (d14apm_this - d14_age_halfwidth) + ...
                               (2*d14_age_halfwidth)*rand;
                    d14wet_k = d14_min_vec(3) + ...
                               (d14_max_vec(3)-d14_min_vec(3))*rand;
                    d14v = [d14ff_k, d14apm_k, d14wet_k];

                    d13v = d13_min_vec + (d13_max_vec-d13_min_vec).*rand(1,3);
                    ethv1 = ethane_min + (ethane_max-ethane_min)*rand;
                    A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                    x_draw = solve_wls_stable(A_k,b_k);
                    tried=tried+1;
                end
            end
            x_pf(k)=max(x_draw(2),0);
        end

        PF_rec_mean(ii,jj) = mean(x_pf,'omitnan');
        PF_rec_p16(ii,jj)  = prctile(x_pf,16);
        PF_rec_p84(ii,jj)  = prctile(x_pf,84);
        PF_width1684_ppb(ii,jj) = PF_rec_p84(ii,jj) - PF_rec_p16(ii,jj);
        PF_width1684_pct(ii,jj) = 100 * PF_width1684_ppb(ii,jj) / max(PF, den_floor);
        Bias_mean(ii,jj)   = PF_rec_mean(ii,jj) - PF;
    end
end


PF_det_meas_pct(PF_det_meas_pct>100)=100; PF_det_meas_pct(PF_det_meas_pct<0)=0;
PF_width1684_pct(PF_width1684_pct>100)=100; PF_width1684_pct(PF_width1684_pct<0)=0;

nexttile 
imagesc(pf_vec, d14apm_vec, (PF_width1684_ppb/2)'); axis xy tight
colormap(colourpalette);
xlabel('Permafrost CH_4 (ppb)'); ylabel('Permafrost \Delta^{14}C (‰)');
%title('MC 16–84% range width [ppb]'); 
%ylabel(cb,'Uncertainty (ppb)');
caxis([0 100]);
set(gca,'FontSize',20);
set(gca, 'XTick', [0 50 100 150 200 250 300 350 400 450 500]);  set(gca, 'XTickLabel', {'0', '','100', '','200','', '300', '','400','', '500'}); 
set(gca, "YTick",[-1000 -900 -800 -700 -600 -500 -400 -300 -200 -100 0])
set(gca, 'YTickLabel', {'-1000', '','-800', '','-600','', '-400', '','-200','', '0'})
grid on

nexttile
hold on
imagesc(pf_vec, d14apm_vec, (PF_width1684_pct/2)'); axis xy tight
colormap(colourpalette); 
caxis([0 105]);
xlabel('Permafrost CH_4 (ppb)');% ylabel('Permafrost \delta^{14}C (‰)');
set(gca,'FontSize',20);
set(gca, 'YTickLabel', []);
set(gca, 'XTick', [0 50 100 150 200 250 300 350 400 450 500]); set(gca, 'XTickLabel', {'0', '','100', '','200','', '300', '','400','', '500'});  
set(gca,'Layer','top');
set(gca, "YTick",[-1000 -900 -800 -700 -600 -500 -400 -300 -200 -100 0])
grid on;
hold off
box on

function mean_half = local_mean_half_ppb(scenario, ...
    pf_vec, frac_vec, fw_vec, d14apm_vec, ...
    bg_CH4, d14bg, d13bg, ethanebg, ...
    d14ff, d14apm_base, d14wet, d13ff, d13apm, d13wt, ethaneff, ...
    U14, U13, Ueth, ...
    Nmeas, ...
    FW_total, FF_fixed, WT_fixed, ...
    K, den_floor, max_resamples, Wsqrt, cond_thresh, lambda_ridge)

    d14vec_mean0 = [d14ff, d14apm_base, d14wet];
    d13vec_mean0 = [d13ff, d13apm,      d13wt];

    switch scenario
        case 'excess'  
            xvec = pf_vec;   yvec = fw_vec;
            fF_ex = 0.5;
            get_cell = @(PF,Y) deal(fF_ex*Y, (1-fF_ex)*Y, d14apm_base); % FF, WT, d14apm_this

        case 'frac'   
            xvec = pf_vec;   yvec = frac_vec;
            get_cell = @(PF,fF) deal(fF*FW_total, (1-fF)*FW_total, d14apm_base);

        case 'age'    
            xvec = pf_vec;   yvec = d14apm_vec;
            get_cell = @(PF,d14this) deal(FF_fixed, WT_fixed, d14this);

        otherwise, error('Unknown scenario');
    end

    nX=numel(xvec); nY=numel(yvec);
    hr_map = nan(nX,nY);

    for ii=1:nX
        PF = xvec(ii);
        for jj=1:nY
            [FF,WT,d14apm_this] = get_cell(PF, yvec(jj));
            CH4_true = bg_CH4 + FF + PF + WT;

            d14_obs = (d14bg*bg_CH4 + d14ff*FF + d14apm_this*PF + d14wet*WT)/CH4_true;
            d13_obs = (d13bg*bg_CH4 + d13ff*FF + d13apm*PF      + d13wt*WT)/CH4_true;
            eth_obs = ethanebg + 1000*ethaneff*FF;


            d14vec_mean = d14vec_mean0; d14vec_mean(2)=d14apm_this;
            A_mean = [1,1,1; d14vec_mean; d13vec_mean0; 1000*[ethaneff,0,0]];


            x_pf = nan(1,K);
            for k=1:K

                CH4_k = CH4_true + Nmeas(1)*randn;
                d14_k = d14_obs  + Nmeas(2)*randn;
                d13_k = d13_obs  + Nmeas(3)*randn;
                eth_k = eth_obs  + Nmeas(4)*randn;
                b_k = [ CH4_k - bg_CH4;
                        d14_k*CH4_k - d14bg*bg_CH4;
                        d13_k*CH4_k - d13bg*bg_CH4;
                        eth_k - ethanebg ];

                % endmember draws (only those that are "on"; truncated ±3σ; eth ≥0)
                d14v = d14vec_mean + [U14(1) U14(2) U14(3)].*randn(1,3);
                d13v = d13vec_mean0 + [U13(1) U13(2) U13(3)].*randn(1,3);
                d14v = arrayfun(@(x,m,s)max(min(x,m+3*s),m-3*s), d14v, d14vec_mean, [U14(1) U14(2) U14(3)]);
                d13v = arrayfun(@(x,m,s)max(min(x,m+3*s),m-3*s), d13v, d13vec_mean0, [U13(1) U13(2) U13(3)]);
                ethv1 = max( max(min(ethaneff + Ueth*randn, ethaneff+3*Ueth), ethaneff-3*Ueth), 0 );

                A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);
                if any(~isfinite(x_draw)) || any(abs(x_draw)>1e6)
                    tried=1;
                    while tried<max_resamples && (~all(isfinite(x_draw)) || any(abs(x_draw)>1e6))
                        d14v = d14vec_mean + [U14(1) U14(2) U14(3)].*randn(1,3);
                        d13v = d13vec_mean0 + [U13(1) U13(2) U13(3)].*randn(1,3);
                        d14v = arrayfun(@(x,m,s)max(min(x,m+3*s),m-3*s), d14v, d14vec_mean, [U14(1) U14(2) U14(3)]);
                        d13v = arrayfun(@(x,m,s)max(min(x,m+3*s),m-3*s), d13v, d13vec_mean0, [U13(1) U13(2) U13(3)]);
                        ethv1 = max( max(min(ethaneff + Ueth*randn, ethaneff+3*Ueth), ethaneff-3*Ueth), 0 );
                        A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                        x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);
                        tried=tried+1;
                    end
                end
                x_pf(k) = max(x_draw(2),0);
            end

            p16 = prctile(x_pf,16); p84 = prctile(x_pf,84);
            hr_map(ii,jj) = 0.5*(p84 - p16); % [ppb]
        end
    end

    mean_half = mean(hr_map(:),'omitnan');
end



%
function y = clip3(x,mu,sig)
    y = max(min(x, mu+3*sig), mu-3*sig);
end
function x = stable_wls(A,b,Wsqrt,cond_thresh,lambda_ridge)
    Aw = Wsqrt*A; bw = Wsqrt*b;
    c = rcond(Aw'*Aw);
    if ~(isfinite(c)) || c < cond_thresh
        x = (Aw'*Aw + lambda_ridge*eye(size(A,2))) \ (Aw'*bw);
    else
        x = Aw \ bw;
    end
end
function [half_ppb, half_pct] = run_frac_plane(FW_total, pf_vec, frac_vec, ...
        bg_CH4, d14bg, d13bg, ethanebg, ...
        d14ff, d14apm, d14wet, d13ff, d13apm, d13wt, ethaneff, ...
        unc_d14vec, unc_d13vec, unc_ethaneff, ...
        noise_ch4, noise_d14C, noise_d13C, noise_eth, ...
        K, den_floor, max_resamples, Wsqrt, cond_thresh, lambda_ridge)

    nPF   = numel(pf_vec);
    nFrac = numel(frac_vec);
    half_ppb = nan(nPF,nFrac);
    half_pct = nan(nPF,nFrac);

    for ii = 1:nPF
        PF = pf_vec(ii);
        for jj = 1:nFrac
            fF = frac_vec(jj);
            FF = fF * FW_total;
            WT = (1-fF) * FW_total;

            CH4_true = bg_CH4 + FF + PF + WT;

            d14_obs = (d14bg*bg_CH4 + d14ff*FF + d14apm*PF + d14wet*WT) / CH4_true;
            d13_obs = (d13bg*bg_CH4 + d13ff*FF + d13apm*PF + d13wt*WT) / CH4_true;
            eth_obs = ethanebg + 1000*ethaneff*FF;

            d14vec_mean = [d14ff, d14apm, d14wet];
            d13vec_mean = [d13ff, d13apm, d13wt];

            x_pf = nan(1,K);
            for k = 1:K
                CH4_k = CH4_true + noise_ch4*randn;
                d14_k = d14_obs  + noise_d14C*randn;
                d13_k = d13_obs  + noise_d13C*randn;
                eth_k = eth_obs  + noise_eth*randn;

                b_k = [ CH4_k - bg_CH4;
                        d14_k*CH4_k - d14bg*bg_CH4;
                        d13_k*CH4_k - d13bg*bg_CH4;
                        eth_k - ethanebg ];

                d14v = d14_min_vec + (d14_max_vec - d14_min_vec).*rand(1,3);
                d13v = d13_min_vec + (d13_max_vec - d13_min_vec).*rand(1,3);
                ethv1 = ethane_min + (ethane_max - ethane_min)*rand;

                A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);

                if any(~isfinite(x_draw)) || any(abs(x_draw) > 1e6)
                    tried = 1;
                    while tried < max_resamples && ...
                          (~all(isfinite(x_draw)) || any(abs(x_draw) > 1e6))
                        d14v = d14_min_vec + (d14_max_vec - d14_min_vec).*rand(1,3);
                        d13v = d13_min_vec + (d13_max_vec - d13_min_vec).*rand(1,3);
                        ethv1 = ethane_min + (ethane_max - ethane_min)*rand;
                        A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                        x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);
                        tried = tried + 1;
                    end
                end
                x_pf(k) = max(x_draw(2),0);
            end

            p16 = prctile(x_pf,16);
            p84 = prctile(x_pf,84);
            hr_ppb = 0.5*(p84 - p16);

            half_ppb(ii,jj) = hr_ppb;
            half_pct(ii,jj) = 100 * hr_ppb / max(PF, den_floor);
        end
    end

    half_pct(half_pct>100) = 100;
    half_pct(half_pct<0)   = 0;
end


function [half_ppb, half_pct] = run_age_plane(FF_fixed, WT_fixed, ...
        pf_vec, d14apm_vec, ...
        bg_CH4, d14bg, d13bg, ethanebg, ...
        d14ff, d14wet, d13ff, d13apm, d13wt, ethaneff, ...
        unc_d14vec, unc_d13vec, unc_ethaneff, ...
        noise_ch4, noise_d14C, noise_d13C, noise_eth, ...
        K, den_floor, max_resamples, Wsqrt, cond_thresh, lambda_ridge)

    nPF  = numel(pf_vec);
    nD14 = numel(d14apm_vec);
    half_ppb = nan(nPF,nD14);
    half_pct = nan(nPF,nD14);

    for ii = 1:nPF
        PF = pf_vec(ii);
        for jj = 1:nD14
            d14apm_this = d14apm_vec(jj);

            FF = FF_fixed;
            WT = WT_fixed;
            CH4_true = bg_CH4 + FF + PF + WT;

            d14_obs = (d14bg*bg_CH4 + d14ff*FF + d14apm_this*PF + d14wet*WT) / CH4_true;
            d13_obs = (d13bg*bg_CH4 + d13ff*FF + d13apm*PF       + d13wt*WT) / CH4_true;
            eth_obs = ethanebg + 1000*ethaneff*FF;

            d14vec_mean = [d14ff, d14apm_this, d14wet];
            d13vec_mean = [d13ff, d13apm,      d13wt];

            x_pf = nan(1,K);
            for k = 1:K
                CH4_k = CH4_true + noise_ch4*randn;
                d14_k = d14_obs  + noise_d14C*randn;
                d13_k = d13_obs  + noise_d13C*randn;
                eth_k = eth_obs  + noise_eth*randn;

                b_k = [ CH4_k - bg_CH4;
                        d14_k*CH4_k - d14bg*bg_CH4;
                        d13_k*CH4_k - d13bg*bg_CH4;
                        eth_k - ethanebg ];

                d14v = d14vec_mean + unc_d14vec.*randn(1,3);
                d13v = d13vec_mean + unc_d13vec.*randn(1,3);
                d14v = arrayfun(@(x,m,s) clip3(x,m,s), d14v, d14vec_mean, unc_d14vec);
                d13v = arrayfun(@(x,m,s) clip3(x,m,s), d13v, d13vec_mean, unc_d13vec);
                ethv1 = max( clip3(ethaneff + unc_ethaneff*randn, ethaneff, unc_ethaneff), 0 );

                A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);

                if any(~isfinite(x_draw)) || any(abs(x_draw) > 1e6)
                    tried = 1;
                    while tried < max_resamples && ...
                          (~all(isfinite(x_draw)) || any(abs(x_draw) > 1e6))
                        d14v = d14_min_vec + (d14_max_vec - d14_min_vec).*rand(1,3);
                        d13v = d13_min_vec + (d13_max_vec - d13_min_vec).*rand(1,3);
                        ethv1 = ethane_min + (ethane_max - ethane_min)*rand;
                        A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                        x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);
                        tried = tried + 1;
                    end
                end
                x_pf(k) = max(x_draw(2),0);
            end

            p16 = prctile(x_pf,16);
            p84 = prctile(x_pf,84);
            hr_ppb = 0.5*(p84 - p16);

            half_ppb(ii,jj) = hr_ppb;
            half_pct(ii,jj) = 100 * hr_ppb / max(PF, den_floor);
        end
    end

    half_pct(half_pct>100) = 100;
    half_pct(half_pct<0)   = 0;
end
function [half_ppb, half_pct] = run_one_fraction_total(fF, pf_vec, EX_vec, ...
        bg_CH4, d14bg, d13bg, ethanebg, ...
        d14vec_mean, d13vec_mean, ethaneff, ...
        d14_min_vec, d14_max_vec, d13_min_vec, d13_max_vec, ethane_min, ethane_max, ...
        noise_ch4, noise_d14C, noise_d13C, noise_eth, ...
        K, den_floor, max_resamples, Wsqrt, cond_thresh, lambda_ridge)

    nPF = numel(pf_vec); 
    nEX = numel(EX_vec);
    half_ppb = nan(nPF,nEX);
    half_pct = nan(nPF,nEX);

    for ii = 1:nPF
        PF = pf_vec(ii);
        for jj = 1:nEX
            EX = EX_vec(jj);      
            if EX < PF

                continue;
            end

            FW_total = EX - PF;
            FF = fF * FW_total;
            WT = (1 - fF) * FW_total;
            CH4_true = bg_CH4 + EX;

            d14_obs = (d14bg*bg_CH4 + d14vec_mean(1)*FF + ...
                       d14vec_mean(2)*PF + d14vec_mean(3)*WT) / CH4_true;
            d13_obs = (d13bg*bg_CH4 + d13vec_mean(1)*FF + ...
                       d13vec_mean(2)*PF + d13vec_mean(3)*WT) / CH4_true;
            eth_obs = ethanebg + 1000*ethaneff*FF;


            x_pf = nan(1,K);
            for k=1:K

                CH4_k = CH4_true + noise_ch4*randn;
                d14_k = d14_obs  + noise_d14C*randn;
                d13_k = d13_obs  + noise_d13C*randn;
                eth_k = eth_obs  + noise_eth*randn;

                b_k = [ CH4_k - bg_CH4;
                        d14_k*CH4_k - d14bg*bg_CH4;
                        d13_k*CH4_k - d13bg*bg_CH4;
                        eth_k - ethanebg ];


                d14v = d14_min_vec + (d14_max_vec - d14_min_vec).*rand(1,3);
                d13v = d13_min_vec + (d13_max_vec - d13_min_vec).*rand(1,3);
                ethv1 = ethane_min + (ethane_max - ethane_min)*rand;

                A_k  = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);

                if any(~isfinite(x_draw)) || any(abs(x_draw)>1e6)
                    tried = 1;
                    while tried < max_resamples && ...
                          (~all(isfinite(x_draw)) || any(abs(x_draw)>1e6))
                        d14v = d14_min_vec + (d14_max_vec - d14_min_vec).*rand(1,3);
                        d13v = d13_min_vec + (d13_max_vec - d13_min_vec).*rand(1,3);
                        ethv1 = ethane_min + (ethane_max - ethane_min)*rand;
                        A_k  = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                        x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);
                        tried = tried + 1;
                    end
                end

                x_pf(k) = max(x_draw(2),0);
            end

            p16   = prctile(x_pf,16);
            p84   = prctile(x_pf,84);
            hr_ppb = 0.5*(p84 - p16);
            hr_pct = 100 * hr_ppb / max(PF, den_floor);

            half_ppb(ii,jj) = hr_ppb;
            half_pct(ii,jj) = min(max(hr_pct,0),100);
        end
    end
end

