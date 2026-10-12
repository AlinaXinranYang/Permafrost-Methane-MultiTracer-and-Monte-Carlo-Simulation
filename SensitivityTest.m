%% Gaussian Monte Carlo contribution + tracer-dependence sensitivity
% Contribution test matched to the three sensitivity planes in the main STD simulation.
% For each uncertainty scenario, the code separately evaluates:
%   (1) PF x total excess CH4, with fossil fraction fixed at 20%
%   (2) PF x non-permafrost fossil fraction, with F+W fixed at 100 ppb
%   (3) PF x permafrost Delta14C, with FF = 20 ppb and WT = 80 ppb
% The mean MC 16-84% half-range is calculated for each plane, then the three
% plane means are averaged with equal weight. Source-signature uncertainty is
% Gaussian and truncated at +/-3 sigma, matching the main STD simulation.

%% ===== Editable parameters =====
dDff = -197; dDwt = -358; dDapm = -387;

bg_CH4 = 1922.6; d14bg = 350; d13bg = -47.6; dDbg = -99; ethanebg = 1.5435;

d14ff = -1000; d14apm = -799; d14wet = 37.7;
d13ff = -44.8; d13apm = -75.8; d13wt = -67.8;
ethaneff = 0.0652;

% Gaussian source-signature standard deviations
unc_d14ff = 0;   unc_d14apm = 125;  unc_d14wet = 38;
unc_d13ff = 5.1; unc_d13apm = 11.5; unc_d13wet = 8.1;
unc_ethaneff = 0.055;

% Measurement unc used in the main Gaussian simulation
noise_ch4 = 3; noise_d14C = 5; noise_d13C = 0.2; noise_eth = 0.07;

K = 500; rng(1);
max_resamples = 5;
cond_thresh = 1e-10;
lambda_ridge = 1e-6;

U_d14 = [unc_d14ff, unc_d14apm, unc_d14wet];   % [fossil PF wetland]
U_d13 = [unc_d13ff, unc_d13apm, unc_d13wet];
U_eth = unc_ethaneff;
N_meas = [noise_ch4, noise_d14C, noise_d13C, noise_eth];

ZERO14 = [0 0 0];
ZERO13 = [0 0 0];
ZEROmeas = [0 0 0 0];

w = [1/noise_ch4^2; 1/noise_d14C^2; 1/noise_d13C^2; 1/noise_eth^2];
Wsqrt = diag(sqrt(w));


pf_contrib_vec = (0:10:500)';

% Plane 1: PF x total excess CH4, fossil fraction fixed at 20%
EX_contrib_vec = (0:10:700)';
fF_total_fixed = 0.20;

% Plane 2: PF x fossil fraction, total non-permafrost CH4 fixed at 100 ppb
fF_contrib_vec = (0:0.02:1)';
FW_frac_fixed = 100;

% Plane 3: PF x permafrost Delta14C, FF=20 ppb and WT=80 ppb
d14_contrib_vec = (-1000:30:0)';
FF_age_fixed = 20;
WT_age_fixed = 80;


K_contrib = 200;

src_names = { ...
    'Meas unc CH_4','Meas unc \Delta^{14}C','Meas unc \delta^{13}C','Meas unc C_2H_6', ...
    '\Delta^{14}C_{FF}','\Delta^{14}C_{PM}','\Delta^{14}C_{WT}', ...
    '\delta^{13}C_{FF}','\delta^{13}C_{PM}','\delta^{13}C_{WT}', ...
    'C_2H_6:CH_4_{FF}'};

nSrc = numel(src_names);


run_three_planes = @(U14,U13,Ueth,Nm) local_mean_half_three_planes( ...
    pf_contrib_vec,EX_contrib_vec,fF_contrib_vec,d14_contrib_vec, ...
    fF_total_fixed,FW_frac_fixed,FF_age_fixed,WT_age_fixed, ...
    bg_CH4,d14bg,d13bg,ethanebg, ...
    d14ff,d14apm,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
    U14,U13,Ueth,Nm,K_contrib,max_resamples,Wsqrt,cond_thresh,lambda_ridge);

%% ===== TAKE-OFF-ONE uncertainty contribution test =====
% Start with ALL uncertainties ON, then remove one uncertainty at a time.
% Contribution = reduction in the equal-weight three-plane MC half-range metric.
rng(1);
total_metric = run_three_planes(U_d14,U_d13,U_eth,N_meas);
metric_without = nan(nSrc,1);

rng(1); metric_without(1)  = run_three_planes(U_d14,U_d13,U_eth,[0 N_meas(2) N_meas(3) N_meas(4)]);
rng(1); metric_without(2)  = run_three_planes(U_d14,U_d13,U_eth,[N_meas(1) 0 N_meas(3) N_meas(4)]);
rng(1); metric_without(3)  = run_three_planes(U_d14,U_d13,U_eth,[N_meas(1) N_meas(2) 0 N_meas(4)]);
rng(1); metric_without(4)  = run_three_planes(U_d14,U_d13,U_eth,[N_meas(1) N_meas(2) N_meas(3) 0]);

rng(1); metric_without(5)  = run_three_planes([0 U_d14(2) U_d14(3)],U_d13,U_eth,N_meas);
rng(1); metric_without(6)  = run_three_planes([U_d14(1) 0 U_d14(3)],U_d13,U_eth,N_meas);
rng(1); metric_without(7)  = run_three_planes([U_d14(1) U_d14(2) 0],U_d13,U_eth,N_meas);

rng(1); metric_without(8)  = run_three_planes(U_d14,[0 U_d13(2) U_d13(3)],U_eth,N_meas);
rng(1); metric_without(9)  = run_three_planes(U_d14,[U_d13(1) 0 U_d13(3)],U_eth,N_meas);
rng(1); metric_without(10) = run_three_planes(U_d14,[U_d13(1) U_d13(2) 0],U_eth,N_meas);

rng(1); metric_without(11) = run_three_planes(U_d14,U_d13,0,N_meas);

contrib_metric = total_metric - metric_without;
% Tiny negative values can occur from Monte Carlo sampling noise.
contrib_metric(contrib_metric < 0) = 0;
contrib_pct = 100 * contrib_metric ./ sum(contrib_metric,'omitnan');

fprintf('\n=== TAKE-OFF-ONE uncertainty sensitivity: average of three Figure-5 planes ===\n');
fprintf('All-on three-plane mean half-range = %.3f ppb\n',total_metric);
for i = 1:nSrc
    fprintf('%-32s %6.2f %%\n',src_names{i},contrib_pct(i));
end


noise_ch4_tracer    = 1.7;
noise_d14C_tracer   = 5;
noise_d13C_tracer   = 0.2;
noise_dD_tracer     = 1.7;
noise_ethane_tracer = 0.07;

pf_vec  = (0:10:500)';
FW_vec  = (0:10:200)';
fF_list = 0:0.2:1;

nPF = numel(pf_vec); nFW = numel(FW_vec); nF = numel(fF_list);
K_rep = 10;
N = nPF * nFW * nF * K_rep;

d14vec    = [d14ff, d14apm, d14wet];
d13vec    = [d13ff, d13apm, d13wt];
dDvec     = [dDff, dDapm, dDwt];
ethanevec = [ethaneff, 0, 0];

A_full = [1,1,1; d14vec; d13vec; dDvec; ethanevec];

err_all = nan(N,1); err_no_d14 = nan(N,1); err_no_d13 = nan(N,1);
err_no_dD = nan(N,1); err_no_eth = nan(N,1); ref_PF = nan(N,1);
idx = 0;

for iPF = 1:nPF
    PF = pf_vec(iPF);
    for iFW = 1:nFW
        FW_total = FW_vec(iFW);
        for iF = 1:nF
            fF = fF_list(iF);
            FF = fF * FW_total;
            WT = (1-fF) * FW_total;
            for k = 1:K_rep
                idx = idx + 1;
                ref_PF(idx) = PF;
                ch4_true = bg_CH4 + FF + PF + WT;
                if ch4_true <= 0, continue; end

                d14_true = (d14bg*bg_CH4 + d14ff*FF + d14apm*PF + d14wet*WT) / ch4_true;
                d13_true = (d13bg*bg_CH4 + d13ff*FF + d13apm*PF + d13wt*WT) / ch4_true;
                dD_true  = (dDbg*bg_CH4 + dDff*FF + dDapm*PF + dDwt*WT) / ch4_true;
                eth_true = ethanebg + ethaneff*FF;

                ch4_w = ch4_true + noise_ch4_tracer*randn;
                d14_w = d14_true + noise_d14C_tracer*randn;
                d13_w = d13_true + noise_d13C_tracer*randn;
                dD_w  = dD_true  + noise_dD_tracer*randn;
                eth_w = eth_true + noise_ethane_tracer*randn;

                b_full = [ch4_w-bg_CH4; ...
                          d14_w*ch4_w-d14bg*bg_CH4; ...
                          d13_w*ch4_w-d13bg*bg_CH4; ...
                          dD_w*ch4_w-dDbg*bg_CH4; ...
                          eth_w-ethanebg];

                sig_ch4 = noise_ch4_tracer;
                sig_d14 = noise_d14C_tracer*ch4_true;
                sig_d13 = noise_d13C_tracer*ch4_true;
                sig_dD  = noise_dD_tracer*ch4_true;
                sig_eth = noise_ethane_tracer;
                sigma = [sig_ch4; sig_d14; sig_d13; sig_dD; sig_eth];
                w_all = 1 ./ (sigma.^2);

                x_all = lscov(A_full,b_full,w_all);
                err_all(idx) = x_all(2)-PF;

                keep_rows = [1 3 4 5];
                x_no_d14 = lscov(A_full(keep_rows,:),b_full(keep_rows),w_all(keep_rows));
                err_no_d14(idx) = x_no_d14(2)-PF;

                keep_rows = [1 2 4 5];
                x_no_d13 = lscov(A_full(keep_rows,:),b_full(keep_rows),w_all(keep_rows));
                err_no_d13(idx) = x_no_d13(2)-PF;

                keep_rows = [1 2 3 5];
                x_no_dD = lscov(A_full(keep_rows,:),b_full(keep_rows),w_all(keep_rows));
                err_no_dD(idx) = x_no_dD(2)-PF;

                keep_rows = [1 2 3 4];
                x_no_eth = lscov(A_full(keep_rows,:),b_full(keep_rows),w_all(keep_rows));
                err_no_eth(idx) = x_no_eth(2)-PF;
            end
        end
    end
end

mask = isfinite(ref_PF) & isfinite(err_all);
ref_PF = ref_PF(mask); err_all = err_all(mask); err_no_d14 = err_no_d14(mask);
err_no_d13 = err_no_d13(mask); err_no_dD = err_no_dD(mask); err_no_eth = err_no_eth(mask);

abs_all = abs(err_all); abs_no_d14 = abs(err_no_d14); abs_no_d13 = abs(err_no_d13);
abs_no_dD = abs(err_no_dD); abs_no_eth = abs(err_no_eth);

MEANABS_all = mean(abs_all,'omitnan');
MEANABS_no_d14 = mean(abs_no_d14,'omitnan');
MEANABS_no_d13 = mean(abs_no_d13,'omitnan');
MEANABS_no_dD = mean(abs_no_dD,'omitnan');
MEANABS_no_eth = mean(abs_no_eth,'omitnan');

fprintf('\nMean absolute PF error (ppb), tracer-dependence test:\n');
fprintf('  All tracers : %.2f ppb\n',MEANABS_all);
fprintf('  No Delta14C : %.2f ppb\n',MEANABS_no_d14);
fprintf('  No delta13C : %.2f ppb\n',MEANABS_no_d13);
fprintf('  No deltaD   : %.2f ppb\n',MEANABS_no_dD);
fprintf('  No ethane   : %.2f ppb\n',MEANABS_no_eth);

%% ===== Two-tile figure  =====
figure('Color','w');
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')

clr_noise = [0.55 0.70 1.00; 0.35 0.55 0.90; 0.15 0.35 0.75; 0.00 0.20 0.55];
clr_d14   = [0.40 0.80 0.40; 0.25 0.60 0.25; 0.10 0.40 0.10];
clr_d13   = [1.00 0.70 0.40; 0.95 0.55 0.10; 0.80 0.35 0.00];
clr_eth   = [0.70 0.50 0.80];
bar_colors = [clr_noise; clr_d14; clr_d13; clr_eth];

% ----- TILE 1: stacked percentage contributions -----
ax = nexttile; hold(ax,'on');
[contrib_pct_sorted,sortIdx] = sort(contrib_pct,'ascend');
src_names_sorted = src_names(sortIdx);
bar_colors_sorted = bar_colors(sortIdx,:);

xpos = 1;
b = bar(ax,xpos,contrib_pct_sorted','stacked','BarWidth',0.65);
for i = 1:numel(b)
    set(b(i),'FaceColor',bar_colors_sorted(i,:),'EdgeColor','none');
end

legend_labels = cell(size(src_names_sorted));
for i = 1:numel(src_names_sorted)
    p = contrib_pct_sorted(i);
    if p == 0
        pct_str = '0%';
    elseif p < 0.1
        pct_str = sprintf('%.1g%%',p);
    else
        pct_str = sprintf('%.1f%%',p);
    end
    legend_labels{i} = sprintf('%s (%s)',src_names_sorted{i},pct_str);
end

lg1 = legend(ax,flip(b),flip(legend_labels),'Location','northeast');
set(lg1,'FontSize',22); lg1.Box = 'on';
ylabel(ax,'Contribution to total uncertainty (%)');
title(ax,'Uncertainty Sensitivity Test');
set(ax,'FontSize',26); grid(ax,'on'); ylim(ax,[0 100]); set(ax,'XTick',[]); xlim(ax,[0 2]);
hold(ax,'off');

% ----- TILE 2: tracer-dependence boxchart of absolute error -----
ax2 = nexttile; hold(ax2,'on');
vals_cell = { ...
    abs_all,      "All tracers"; ...
    abs_no_d14,   "No \Delta^{14}C"; ...
    abs_no_d13,   "No \delta^{13}C"; ...
    abs_no_dD,    "No \deltaD"; ...
    abs_no_eth,   "No C_2H_6"};

C = lines(5); h_mean = [];
for k = 1:5
    vals_k = vals_cell{k,1}; name_k = vals_cell{k,2};
    bc = boxchart(categorical(repmat(name_k,numel(vals_k),1)),vals_k, ...
        'MarkerStyle','none','BoxWidth',0.5,'BoxFaceAlpha',0.0);
    edgeCol = C(k,:); faceCol = 0.7*C(k,:) + 0.3;
    bc.BoxEdgeColor = edgeCol; bc.WhiskerLineColor = edgeCol;
    bc.LineWidth = 1.8; bc.BoxFaceColor = faceCol; bc.BoxFaceAlpha = 0.5;
    mean_k = mean(vals_k,'omitnan');
    h = plot(categorical(name_k),mean_k,'x','MarkerSize',20,'LineWidth',2.5,'Color','k');
    if isempty(h_mean), h_mean = h; end
    text(categorical(name_k),mean_k+2,sprintf('%.1f',mean_k), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',25,'Color','k','FontWeight','bold');
end

ylabel(ax2,{'Permafrost CH_4','Difference from Truth (ppb)'});
title(ax2,'Tracer-Dependence Sensitivity Test');
set(ax2,'FontSize',26); grid(ax2,'on'); ylim(ax2,[0 80]);
lg2 = legend(ax2,h_mean,'Mean'); set(lg2,'FontSize',22); lg2.Box = 'on';
hold(ax2,'off');

%% Bias diagnostics for tile 2
BIAS_all = mean(err_all,'omitnan'); BIAS_no_d14 = mean(err_no_d14,'omitnan');
BIAS_no_d13 = mean(err_no_d13,'omitnan'); BIAS_no_dD = mean(err_no_dD,'omitnan');
BIAS_no_eth = mean(err_no_eth,'omitnan');
fprintf('\nMean bias in PF estimate (ppb), tracer-dependence test:\n');
fprintf('  All tracers : %+6.2f\n',BIAS_all);
fprintf('  No Delta14C : %+6.2f\n',BIAS_no_d14);
fprintf('  No delta13C : %+6.2f\n',BIAS_no_d13);
fprintf('  No deltaD   : %+6.2f\n',BIAS_no_dD);
fprintf('  No ethane   : %+6.2f\n',BIAS_no_eth);



%% 
function mean_three = local_mean_half_three_planes( ...
        pf_vec,EX_vec,fF_vec,d14apm_vec, ...
        fF_total_fixed,FW_frac_fixed,FF_age_fixed,WT_age_fixed, ...
        bg_CH4,d14bg,d13bg,ethanebg, ...
        d14ff,d14apm_base,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
        U14,U13,Ueth,Nmeas,K,max_resamples,Wsqrt,cond_thresh,lambda_ridge)

    hr_total = nan(numel(pf_vec),numel(EX_vec));
    for ii = 1:numel(pf_vec)
        PF = pf_vec(ii);
        for jj = 1:numel(EX_vec)
            EX = EX_vec(jj);
            if EX < PF
                continue
            end
            FW_total = EX - PF;
            FF = fF_total_fixed * FW_total;
            WT = (1-fF_total_fixed) * FW_total;
            hr_total(ii,jj) = local_cell_half_range(PF,FF,WT,d14apm_base, ...
                bg_CH4,d14bg,d13bg,ethanebg,d14ff,d14wet, ...
                d13ff,d13apm,d13wt,ethaneff,U14,U13,Ueth,Nmeas,K, ...
                max_resamples,Wsqrt,cond_thresh,lambda_ridge);
        end
    end

    hr_frac = nan(numel(pf_vec),numel(fF_vec));
    for ii = 1:numel(pf_vec)
        PF = pf_vec(ii);
        for jj = 1:numel(fF_vec)
            fF = fF_vec(jj);
            FF = fF * FW_frac_fixed;
            WT = (1-fF) * FW_frac_fixed;
            hr_frac(ii,jj) = local_cell_half_range(PF,FF,WT,d14apm_base, ...
                bg_CH4,d14bg,d13bg,ethanebg,d14ff,d14wet, ...
                d13ff,d13apm,d13wt,ethaneff,U14,U13,Ueth,Nmeas,K, ...
                max_resamples,Wsqrt,cond_thresh,lambda_ridge);
        end
    end

    hr_age = nan(numel(pf_vec),numel(d14apm_vec));
    for ii = 1:numel(pf_vec)
        PF = pf_vec(ii);
        for jj = 1:numel(d14apm_vec)
            d14apm_this = d14apm_vec(jj);
            hr_age(ii,jj) = local_cell_half_range(PF,FF_age_fixed,WT_age_fixed,d14apm_this, ...
                bg_CH4,d14bg,d13bg,ethanebg,d14ff,d14wet, ...
                d13ff,d13apm,d13wt,ethaneff,U14,U13,Ueth,Nmeas,K, ...
                max_resamples,Wsqrt,cond_thresh,lambda_ridge);
        end
    end

    mean_total = mean(hr_total(:),'omitnan');
    mean_frac  = mean(hr_frac(:),'omitnan');
    mean_age   = mean(hr_age(:),'omitnan');

    mean_three = mean([mean_total,mean_frac,mean_age],'omitnan');
end

function hr = local_cell_half_range(PF,FF,WT,d14apm_this, ...
        bg_CH4,d14bg,d13bg,ethanebg,d14ff,d14wet, ...
        d13ff,d13apm,d13wt,ethaneff,U14,U13,Ueth,Nmeas,K, ...
        max_resamples,Wsqrt,cond_thresh,lambda_ridge)

    CH4_true = bg_CH4 + PF + FF + WT;
    d14_obs = (d14bg*bg_CH4 + d14ff*FF + d14apm_this*PF + d14wet*WT) / CH4_true;
    d13_obs = (d13bg*bg_CH4 + d13ff*FF + d13apm*PF + d13wt*WT) / CH4_true;
    eth_obs = ethanebg + 1000*ethaneff*FF;

    d14_mean = [d14ff,d14apm_this,d14wet];
    d13_mean = [d13ff,d13apm,d13wt];
    x_pf = nan(1,K);

    for r = 1:K
        CH4_k = CH4_true + Nmeas(1)*randn;
        d14_k = d14_obs  + Nmeas(2)*randn;
        d13_k = d13_obs  + Nmeas(3)*randn;
        eth_k = eth_obs  + Nmeas(4)*randn;

        b_k = [CH4_k-bg_CH4; ...
               d14_k*CH4_k-d14bg*bg_CH4; ...
               d13_k*CH4_k-d13bg*bg_CH4; ...
               eth_k-ethanebg];

        d14v = d14_mean + U14.*randn(1,3);
        d13v = d13_mean + U13.*randn(1,3);
        d14v = arrayfun(@(x,m,s)clip3(x,m,s),d14v,d14_mean,U14);
        d13v = arrayfun(@(x,m,s)clip3(x,m,s),d13v,d13_mean,U13);
        ethv1 = max(clip3(ethaneff + Ueth*randn,ethaneff,Ueth),0);

        A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
        x_draw = stable_wls(A_k,b_k,Wsqrt,cond_thresh,lambda_ridge);

        if any(~isfinite(x_draw)) || any(abs(x_draw)>1e6)
            tried = 1;
            while tried < max_resamples && (~all(isfinite(x_draw)) || any(abs(x_draw)>1e6))
                d14v = d14_mean + U14.*randn(1,3);
                d13v = d13_mean + U13.*randn(1,3);
                d14v = arrayfun(@(x,m,s)clip3(x,m,s),d14v,d14_mean,U14);
                d13v = arrayfun(@(x,m,s)clip3(x,m,s),d13v,d13_mean,U13);
                ethv1 = max(clip3(ethaneff + Ueth*randn,ethaneff,Ueth),0);
                A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                x_draw = stable_wls(A_k,b_k,Wsqrt,cond_thresh,lambda_ridge);
                tried = tried + 1;
            end
        end
        x_pf(r) = max(x_draw(2),0);
    end

    p16 = prctile(x_pf,16);
    p84 = prctile(x_pf,84);
    hr = 0.5*(p84-p16);
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

function y = clip3(x,mu,sig)
    y = max(min(x,mu+3*sig),mu-3*sig);
end
