
% ===== Editable parameters =====
pf_max = 500; pf_step = 10; frac_step = 0.02; FW_total = 200;
bg_CH4 = 1950; d14bg = 350; d13bg = -47.6; dDbg = -99; ethanebg = 1.194; % ppb
d14ff = -1000; d14apm = -799; d14wet = 37.7;
d13ff = -44.8;  d13apm = -75.8; d13wt = -67.8;
ethaneff = 0.0652;                     % ethane per ppb CH4 fossil
dDff = -197;dDwt = -358;dDapm = -387;
%new unc trial 
unc_d14ff=0; unc_d14apm=99; unc_d14wet=3;
unc_d13ff=3; unc_d13apm=11.5; unc_d13wet=12.8;
unc_ethaneff=0.035;

%unc_d14ff=0; unc_d14apm=39; unc_d14wet=3;
%unc_d13ff=3; unc_d13apm=11.5; unc_d13wet=12.8;
%unc_ethaneff=0.005;

%noise
noise_ch4=1.7; noise_d14C=5; noise_d13C=0.2; noise_eth=0.07;
K = 500; rng(1);
den_floor = 10;            % ppb floor for % metrics near PF≈0
max_resamples = 5;         % retries if a draw is unstable
cond_thresh = 1e-10;       % conditioning guard for WLS
lambda_ridge = 1e-6;   



% ---- 1. Source names ----
src_names = { ...
  'Noise CH_4','Noise \Delta^{14}C','Noise \delta^{13}C','Noise ethane', ...
  '\Delta^{14}C_{fossil}','\Delta^{14}C_{permafrost}','\Delta^{14}C_{wetland}', ...
  '\delta^{13}C_{fossil}','\delta^{13}C_{permafrost}','\delta^{13}C_{wetland}', ...
  'C_2H_6:CH_4_{fossil}'};

% ---- 2. Cache baseline uncertainties (so we can toggle them) ----
U_d14  = [unc_d14ff, unc_d14apm, unc_d14wet];  % [ff PF wet]
U_d13  = [unc_d13ff, unc_d13apm, unc_d13wet];
U_eth  =  unc_ethaneff;
N_meas = [noise_ch4, noise_d14C, noise_d13C, noise_eth];

% measurement weights (assumed fixed; we only toggle the actual noise added)
w = [ 1/noise_ch4^2; 1/noise_d14C^2; 1/noise_d13C^2; 1/noise_eth^2 ];
Wsqrt = diag(sqrt(w));

% ---- 3. Grids: PF 0–500 ppb, F+W 0–200 ppb, fossil fraction 0–1 ----
pf_max   = 500;   pf_step   = 10;
FW_min   = 0;     FW_max    = 200;   FW_step = 20;
fF_min   = 0;     fF_max    = 1.0;   fF_step = 0.2;  % 0,0.2,...,1

pf_vec   = (0:pf_step:pf_max)';      % permafrost CH4 (ppb)
FW_vec   = (FW_min:FW_step:FW_max)'; % total F+W CH4 (ppb)
fF_vec   = (fF_min:fF_step:fF_max)'; % fossil fraction of F+W

% ---- 4. Convenience zero-vectors for toggling ----
ZERO14  = [0 0 0];
ZERO13  = [0 0 0];
ZEROmeas = [0 0 0 0];

% ---- 5. Helper: run one global scenario and return mean half-range [ppb] ----
run_mean_half_ppb_global = @(U14, U13, Ueth, Nmeas) ...
    local_mean_half_ppb_global( ...
        pf_vec, FW_vec, fF_vec, ...
        bg_CH4, d14bg, d13bg, ethanebg, ...
        d14ff, d14apm, d14wet, ...
        d13ff, d13apm, d13wt, ...
        ethaneff, ...
        U14, U13, Ueth, ...
        Nmeas, ...
        K, den_floor, max_resamples, ...
        Wsqrt, cond_thresh, lambda_ridge);

% ---- 6. Compute global contributions (one bar) ----
nSrc        = numel(src_names);
contrib_ppb = zeros(nSrc,1);

% "ALL" scenario: everything on
total_ppb = run_mean_half_ppb_global(U_d14, U_d13, U_eth, N_meas);

% Noise-only scenarios
contrib_ppb(1)  = run_mean_half_ppb_global(ZERO14, ZERO13, 0, [N_meas(1) 0 0 0]); % Noise CH4
contrib_ppb(2)  = run_mean_half_ppb_global(ZERO14, ZERO13, 0, [0 N_meas(2) 0 0]); % Noise d14
contrib_ppb(3)  = run_mean_half_ppb_global(ZERO14, ZERO13, 0, [0 0 N_meas(3) 0]); % Noise d13
contrib_ppb(4)  = run_mean_half_ppb_global(ZERO14, ZERO13, 0, [0 0 0 N_meas(4)]); % Noise ethane

% Δ14C endmembers
contrib_ppb(5)  = run_mean_half_ppb_global([U_d14(1) 0 0], ZERO13, 0, ZEROmeas);   % Δ14C_ff
contrib_ppb(6)  = run_mean_half_ppb_global([0 U_d14(2) 0], ZERO13, 0, ZEROmeas);   % Δ14C_PF
contrib_ppb(7)  = run_mean_half_ppb_global([0 0 U_d14(3)], ZERO13, 0, ZEROmeas);   % Δ14C_wet

% δ13C endmembers
contrib_ppb(8)  = run_mean_half_ppb_global(ZERO14, [U_d13(1) 0 0], 0, ZEROmeas);   % δ13C_ff
contrib_ppb(9)  = run_mean_half_ppb_global(ZERO14, [0 U_d13(2) 0], 0, ZEROmeas);   % δ13C_PF
contrib_ppb(10) = run_mean_half_ppb_global(ZERO14, [0 0 U_d13(3)], 0, ZEROmeas);   % δ13C_wet

% Ethane ratio
contrib_ppb(11) = run_mean_half_ppb_global(ZERO14, ZERO13, U_eth, ZEROmeas);       % ethane ratio

% ==== Leave-one-out RMS test over PF–(F+W)–fossil-fraction space ====
% Assumes the following are already defined in your script:
%   bg_CH4, d14bg, d13bg, dDbg, ethanebg   (background)
%   d14ff, d14apm, d14wet
%   d13ff, d13apm, d13wt
%   dDff,  dDapm,  dDwt
%   ethaneff
%
% and measurement noise:
noise_ch4    = 5;
noise_d14C   = 5;
noise_d13C   = 0.2;
noise_dD     = 1.7;
noise_ethane = 0.07;

% ----- Grids -----
pf_vec   = (0:10:500)';   % permafrost CH4 (ppb)
FW_vec   = (0:10:200)';   % total F+W excess CH4 (ppb)
fF_list  = 0:0.2:1;       % fossil fraction of (F+W): 0,0.2,...,1.0

nPF   = numel(pf_vec);
nFW   = numel(FW_vec);
nF    = numel(fF_list);
K_rep = 10;               % number of noise realisations per grid point

% Total number of samples
N = nPF * nFW * nF * K_rep;

% ----- Endmember vectors for WLS -----
d14vec    = [d14ff,  d14apm,  d14wet];
d13vec    = [d13ff,  d13apm,  d13wt ];
dDvec     = [dDff,   dDapm,   dDwt  ];
ethanevec = [ethaneff, 0, 0];   % ethane [ppb] = ethaneff * FF

% A_full rows:
% 1: CH4 excess  (CH4exc)
% 2: Δ14C excess (d14exc)
% 3: δ13C excess (d13exc)
% 4: δD excess   (dDexc)
% 5: ethane exc  (ethane_exc)
A_full = [ 1, 1, 1;
           d14vec;
           d13vec;
           dDvec;
           ethanevec ];

% Storage for errors
err_all    = nan(N,1);
err_no_d14 = nan(N,1);
err_no_d13 = nan(N,1);
err_no_dD  = nan(N,1);
err_no_eth = nan(N,1);
ref_PF     = nan(N,1);

idx = 0;

for iPF = 1:nPF
    PF = pf_vec(iPF);
    for iFW = 1:nFW
        FW_total = FW_vec(iFW);

        for iF = 1:nF
            fF = fF_list(iF);
            FF = fF * FW_total;
            WT = (1 - fF) * FW_total;

            for k = 1:K_rep
                idx = idx + 1;
                ref_PF(idx) = PF;

                % ----- True (noiseless) total CH4 & isotopes -----
                ch4_true = bg_CH4 + FF + PF + WT;

                if ch4_true <= 0
                    % shouldn't happen, but guard anyway
                    err_all(idx)    = NaN;
                    err_no_d14(idx) = NaN;
                    err_no_d13(idx) = NaN;
                    err_no_dD(idx)  = NaN;
                    err_no_eth(idx) = NaN;
                    continue;
                end

                % true mixed isotopic signatures
                d14_true = (d14bg*bg_CH4 + d14ff*FF  + d14apm*PF  + d14wet*WT) / ch4_true;
                d13_true = (d13bg*bg_CH4 + d13ff*FF  + d13apm*PF  + d13wt*WT ) / ch4_true;
                dD_true  = (dDbg *bg_CH4 + dDff *FF  + dDapm*PF   + dDwt*WT  ) / ch4_true;

                % true ethane
                eth_true = ethanebg + ethaneff * FF;

                % ----- Add measurement noise (as in your main code) -----
                ch4_w = ch4_true + noise_ch4    * randn;
                d14_w = d14_true + noise_d14C   * randn;
                d13_w = d13_true + noise_d13C   * randn;
                dD_w  = dD_true  + noise_dD     * randn;
                eth_w = eth_true + noise_ethane * randn;

                % ----- Convert to "excess" quantities -----
                ch4exc    = ch4_w - bg_CH4;
                d14exc    = d14_w * ch4_w - d14bg * bg_CH4;
                d13exc    = d13_w * ch4_w - d13bg * bg_CH4;
                dDexc     = dD_w  * ch4_w - dDbg  * bg_CH4;
                ethaneexc = eth_w - ethanebg;

                b_full = [ch4exc;
                          d14exc;
                          d13exc;
                          dDexc;
                          ethaneexc];

                % ----- Row-wise sigmas for weighting -----
                % Approximate σ for each equation
                sig_ch4 = noise_ch4;
                sig_d14 = noise_d14C * ch4_true;
                sig_d13 = noise_d13C * ch4_true;
                sig_dD  = noise_dD   * ch4_true;
                sig_eth = noise_ethane;

                sigma = [sig_ch4; sig_d14; sig_d13; sig_dD; sig_eth];
                w_all = 1 ./ (sigma.^2);

                % ========== 1) All tracers ==========
                x_all = lscov(A_full, b_full, w_all);
                err_all(idx) = x_all(2) - PF;

                % ========== Leave-one-out cases ==========
                % Row indices: 1=CH4, 2=Δ14C, 3=δ13C, 4=δD, 5=ethane
                % Always keep CH4 row (1), drop only the tracer row indicated.

                % Drop Δ14C (row 2)
                keep_rows = [1 3 4 5];
                x_no_d14 = lscov(A_full(keep_rows,:), b_full(keep_rows), w_all(keep_rows));
                err_no_d14(idx) = x_no_d14(2) - PF;

                % Drop δ13C (row 3)
                keep_rows = [1 2 4 5];
                x_no_d13 = lscov(A_full(keep_rows,:), b_full(keep_rows), w_all(keep_rows));
                err_no_d13(idx) = x_no_d13(2) - PF;

                % Drop δD (row 4)
                keep_rows = [1 2 3 5];
                x_no_dD = lscov(A_full(keep_rows,:), b_full(keep_rows), w_all(keep_rows));
                err_no_dD(idx) = x_no_dD(2) - PF;

                % Drop ethane (row 5)
                keep_rows = [1 2 3 4];
                x_no_eth = lscov(A_full(keep_rows,:), b_full(keep_rows), w_all(keep_rows));
                err_no_eth(idx) = x_no_eth(2) - PF;
            end
        end
    end
end

% ----- Clean NaNs (if any) -----
mask = isfinite(ref_PF) & isfinite(err_all);
ref_PF     = ref_PF(mask);
err_all    = err_all(mask);
err_no_d14 = err_no_d14(mask);
err_no_d13 = err_no_d13(mask);
err_no_dD  = err_no_dD(mask);
err_no_eth = err_no_eth(mask);

% Use absolute errors for boxchart
abs_all    = abs(err_all);
abs_no_d14 = abs(err_no_d14);
abs_no_d13 = abs(err_no_d13);
abs_no_dD  = abs(err_no_dD);
abs_no_eth = abs(err_no_eth);

% ----- RMS diagnostics (ppb) -----
RMS_all    = sqrt(mean(err_all.^2));
RMS_no_d14 = sqrt(mean(err_no_d14.^2));
RMS_no_d13 = sqrt(mean(err_no_d13.^2));
RMS_no_dD  = sqrt(mean(err_no_dD.^2));
RMS_no_eth = sqrt(mean(err_no_eth.^2));

fprintf('\nRMS PF error (ppb) over PF=0–500, F+W=0–200, fF=0–1 (with noise):\n');
fprintf('  All tracers      : %.2f ppb\n', RMS_all);
fprintf('  No Δ14C          : %.2f ppb\n', RMS_no_d14);
fprintf('  No δ13C          : %.2f ppb\n', RMS_no_d13);
fprintf('  No δD            : %.2f ppb\n', RMS_no_dD);
fprintf('  No ethane        : %.2f ppb\n\n', RMS_no_eth);

% ----- Boxchart of |PF_est - PF_true| -----
vals = [abs_all; abs_no_d14; abs_no_d13; abs_no_dD; abs_no_eth];
grp  = [repmat("All tracers", numel(abs_all),    1); ...
        repmat("No Δ^{14}C",  numel(abs_no_d14), 1); ...
        repmat("No δ^{13}C",  numel(abs_no_d13), 1); ...
        repmat("No δD",       numel(abs_no_dD),  1); ...
        repmat("No ethane",   numel(abs_no_eth), 1)];
%% updated plotting in percentage 

figure('Color','w');
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')

%===================== TILE 1: stacked % contributions (numeric axis, flipped legend) =====================
ax = nexttile;
hold(ax,'on');

% --- Convert ppb contributions -> percentage of SUM(individual contributions) ---
tot_ppb     = sum(contrib_ppb, 'omitnan');   % scalar
contrib_pct = 100 * (contrib_ppb ./ tot_ppb); % [nSources x 1]

% --- Put stacked bar at numeric x = 1 (so xlim works) ---
xpos = 1;
b = bar(ax, xpos, contrib_pct', 'stacked', 'BarWidth', 0.65);

ylabel(ax,'Contribution to total uncertainty (%)');
title(ax,'Uncertainty Sensitivity Test');
set(ax,'FontSize',26);   % bigger axis ticks (even though we hide XTicks)

% --- Colours (unchanged) ---
clr_noise = [0.55 0.70 1.00; 0.35 0.55 0.90; 0.15 0.35 0.75; 0.00 0.20 0.55];
clr_d14   = [0.40 0.80 0.40; 0.25 0.60 0.25; 0.10 0.40 0.10];
clr_d13   = [1.00 0.70 0.40; 0.95 0.55 0.10; 0.80 0.35 0.00];
clr_eth   = [0.70 0.50 0.80];
bar_colors = [clr_noise; clr_d14; clr_d13; clr_eth];

for i = 1:numel(b)
    set(b(i),'FaceColor',bar_colors(i,:),'EdgeColor','none');
end

src_names = { ...
  'Meas Unc CH_4','Meas Unc \Delta^{14}C','Meas Unc \delta^{13}C','Meas Unc ethane', ...
  '\Delta^{14}C_{ff}','\Delta^{14}C_{pm}','\Delta^{14}C_{wt}', ...
  '\delta^{13}C_{ff}','\delta^{13}C_{pm}','\delta^{13}C_{wt}', ...
  'C_2H_6:CH_4_{ff}'};

% --- Flip legend order to match visual top-to-bottom of the stack ---
src_names_flip = flip(src_names);
b_flip         = flip(b);

lg1 = legend(ax, b_flip, src_names_flip, 'Location','northeast');
set(lg1,'FontSize',22);     % bigger legend text
lg1.Box = 'on';

grid(ax,'on');
ylim(ax,[0 100]);

% --- Remove x ticks entirely ---
set(ax,'XTick',[]);

% --- Control horizontal space manually (room for legend on right) ---
xlim(ax,[0 2]);

hold(ax,'off');


%===================== TILE 2: absolute PF error in ppb (bigger text) =====================
ax2 = nexttile;
hold(ax2,'on');

vals_cell = {
    abs_all,      "All tracers"
    abs_no_d14,   "No \Delta^{14}C"
    abs_no_d13,   "No \delta^{13}C"
    abs_no_dD,    "No \deltaD"
    abs_no_eth,   "No C_2H_6"
};

C = lines(5);
h_rms = [];

for k = 1:5
    vals_k = vals_cell{k,1};
    name_k = vals_cell{k,2};

    bc = boxchart( ...
        categorical(repmat(name_k, numel(vals_k), 1)), ...
        vals_k, ...
        'MarkerStyle','none', ...
        'BoxWidth',0.5, ...
        'BoxFaceAlpha',0.0 );

    edgeCol = C(k,:);
    faceCol = 0.7*C(k,:) + 0.3;

    bc.BoxEdgeColor     = edgeCol;
    bc.WhiskerLineColor = edgeCol;
    bc.LineWidth        = 1.8;
    bc.BoxFaceColor     = faceCol;
    bc.BoxFaceAlpha     = 0.5;

    rms_k = sqrt(mean(vals_k.^2,'omitnan'));

    h = plot(categorical(name_k), rms_k, 'x', ...
        'MarkerSize',20, ...
        'LineWidth',2.5, ...
        'Color','k');

    if isempty(h_rms)
        h_rms = h;
    end

    text(categorical(name_k), rms_k + 2, sprintf('%.1f', rms_k), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',25, ...
        'Color','k', ...
        'FontWeight','bold');
end

ylabel(ax2, {'Permafrost CH_4','Difference from Truth (ppb)'});
title(ax2,'Tracer-Dependence Sensitivity Test');
set(ax2,'FontSize',26);
grid(ax2,'on');

ylim(ax2,[0 80]);  % keep your original scale

lg2 = legend(ax2, h_rms, 'RMS');
set(lg2,'FontSize',22);
lg2.Box = 'on';

hold(ax2,'off');


%% plotting in ppb


figure('Color','w');
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')

% ---- FIRST TILE: stacked contributions ----
ax = nexttile;    % <-- use the tile axes directly
hold(ax,'on');


cats = categorical({'Global (PF 0–500, F+W 0–200, fossil 0–100%)'});
b = bar(ax, cats, contrib_ppb','stacked');
ylabel(ax,'Permafrost CH_4 uncertainty (ppb)');
title(ax,'Uncertainty Sensitivity Test');
set(ax,'FontSize',20);

% colour groups (same style as before)
clr_noise = [0.55 0.70 1.00; 0.35 0.55 0.90; 0.15 0.35 0.75; 0.00 0.20 0.55]; % blues
clr_d14   = [0.40 0.80 0.40; 0.25 0.60 0.25; 0.10 0.40 0.10];                 % greens
clr_d13   = [1.00 0.70 0.40; 0.95 0.55 0.10; 0.80 0.35 0.00];                 % oranges
clr_eth   = [0.70 0.50 0.80];                                                 % purple
bar_colors = [clr_noise; clr_d14; clr_d13; clr_eth];

for i = 1:numel(b)
    set(b(i),'FaceColor',bar_colors(i,:),'EdgeColor','none');
end
src_names = { ...
  'Meas Unc CH_4','Meas Unc \Delta^{14}C','Meas Unc \delta^{13}C','Meas Unc ethane', ...
  '\Delta^{14}C_{ff}','\Delta^{14}C_{pm}','\Delta^{14}C_{wt}', ...
  '\delta^{13}C_{ff}','\delta^{13}C_{pm}','\delta^{13}C_{wt}', ...
  'C_2H_6:CH_4_{ff}'};
legend(ax, src_names, 'Location','southoutside');

% Overlay yellow dot for "ALL" total
%plot(ax, cats, total_ppb, 'o', 'MarkerFaceColor',[1 1 0], ...
     %'MarkerEdgeColor','k', 'MarkerSize',7);

grid(ax,'on');
hold(ax,'off');
set(gca,'xtick',[])

nexttile
vals_cell = {
    abs_all,      "All tracers"
    abs_no_d14,   "No \Delta^{14}C"
    abs_no_d13,   "No \delta^{13}C"
    abs_no_dD,    "No \deltaD"
    abs_no_eth,   "No C_2H_6"
};

C = lines(5);   % colours for boxes

hold on

% Store handles for legend
h_rms = [];

for k = 1:5
    vals_k = vals_cell{k,1};
    name_k = vals_cell{k,2};

    % ------ Draw box -------
    bc = boxchart( ...
        categorical(repmat(name_k, numel(vals_k), 1)), ...
        vals_k, ...
        'MarkerStyle','none', ...
        'BoxWidth',0.5, ...
        'BoxFaceAlpha',0.0 );

    edgeCol = C(k,:);
    faceCol = 0.7*C(k,:) + 0.3;

    bc.BoxEdgeColor     = edgeCol;
    bc.WhiskerLineColor = edgeCol;
    bc.LineWidth        = 1.5;
    bc.BoxFaceColor     = faceCol;
    bc.BoxFaceAlpha     = 0.5;

    % ------- RMS value -------
    rms_k = sqrt(mean(vals_k.^2,'omitnan'));

    % add black cross
    h = plot(categorical(name_k), rms_k, 'x', ...
        'MarkerSize',12, ...
        'LineWidth',2, ...
        'Color','k');

    % Save one handle for legend
    if isempty(h_rms)
        h_rms = h;
    end

    % label
    text(categorical(name_k), rms_k + 3, sprintf('%.1f', rms_k), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',20, ...
        'Color','k', ...
        'FontWeight','bold');
end

ylabel({'Permafrost CH_4', 'Difference from Truth (ppb)'});
title('Tracer-Dependence Sensitivity Test');
set(gca,'FontSize',20);
grid on;
ylim([0 80]);

% --- Legend ---
legend(h_rms, 'RMS', 'FontSize',20);

hold off;




%% ===== Local helper: global mean half-range [ppb] =====
function mean_half = local_mean_half_ppb_global( ...
        pf_vec, FW_vec, fF_vec, ...
        bg_CH4, d14bg, d13bg, ethanebg, ...
        d14ff, d14apm, d14wet, ...
        d13ff, d13apm, d13wt, ...
        ethaneff, ...
        U14, U13, Ueth, ...
        Nmeas, ...
        K, den_floor, max_resamples, ...
        Wsqrt, cond_thresh, lambda_ridge)

    nPF  = numel(pf_vec);
    nFW  = numel(FW_vec);
    nFFr = numel(fF_vec);

    half_ppb_all = nan(nPF, nFW, nFFr);

    for ii = 1:nPF
        PF = pf_vec(ii);

        for jj = 1:nFW
            FW_total = FW_vec(jj);

            for kk = 1:nFFr
                fF = fF_vec(kk);        % fossil fraction of F+W
                FF = fF * FW_total;
                WT = (1 - fF) * FW_total;

                CH4_true = bg_CH4 + PF + FF + WT;

                % Mean (noiseless) mixed signals (using fixed d14apm)
                d14_obs = (d14bg*bg_CH4 + d14ff*FF + d14apm*PF + d14wet*WT) / CH4_true;
                d13_obs = (d13bg*bg_CH4 + d13ff*FF + d13apm*PF + d13wt*WT) / CH4_true;
                eth_obs = ethanebg + 1000*ethaneff*FF;

                x_pf = nan(1,K);

                for k = 1:K
                    % --- measurement draws (using scenario Nmeas) ---
                    CH4_k = CH4_true + Nmeas(1)*randn;
                    d14_k = d14_obs  + Nmeas(2)*randn;
                    d13_k = d13_obs  + Nmeas(3)*randn;
                    eth_k = eth_obs  + Nmeas(4)*randn;

                    b_k = [ CH4_k - bg_CH4;
                            d14_k*CH4_k - d14bg*bg_CH4;
                            d13_k*CH4_k - d13bg*bg_CH4;
                            eth_k - ethanebg ];

                    % --- endmember draws (truncated ±3σ; ethane ≥ 0) ---
                    d14v = [d14ff, d14apm, d14wet] + U14.*randn(1,3);
                    d13v = [d13ff, d13apm, d13wt]  + U13.*randn(1,3);

                    d14v = arrayfun(@(x,m,s)clip3(x,m,s), d14v, ...
                                    [d14ff, d14apm, d14wet], U14);
                    d13v = arrayfun(@(x,m,s)clip3(x,m,s), d13v, ...
                                    [d13ff, d13apm, d13wt],  U13);

                    ethv1 = max( clip3(ethaneff + Ueth*randn, ethaneff, Ueth), 0 );

                    A_k = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];

                    x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);

                    if any(~isfinite(x_draw)) || any(abs(x_draw)>1e6)
                        tried=1;
                        while tried<max_resamples && ...
                              (~all(isfinite(x_draw)) || any(abs(x_draw)>1e6))
                            d14v = [d14ff, d14apm, d14wet] + U14.*randn(1,3);
                            d13v = [d13ff, d13apm, d13wt]  + U13.*randn(1,3);

                            d14v = arrayfun(@(x,m,s)clip3(x,m,s), d14v, ...
                                            [d14ff, d14apm, d14wet], U14);
                            d13v = arrayfun(@(x,m,s)clip3(x,m,s), d13v, ...
                                            [d13ff, d13apm, d13wt],  U13);

                            ethv1 = max( clip3(ethaneff + Ueth*randn, ethaneff, Ueth), 0 );
                            A_k   = [1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
                            x_draw = stable_wls(A_k, b_k, Wsqrt, cond_thresh, lambda_ridge);
                            tried = tried + 1;
                        end
                    end

                    x_pf(k) = max(x_draw(2),0);  % PF estimate
                end

                p16 = prctile(x_pf,16);
                p84 = prctile(x_pf,84);
                hr_ppb = 0.5*(p84 - p16);

                half_ppb_all(ii,jj,kk) = hr_ppb;
            end
        end
    end

    mean_half = mean(half_ppb_all(:), 'omitnan');
end

%% ===== stable weighted least squares =====
function x = stable_wls(A,b,Wsqrt,cond_thresh,lambda_ridge)
    Aw = Wsqrt*A;
    bw = Wsqrt*b;
    c = rcond(Aw'*Aw);
    if ~(isfinite(c)) || c < cond_thresh
        x = (Aw'*Aw + lambda_ridge*eye(size(A,2))) \ (Aw'*bw);
    else
        x = Aw \ bw;
    end
end

%% ===== clip helper =====
function y = clip3(x,mu,sig)
    y = max(min(x, mu+3*sig), mu-3*sig);
end

