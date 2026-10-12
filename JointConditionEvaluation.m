%% Joint detection boundary with fixed non-permafrost CH4 = 100 ppb
% Gaussian/SD and Uniform source-signature uncertainty in the same script.
%
% Scenario dimensions:
%   PF_true x TOTAL fossil fraction of total excess x PF d14C
%
% Fixed condition:
%   FF + WT = 100 ppb
%
% Total fossil fraction is defined as:
%   fTot = FF / (PF + FF + WT) = FF / (PF + 100)
% Therefore, for each PF and fTot:
%   FF = fTot*(PF + 100)
%   WT = 100 - FF
%
% A combination is physically valid only when FF <= 100 (WT >= 0).
% Thus the feasible total-fossil-fraction range decreases as PF increases.
%
% Detection edges:
%   percentage uncertainty = 100 * [(P84-P16)/2] / PF_true
% Thresholds are 15% and 30%.
%
% Outputs:
%   JointDetection_FixedNonPF100_Gaussian.csv
%   JointDetection_FixedNonPF100_Uniform.csv
%   JointDetection_FixedNonPF100_results.mat
%
% Figures:
%   1) explicit total-fossil-fraction slices (0,10,20,30,50%) vs PF d14C
%   2) four stepped-colour heatmaps (Gaussian/Uniform x 15/30%)
clear; clc;
rng(1);
%% ===================== PARAMETERS =====================
bg_CH4   = 1922.6;
d14bg    = 350;
d13bg    = -47.6;
ethanebg = 1.5435;
d14ff  = -1000;
d14wet = 37.7;
d13ff  = -44.8;
d13apm = -75.8;
d13wt  = -67.8;
ethaneff = 0.0652;
% Gaussian source-signature SDs: [fossil, permafrost, wetland]
unc_d14 = [0, 125, 38];
unc_d13 = [5.1, 11.5, 8.1];
unc_eth = 0.055
unc_eth = 0.055/2;
% Uniform source-signature ranges
% PF d14C is handled around the scenario-centre value below.
d14wet_min = -28.1;
d14wet_max = 103.5;
d13_min = [-51, -81.6, -83];
d13_max = [-35, -61.2, -55];
ethane_min = 0.009
ethane_max = 0.20
ethane_min = 0.009+(0.2-0.009)/2;
ethane_max = 0.20-(0.2-0.009)/2;
pf_d14_uniform_halfwidth = sqrt(3)*125;
% Measurement noise stays Gaussian in both treatments
noise_ch4  = 3;
noise_d14C = 5;
noise_d13C = 0.2;
noise_eth  = 0.07;
K = 200;                         % change to 500 for final run
max_resamples = 5;
cond_thresh = 1e-10;
lambda_ridge = 1e-6;
w = [1/noise_ch4^2; 1/noise_d14C^2; 1/noise_d13C^2; 1/noise_eth^2];
Wsqrt = diag(sqrt(w));
%% ===================== JOINT GRID =====================
FW_fixed = 100;                  % fixed non-PF CH4 = FF + WT, ppb
fTot_vec = sort(unique([(0:0.02:0.50), 0.05, 0.25]))';  % 2% grid + explicit 5% and 25% slices
D14_vec  = (-1000:50:0)';        % PF d14C scenario centre, permil
PF_coarse = (10:10:500)';
PF_refine_step = 2;
thresholds_pct = [15 30];
fprintf('\nFixed non-PF joint detection search\n');
fprintf('FF+WT fixed at %.0f ppb\n',FW_fixed);
fprintf('Total fossil-fraction points = %d, d14C points = %d\n',numel(fTot_vec),numel(D14_vec));
fprintf('MC draws per tested PF = %d\n',K);
%% ===================== GAUSSIAN =====================
fprintf('\nRunning Gaussian / SD boundary...\n');
T_gauss = run_boundary_fixed_nonpf('gaussian', ...
    FW_fixed, fTot_vec, D14_vec, PF_coarse, PF_refine_step, thresholds_pct, ...
    bg_CH4,d14bg,d13bg,ethanebg, ...
    d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
    unc_d14,unc_d13,unc_eth, ...
    d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
    pf_d14_uniform_halfwidth, ...
    noise_ch4,noise_d14C,noise_d13C,noise_eth, ...
    K,max_resamples,Wsqrt,cond_thresh,lambda_ridge);
writetable(T_gauss,'JointDetection_FixedNonPF100_Gaussian.csv');
%% ===================== UNIFORM =====================
fprintf('\nRunning Uniform boundary...\n');
T_uniform = run_boundary_fixed_nonpf('uniform', ...
    FW_fixed, fTot_vec, D14_vec, PF_coarse, PF_refine_step, thresholds_pct, ...
    bg_CH4,d14bg,d13bg,ethanebg, ...
    d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
    unc_d14,unc_d13,unc_eth, ...
    d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
    pf_d14_uniform_halfwidth, ...
    noise_ch4,noise_d14C,noise_d13C,noise_eth, ...
    K,max_resamples,Wsqrt,cond_thresh,lambda_ridge);
writetable(T_uniform,'JointDetection_FixedNonPF100_Uniform.csv');
%% ===================== PRINT SUMMARY =====================
print_summary(T_gauss,'Gaussian');
print_summary(T_uniform,'Uniform');
%% ===================== FIGURE 1: EXPLICIT FOSSIL-FRACTION SLICES =====================
ff_show = [0 5 10 20 25 30 50];
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
plot_age_slices(T_gauss,T_uniform,ff_show,15);
xlabel('Permafrost \Delta^{14}C (‰)');
ylabel('Minimum detectable permafrost CH_4 (ppb)');
title('15% uncertainty edge');
grid on; box on; set(gca,'FontSize',16);
nexttile; hold on;
plot_age_slices(T_gauss,T_uniform,ff_show,30);
xlabel('Permafrost \Delta^{14}C (‰)');
ylabel('Minimum detectable permafrost CH_4 (ppb)');
title('30% uncertainty edge');
grid on; box on; set(gca,'FontSize',16);
sgtitle('Fixed non-permafrost CH_4 = 100 ppb; total fossil fraction');
%% ===================== FIGURE 2: STEPPED HEATMAPS =====================
figure('Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
plot_boundary_heatmap_step(T_gauss,15,'Gaussian, 15%');
plot_boundary_heatmap_step(T_uniform,15,'Uniform, 15%');
plot_boundary_heatmap_step(T_gauss,30,'Gaussian, 30%');
plot_boundary_heatmap_step(T_uniform,30,'Uniform, 30%');
sgtitle('Minimum permafrost CH_4 required; non-permafrost CH_4 fixed at 100 ppb');

%% ===================== FIGURE 3: FIXED-PF BOUNDARY SLICES =====================
% For selected true permafrost CH4 contributions, show the maximum total
% fossil fraction that still satisfies the detection threshold as a function
% of permafrost Delta14C. Gaussian = solid, Uniform = dashed.

pf_levels_show = [50 60 70 80 90 100 125 150 175 200 250];

figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot_pf_level_boundaries(T_gauss,T_uniform,pf_levels_show,15);
xlabel('Permafrost \Delta^{14}C (‰)');
ylabel('Maximum total fossil fraction (%)');
title('15% uncertainty boundary');
grid on; box on; set(gca,'FontSize',16);

nexttile; hold on;
plot_pf_level_boundaries(T_gauss,T_uniform,pf_levels_show,30);
xlabel('Permafrost \Delta^{14}C (‰)');
ylabel('Maximum total fossil fraction (%)');
title('30% uncertainty boundary');
grid on; box on; set(gca,'FontSize',16);

sgtitle('Fixed permafrost CH_4 contribution; non-permafrost CH_4 = 100 ppb');

%% ===================== SAVE =====================
save('JointDetection_FixedNonPF100_results.mat', ...
    'T_gauss','T_uniform','FW_fixed','fTot_vec','D14_vec', ...
    'PF_coarse','PF_refine_step','thresholds_pct','K');
fprintf('\nFinished.\n');
%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================
function T = run_boundary_fixed_nonpf(mode, ...
    FW_fixed, fTot_vec, D14_vec, PF_coarse, PF_refine_step, thresholds_pct, ...
    bg_CH4,d14bg,d13bg,ethanebg, ...
    d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
    unc_d14,unc_d13,unc_eth, ...
    d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
    pf_d14_uniform_halfwidth, ...
    noise_ch4,noise_d14C,noise_d13C,noise_eth, ...
    K,max_resamples,Wsqrt,cond_thresh,lambda_ridge)
    nF = numel(fTot_vec);
    nD = numel(D14_vec);
    N = nF*nD;
    fTot_out = nan(N,1);
    D14_out = nan(N,1);
    PF15_out = nan(N,1); U15_out = nan(N,1);
    PF30_out = nan(N,1); U30_out = nan(N,1);
    FF15_out = nan(N,1); WT15_out = nan(N,1); EX15_out = nan(N,1);
    FF30_out = nan(N,1); WT30_out = nan(N,1); EX30_out = nan(N,1);
    row = 0;
    for iF=1:nF
        fTot = fTot_vec(iF);
        for iD=1:nD
            d14_this = D14_vec(iD);
            row=row+1;
            fTot_out(row)=100*fTot;
            D14_out(row)=d14_this;
            rng(100000 + iF*1000 + iD);
            Ucoarse = nan(size(PF_coarse));
            validPF = false(size(PF_coarse));
            for ip=1:numel(PF_coarse)
                PF=PF_coarse(ip);
                EX=PF+FW_fixed;
                FF=fTot*EX;
                WT=FW_fixed-FF;
                if WT < -1e-12
                    continue
                end
                validPF(ip)=true;
                Ucoarse(ip)=mc_half_pct_fixed_nonpf(PF,fTot,d14_this,mode,FW_fixed, ...
                    bg_CH4,d14bg,d13bg,ethanebg, ...
                    d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
                    unc_d14,unc_d13,unc_eth, ...
                    d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
                    pf_d14_uniform_halfwidth, ...
                    noise_ch4,noise_d14C,noise_d13C,noise_eth, ...
                    K,max_resamples,Wsqrt,cond_thresh,lambda_ridge);
            end
            idx15=find(validPF & Ucoarse<=thresholds_pct(1),1,'first');
            idx30=find(validPF & Ucoarse<=thresholds_pct(2),1,'first');
            if ~isempty(idx15)
                [PF15_out(row),U15_out(row)] = refine_fixed_nonpf(thresholds_pct(1),idx15, ...
                    PF_coarse,Ucoarse,PF_refine_step,fTot,d14_this,mode,FW_fixed, ...
                    bg_CH4,d14bg,d13bg,ethanebg, ...
                    d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
                    unc_d14,unc_d13,unc_eth, ...
                    d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
                    pf_d14_uniform_halfwidth, ...
                    noise_ch4,noise_d14C,noise_d13C,noise_eth, ...
                    K,max_resamples,Wsqrt,cond_thresh,lambda_ridge);
                EX15_out(row)=PF15_out(row)+FW_fixed;
                FF15_out(row)=fTot*EX15_out(row);
                WT15_out(row)=FW_fixed-FF15_out(row);
            end
            if ~isempty(idx30)
                [PF30_out(row),U30_out(row)] = refine_fixed_nonpf(thresholds_pct(2),idx30, ...
                    PF_coarse,Ucoarse,PF_refine_step,fTot,d14_this,mode,FW_fixed, ...
                    bg_CH4,d14bg,d13bg,ethanebg, ...
                    d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
                    unc_d14,unc_d13,unc_eth, ...
                    d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
                    pf_d14_uniform_halfwidth, ...
                    noise_ch4,noise_d14C,noise_d13C,noise_eth, ...
                    K,max_resamples,Wsqrt,cond_thresh,lambda_ridge);
                EX30_out(row)=PF30_out(row)+FW_fixed;
                FF30_out(row)=fTot*EX30_out(row);
                WT30_out(row)=FW_fixed-FF30_out(row);
            end
            if mod(row,50)==0 || row==N
                fprintf('%s: %d / %d combinations complete\n',mode,row,N);
            end
        end
    end
    T=table(fTot_out,D14_out,PF15_out,U15_out,FF15_out,WT15_out,EX15_out, ...
        PF30_out,U30_out,FF30_out,WT30_out,EX30_out, ...
        'VariableNames',{'TotalFossilFraction_pct','PF_d14C_permille', ...
        'PFmin_15pct_ppb','AchievedUnc15_pct','Fossil15_ppb','Wetland15_ppb','TotalExcess15_ppb', ...
        'PFmin_30pct_ppb','AchievedUnc30_pct','Fossil30_ppb','Wetland30_ppb','TotalExcess30_ppb'});
end
function [PFcross,Ucross] = refine_fixed_nonpf(target,idx,PF_coarse,Ucoarse,step, ...
    fTot,d14_this,mode,FW_fixed, ...
    bg_CH4,d14bg,d13bg,ethanebg, ...
    d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
    unc_d14,unc_d13,unc_eth, ...
    d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
    pf_d14_uniform_halfwidth, ...
    noise_ch4,noise_d14C,noise_d13C,noise_eth, ...
    K,max_resamples,Wsqrt,cond_thresh,lambda_ridge)
    if idx==1
        lowPF=max(step,PF_coarse(1)-20);
    else
        lowPF=PF_coarse(idx-1)+step;
    end
    highPF=PF_coarse(idx);
    refine_vec=(lowPF:step:highPF)';
    PFcross=highPF;
    Ucross=Ucoarse(idx);
    for j=1:numel(refine_vec)
        PF=refine_vec(j);
        EX=PF+FW_fixed;
        FF=fTot*EX;
        WT=FW_fixed-FF;
        if WT < -1e-12
            continue
        end
        U=mc_half_pct_fixed_nonpf(PF,fTot,d14_this,mode,FW_fixed, ...
            bg_CH4,d14bg,d13bg,ethanebg, ...
            d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
            unc_d14,unc_d13,unc_eth, ...
            d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
            pf_d14_uniform_halfwidth, ...
            noise_ch4,noise_d14C,noise_d13C,noise_eth, ...
            K,max_resamples,Wsqrt,cond_thresh,lambda_ridge);
        if U<=target
            PFcross=PF; Ucross=U; return
        end
    end
end
function half_pct = mc_half_pct_fixed_nonpf(PF,fTot,d14_this,mode,FW_fixed, ...
    bg_CH4,d14bg,d13bg,ethanebg, ...
    d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
    unc_d14,unc_d13,unc_eth, ...
    d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
    pf_d14_uniform_halfwidth, ...
    noise_ch4,noise_d14C,noise_d13C,noise_eth, ...
    K,max_resamples,Wsqrt,cond_thresh,lambda_ridge)
    EX=PF+FW_fixed;
    FF=fTot*EX;
    WT=FW_fixed-FF;
    if PF<=0 || WT < -1e-12
        half_pct=NaN; return
    end
    WT=max(WT,0);
    CH4_true=bg_CH4+EX;
    d14_obs=(d14bg*bg_CH4+d14ff*FF+d14_this*PF+d14wet*WT)/CH4_true;
    d13_obs=(d13bg*bg_CH4+d13ff*FF+d13apm*PF+d13wt*WT)/CH4_true;
    eth_obs=ethanebg+1000*ethaneff*FF;
    x_pf=nan(K,1);
    for k=1:K
        CH4_k=CH4_true+noise_ch4*randn;
        d14_k=d14_obs+noise_d14C*randn;
        d13_k=d13_obs+noise_d13C*randn;
        eth_k=eth_obs+noise_eth*randn;
        b_k=[CH4_k-bg_CH4; ...
             d14_k*CH4_k-d14bg*bg_CH4; ...
             d13_k*CH4_k-d13bg*bg_CH4; ...
             eth_k-ethanebg];
        [d14v,d13v,ethv1]=draw_endmembers(mode,d14_this, ...
            d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
            unc_d14,unc_d13,unc_eth, ...
            d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
            pf_d14_uniform_halfwidth);
        A_k=[1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
        x_draw=stable_wls(A_k,b_k,Wsqrt,cond_thresh,lambda_ridge);
        tried=1;
        while tried<max_resamples && (any(~isfinite(x_draw)) || any(abs(x_draw)>1e6))
            [d14v,d13v,ethv1]=draw_endmembers(mode,d14_this, ...
                d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
                unc_d14,unc_d13,unc_eth, ...
                d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
                pf_d14_uniform_halfwidth);
            A_k=[1,1,1; d14v; d13v; 1000*[ethv1,0,0]];
            x_draw=stable_wls(A_k,b_k,Wsqrt,cond_thresh,lambda_ridge);
            tried=tried+1;
        end
        if all(isfinite(x_draw)) && all(abs(x_draw)<=1e6)
            x_pf(k)=max(x_draw(2),0);
        end
    end
    p16=prctile(x_pf,16);
    p84=prctile(x_pf,84);
    half_pct=100*(0.5*(p84-p16))/PF;
end
function [d14v,d13v,ethv1]=draw_endmembers(mode,d14_this, ...
    d14ff,d14wet,d13ff,d13apm,d13wt,ethaneff, ...
    unc_d14,unc_d13,unc_eth, ...
    d14wet_min,d14wet_max,d13_min,d13_max,ethane_min,ethane_max, ...
    pf_d14_uniform_halfwidth)
    switch lower(mode)
        case 'gaussian'
            mu14=[d14ff,d14_this,d14wet];
            mu13=[d13ff,d13apm,d13wt];
            d14v=mu14+unc_d14.*randn(1,3);
            d13v=mu13+unc_d13.*randn(1,3);
            d14v=arrayfun(@(x,m,s)clip3(x,m,s),d14v,mu14,unc_d14);
            d13v=arrayfun(@(x,m,s)clip3(x,m,s),d13v,mu13,unc_d13);
            ethv1=max(clip3(ethaneff+unc_eth*randn,ethaneff,unc_eth),0);
        case 'uniform'
            d14pf=(d14_this-pf_d14_uniform_halfwidth)+2*pf_d14_uniform_halfwidth*rand;
            d14wet=d14wet_min+(d14wet_max-d14wet_min)*rand;
            d14v=[-1000,d14pf,d14wet];
            d13v=d13_min+(d13_max-d13_min).*rand(1,3);
            ethv1=ethane_min+(ethane_max-ethane_min)*rand;
        otherwise
            error('Unknown mode');
    end
end
function plot_age_slices(Tg,Tu,ff_show,target)
    sets={Tg,Tu}; names={'Gaussian','Uniform'}; styles={'-','--'};
    cols=lines(numel(ff_show));
    for is=1:2
        B=sets{is};
        for jf=1:numel(ff_show)
            f=ff_show(jf);
            sub=B(abs(B.TotalFossilFraction_pct-f)<1e-9,:);
            if isempty(sub), continue; end
            if target==15, y=sub.PFmin_15pct_ppb; else, y=sub.PFmin_30pct_ppb; end
            x=sub.PF_d14C_permille;
            good=isfinite(x)&isfinite(y);
            if ~any(good), continue; end
            [x,ord]=sort(x(good)); yy=y(good); yy=yy(ord);
            plot(x,yy,'LineWidth',2,'LineStyle',styles{is},'Color',cols(jf,:), ...
                'DisplayName',sprintf('%s, %g%% total fossil',names{is},f));
        end
    end
    legend('Location','best');
end
function plot_boundary_heatmap_step(B,target,titleStr)
    ax=nexttile;
    fvals=unique(B.TotalFossilFraction_pct);
    dvals=unique(B.PF_d14C_permille);
    M=nan(numel(dvals),numel(fvals));
    if target==15, vv=B.PFmin_15pct_ppb; else, vv=B.PFmin_30pct_ppb; end
    for i=1:height(B)
        ir=find(dvals==B.PF_d14C_permille(i),1);
        ic=find(fvals==B.TotalFossilFraction_pct(i),1);
        M(ir,ic)=vv(i);
    end
    edges=[0 75 125 175 225 275 325 375 425 475 inf];
    labels={'50','100','150','200','250','300','350','400','450','>500 / none'};
    Mplot=M; Mplot(~isfinite(Mplot))=inf;
    idx=discretize(Mplot,edges);
    basepal=[0.000 0.400 0.400; ...
             0.220 0.664 0.700; ...
             0.600 0.940 0.980; ...
             0.900 1.000 1.000; ...
             1.000 0.793 0.600; ...
             1.000 0.560 0.200; ...
             0.600 0.250 0.000];
    cmap9=interp1(linspace(0,1,size(basepal,1)),basepal,linspace(0,1,9),'linear');
    cmap=[cmap9;0.75 0.75 0.75];
    imagesc(ax,fvals,dvals,idx); axis(ax,'xy');
    xlabel(ax,'Total fossil fraction of excess CH_4 (%)');
    ylabel(ax,'Permafrost \Delta^{14}C (‰)');
    title(ax,titleStr);
    colormap(ax,cmap); clim(ax,[1 10]); set(ax,'FontSize',14);
    cb=colorbar(ax); cb.Ticks=1:10; cb.TickLabels=labels;
    ylabel(cb,'Minimum PF CH_4 (ppb)');
end
function plot_max_fossil_boundary(Tg,Tu,PF_levels,target)
    sets = {Tg,Tu};
    names = {'Gaussian','Uniform'};
    styles = {'-','--'};
    cols = lines(numel(PF_levels));
    major = ismember(PF_levels,[50 100 150 200 250]);
    for is = 1:2
        B = sets{is};
        dvals = unique(B.PF_d14C_permille);
        for ip = 1:numel(PF_levels)
            pfAvail = PF_levels(ip);
            maxFF = nan(size(dvals));
            for id = 1:numel(dvals)
                sub = B(B.PF_d14C_permille==dvals(id),:);
                if target==15
                    req = sub.PFmin_15pct_ppb;
                else
                    req = sub.PFmin_30pct_ppb;
                end
                ok = isfinite(req) & req <= pfAvail;
                if any(ok)
                    maxFF(id) = max(sub.TotalFossilFraction_pct(ok));
                end
            end
            lw = 1.2;
            if major(ip), lw = 2.4; end
            plot(dvals,maxFF,'LineStyle',styles{is},'Color',cols(ip,:), ...
                'LineWidth',lw,'DisplayName',sprintf('%s, %g ppb',names{is},pfAvail));
        end
    end
    ylim([0 50]);
    legend('Location','eastoutside','FontSize',9);
end

function plot_pf_level_boundaries(Tg,Tu,pf_levels_show,target)
    sets = {Tg,Tu};
    names = {'Gaussian','Uniform'};
    styles = {'-','--'};
    cols = lines(numel(pf_levels_show));

    for is = 1:2
        B = sets{is};

        if target == 15
            pfmin = B.PFmin_15pct_ppb;
        else
            pfmin = B.PFmin_30pct_ppb;
        end

        dvals = unique(B.PF_d14C_permille);

        for ip = 1:numel(pf_levels_show)
            PFlevel = pf_levels_show(ip);
            xout = nan(numel(dvals),1);
            yout = nan(numel(dvals),1);

            for id = 1:numel(dvals)
                d = dvals(id);
                sub = B(B.PF_d14C_permille == d, :);
                if isempty(sub), continue; end

                if target == 15
                    pf_req = sub.PFmin_15pct_ppb;
                else
                    pf_req = sub.PFmin_30pct_ppb;
                end

                ok = isfinite(pf_req) & (pf_req <= PFlevel);
                if ~any(ok), continue; end

                xout(id) = d;
                yout(id) = max(sub.TotalFossilFraction_pct(ok));
            end

            good = isfinite(xout) & isfinite(yout);
            if ~any(good), continue; end

            lw = 1.4;
            if ismember(PFlevel,[50 100 150 200 250])
                lw = 2.6;
            end

            plot(xout(good),yout(good), ...
                'LineWidth',lw, ...
                'LineStyle',styles{is}, ...
                'Color',cols(ip,:), ...
                'DisplayName',sprintf('%s, %g ppb PF',names{is},PFlevel));
        end
    end

    legend('Location','best');
end

function print_summary(B,name)
    fprintf('\n=== %s fixed-nonPF boundary ===\n',name);
    for target=[15 30]
        v=B.(sprintf('PFmin_%dpct_ppb',target));
        ok=isfinite(v);
        if ~any(ok)
            fprintf('%d%%: no threshold found.\n',target); continue
        end
        ids=find(ok);
        [mn,imn]=min(v(ok)); [mx,imx]=max(v(ok));
        r1=ids(imn); r2=ids(imx);
        fprintf('%d%% edge: min PF %.1f ppb at total fossil %.0f%%, d14C %.0f permil\n', ...
            target,mn,B.TotalFossilFraction_pct(r1),B.PF_d14C_permille(r1));
        fprintf('          max PF %.1f ppb at total fossil %.0f%%, d14C %.0f permil\n', ...
            mx,B.TotalFossilFraction_pct(r2),B.PF_d14C_permille(r2));
        fprintf('          thresholds found for %.1f%% of tested fossil-age combinations\n',100*mean(ok));
    end
end
function y=clip3(x,mu,sig)
    if sig==0, y=mu; else, y=max(min(x,mu+3*sig),mu-3*sig); end
end
function x=stable_wls(A,b,Wsqrt,cond_thresh,lambda_ridge)
    Aw=Wsqrt*A; bw=Wsqrt*b;
    c=rcond(Aw'*Aw);
    if ~(isfinite(c)) || c<cond_thresh
        x=(Aw'*Aw+lambda_ridge*eye(size(A,2)))\(Aw'*bw);
    else
        x=Aw\bw;
    end
end