%% RUN_CALIBRATION_GRID_REPO
%
% Full calibration-grid run; ensemble count is set below (default 20).
%
% Calibration quantities:
%   1. HPV prevalence age 25-34  (group 3)
%   2. HPV prevalence age 35-44  (group 4)
%   3. Lesion prevalence age 25-49
%   4. Cervical cancer ASR
%
% Saves:
%   - results/calibration_allcomb_trajectories_ens20.mat

% Required model:
%   ABmodelZambia_paper_calibration_only_repo3.m
clearvars;
clc;
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
outputDir = fullfile(repoRoot, 'results');
if ~isfolder(outputDir)
    mkdir(outputDir);
end
%% SETTINGS
ens      = 20;
time     = 81;
burntime = 50;
%% PARAMETER GRIDS
muGrid   = 0.011:0.001:0.014;
PcumGrid = 0.29:0.01:0.33;
a2Grid   = 0.67:0.03:0.74;
a3Grid   = 0.05;
a4Grid   = 0.10;
%% CALIBRATION TARGETS
targetHPV       = [0, 0.35, 0.30, 0.20, 0.15, 0.15]*1.15;
targetLesion    = 35284 / 545581;
targetCancerASR = 65.5;
ageLabels = {'0-14','15-24','25-34','35-44','45-54','55+'};
%% MODEL
calibFcn = @ABmodelZambia_paper_calibration_only_repo3;
if exist(func2str(calibFcn),'file') ~= 2
    error('Cannot find %s.m on the MATLAB path.',func2str(calibFcn));
end
fprintf('Using calibration function: %s\n',func2str(calibFcn));
fprintf('Ensembles per parameter combination: %d\n',ens);
%% FULL PARAMETER GRID
[A2,A3,A4,M,P] = ndgrid(a2Grid,a3Grid,a4Grid,muGrid,PcumGrid);
A2 = A2(:);
A3 = A3(:);
A4 = A4(:);
M  = M(:);
P  = P(:);
nJobs = numel(A2);
fprintf('Number of parameter combinations: %d\n',nJobs);
fprintf('Total stochastic simulations: %d\n',nJobs*ens);
%% PREALLOCATE CALIBRATION OUTPUT
a2_out   = zeros(nJobs,1);
a3_out   = zeros(nJobs,1);
a4_out   = zeros(nJobs,1);
mu_out   = zeros(nJobs,1);
Pcum_out = zeros(nJobs,1);
loss_out = NaN(nJobs,1);
Lhpv_out = NaN(nJobs,1);
Lles_out = NaN(nJobs,1);
Lcan_out = NaN(nJobs,1);
lesion_out = NaN(nJobs,1);
cancer_out = NaN(nJobs,1);
hpv_out    = NaN(nJobs,6);
%% PREALLOCATE ANNUAL TRAJECTORIES
% Ensemble means for each parameter combination.
cancerTrajAll = NaN(time,nJobs);
lesionTrajAll = NaN(time,nJobs);
hpvTrajAll    = NaN(time,6,nJobs);
%% PARALLEL POOL
useParallel = license('test','Distrib_Computing_Toolbox');
if useParallel
    pool = gcp('nocreate');
    if isempty(pool)
        try
            parpool('local');
        catch ME
            warning('Could not start parallel pool: %s',ME.message);
            useParallel = false;
        end
    end
end
%% RUN CALIBRATION GRID
if useParallel
    fprintf('Running grid in parallel...\n');
    parfor j = 1:nJobs
        pars_local = struct( ...
            'a2',   A2(j), ...
            'a3',   A3(j), ...
            'a4',   A4(j), ...
            'mu',   M(j), ...
            'Pcum', P(j));
        [~,~,~,~,~,calib_j] = calibFcn(pars_local,time,ens);
        [L,L_hpv,L_les,L_can] = calibLoss( ...
            calib_j,targetHPV,targetLesion,targetCancerASR);
        a2_out(j)   = A2(j);
        a3_out(j)   = A3(j);
        a4_out(j)   = A4(j);
        mu_out(j)   = M(j);
        Pcum_out(j) = P(j);
        loss_out(j) = L;
        Lhpv_out(j) = L_hpv;
        Lles_out(j) = L_les;
        Lcan_out(j) = L_can;
        lesion_out(j) = calib_j.lesion;
        cancer_out(j) = calib_j.cancer;
        hpv_out(j,:)  = calib_j.hpv;
        cancerTrajAll(:,j) = mean(calib_j.cancerASRByYear,2,'omitnan');
        lesionTrajAll(:,j) = mean(calib_j.lesionByYear,2,'omitnan');
        hpvTrajAll(:,:,j)  = mean(calib_j.hpvByYear,3,'omitnan');
    end
else
    fprintf('Running grid sequentially...\n');
    for j = 1:nJobs
        fprintf('Combination %d/%d\n',j,nJobs);
        pars_local = struct( ...
            'a2',   A2(j), ...
            'a3',   A3(j), ...
            'a4',   A4(j), ...
            'mu',   M(j), ...
            'Pcum', P(j));
        [~,~,~,~,~,calib_j] = calibFcn(pars_local,time,ens);
        [L,L_hpv,L_les,L_can] = calibLoss( ...
            calib_j,targetHPV,targetLesion,targetCancerASR);
        a2_out(j)   = A2(j);
        a3_out(j)   = A3(j);
        a4_out(j)   = A4(j);
        mu_out(j)   = M(j);
        Pcum_out(j) = P(j);
        loss_out(j) = L;
        Lhpv_out(j) = L_hpv;
        Lles_out(j) = L_les;
        Lcan_out(j) = L_can;
        lesion_out(j) = calib_j.lesion;
        cancer_out(j) = calib_j.cancer;
        hpv_out(j,:)  = calib_j.hpv;
        cancerTrajAll(:,j) = mean(calib_j.cancerASRByYear,2,'omitnan');
        lesionTrajAll(:,j) = mean(calib_j.lesionByYear,2,'omitnan');
        hpvTrajAll(:,:,j)  = mean(calib_j.hpvByYear,3,'omitnan');
    end
end
%% PARAMETER TABLE
ParameterTable = table( ...
    (1:nJobs)', ...
    a2_out,a3_out,a4_out,mu_out,Pcum_out, ...
    loss_out,Lhpv_out,Lles_out,Lcan_out, ...
    lesion_out,cancer_out, ...
    hpv_out(:,1),hpv_out(:,2),hpv_out(:,3), ...
    hpv_out(:,4),hpv_out(:,5),hpv_out(:,6), ...
    'VariableNames',{ ...
    'Combination','a2','a3','a4','mu','Pcum', ...
    'L','L_HPV','L_Lesion','L_Cancer', ...
    'LesionPrevalence','ASR', ...
    'HPV_0_14','HPV_15_24','HPV_25_34', ...
    'HPV_35_44','HPV_45_54','HPV_55plus'});
%% RANK PARAMETER COMBINATIONS
if any(~isfinite(loss_out))
    error('ZambiaCalibration:NonfiniteGridLoss', ...
        'At least one grid loss is nonfinite. Inspect outputs before ranking.');
end
[lossSorted,rankAll] = sort(loss_out,'ascend');
best10idx = rankAll(1:min(10,nJobs));
best10    = ParameterTable(best10idx,:);
best10.Rank = (1:height(best10))';
best10 = movevars(best10,'Rank','Before',1);
fprintf('\n============================================\n');
fprintf('BEST 10 PARAMETER COMBINATIONS\n');
fprintf('============================================\n');
disp(best10);
%% BEST PARAMETER SET
bestIdx = best10idx(1);
pars_best = struct( ...
    'a2',   a2_out(bestIdx), ...
    'a3',   a3_out(bestIdx), ...
    'a4',   a4_out(bestIdx), ...
    'mu',   mu_out(bestIdx), ...
    'Pcum', Pcum_out(bestIdx));
fprintf('\n============================================\n');
fprintf('BEST-FITTING PARAMETER SET\n');
fprintf('============================================\n');
disp(ParameterTable(bestIdx,:));
fprintf('\nTargets:\n');
fprintf('HPV 25-34: %.3f\n',targetHPV(3));
fprintf('HPV 35-44: %.3f\n',targetHPV(4));
fprintf('Lesion 25-49: %.4f\n',targetLesion);
fprintf('Cancer ASR: %.2f\n',targetCancerASR);
%% RERUN BEST PARAMETER SET
% Keeps all ens stochastic trajectories for the best-fitting combination.
fprintf('\nRerunning best parameter set to retain ensemble trajectories...\n');
[~,~,~,~,~,bestCalib] = calibFcn(pars_best,time,ens);
%% HISTORICAL PERIOD: 1990-2020
yearsHist = 1990:2020;
histIdx   = (burntime+1):time;
if numel(histIdx) ~= numel(yearsHist)
    error('Historical-period indexing does not match 1990-2020.');
end
cancerTrajHist = cancerTrajAll(histIdx,:);
lesionTrajHist = lesionTrajAll(histIdx,:);
hpvTrajHist    = hpvTrajAll(histIdx,:,:);
%% SAVE
outputFile = fullfile(outputDir, sprintf('calibration_allcomb_trajectories_ens%d.mat', ens));
save(outputFile, ...
    'ParameterTable','best10','best10idx','rankAll','lossSorted', ...
    'pars_best','bestCalib', ...
    'a2Grid','a3Grid','a4Grid','muGrid','PcumGrid', ...
    'targetHPV','targetLesion','targetCancerASR','ageLabels', ...
    'ens','time','burntime','yearsHist','histIdx', ...
    'cancerTrajAll','lesionTrajAll','hpvTrajAll', ...
    'cancerTrajHist','lesionTrajHist','hpvTrajHist', ...
    'hpv_out','lesion_out','cancer_out','loss_out', ...
    'Lhpv_out','Lles_out','Lcan_out', ...
    '-v7.3');
fprintf('\nSaved calibration results to %s\n',outputFile);
%% LOCAL CALIBRATION LOSS
function [L,L_hpv,L_les,L_can] = calibLoss( ...
    metrics,targetHPV,targetLesion,targetCancerASR)
    eps0 = 0.01;
    % HPV calibration:
    % group 3 = age 25-34
    % group 4 = age 35-44
    L_hpv = mean( ...
        ((metrics.hpv(3:4) - targetHPV(3:4)) ./ ...
        (targetHPV(3:4) + eps0)).^2);
    % Lesion prevalence among women aged 25-49
    L_les = ((metrics.lesion - targetLesion) / targetLesion)^2;
    % Cervical cancer age-standardized incidence
    L_can = ((metrics.cancer - targetCancerASR) / targetCancerASR)^2;
    L = L_hpv + L_les + L_can;
end
