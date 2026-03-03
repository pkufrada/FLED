% KRAS_Iter_withAbsolute.m (Clinical-only, absolute-aware EFIRM fitting)
% Function: Optimize fragment proportion fitting using absolute KRAS copy numbers

%% Step 0: Create Output Directory + Backup Current Script
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outdir = fullfile(pwd, ['results_' timestamp]);
if ~exist(outdir, 'dir'); mkdir(outdir); end

% === Automatic script backup as a .txt file ===
currentScriptName = mfilename('fullpath');
[~, scriptBaseName] = fileparts(currentScriptName);
backupName = fullfile(outdir, ['log_' scriptBaseName, '.txt']);
try
    copyfile([currentScriptName, '.m'], backupName);
    fprintf('✔ Current script backed up to: %s\n', backupName);
catch
    warning('⚠ Unable to backup current script; mfilename may not recognize the script.');
end

%% Step 1: Read Clinical Data + Absolute KRAS
%% ==================== original EFIRM reading (-nA) ===================================
% clinicalTable = {
%     {'control', 'NC', 60.7470, 3.3240, 56.3755, 1.3595, 67.4005, 3.0825},...
%     {'control', 'NC', 56.5160, 1.0050, 54.9910, 0.6330, 63.4430, 4.1040},...
%     {'control', 'PC', 77.1750, 3.2840, 59.3625, 0.9725, 75.0805, 2.0945},...
%     {'clinical sample', 'PDAC_12', 57.4970, 1.8780, 57.2720, 5.1100, 65.6880, 1.6400},...
%     {'clinical sample', 'PDAC_3642', 58.6465, 3.9675, 58.7690, 6.3320, 43.9950, 19.6520},...
%     {'clinical sample', 'PDAC_3814', 57.6275, 4.9735, 67.7485, 15.3945, 82.5250, 14.2250},...
%     {'clinical sample', 'PDAC_442', 62.2925, 3.2045, 60.3155, 3.1915, 63.4955, 4.7715},...
%     {'clinical sample', 'PDAC_6944', 70.0390, 4.6780, 67.9255, 7.6875, 76.6985, 4.2275},...
%     {'clinical sample', 'PDAC_7082', 69.7280, 4.3780, 68.5095, 5.2615, 76.0630, 3.4280}
% };
%% ==================== subtracted EFIRM reading (-nA) ===================================
% clinicalTable = {
%     {'clinical sample', 'PDAC_12',    0.981,   0.033229528,  2.281, 0.092924297, 2.245, 0.025849976}, ...
%     {'clinical sample', 'PDAC_3642',  2.1305,  0.051409135,  3.778, 0.106666667, 0,     0.261745726}, ...
%     {'clinical sample', 'PDAC_3814',  1.1115,  0.086500165, 12.7575, 0.26879627, 19.082, 0.216554013}, ...
%     {'clinical sample', 'PDAC_442',   5.7765,  0.054640942,  5.3245, 0.054305842, 0.0525, 0.108455506}, ...
%     {'clinical sample', 'PDAC_6944', 13.523,  0.081176522, 12.9345, 0.113471147, 13.2555, 0.051226901}, ...
%     {'clinical sample', 'PDAC_7082', 13.212,  0.070281334, 13.5185, 0.087232967, 12.62,  0.053988078}
% };
%% ==================== calibrated EFIRM cps/uL ===================================
clinicalTable = {
    {'clinical sample', 'PDAC-12',    0.911112802, 0.111754573,  1.706672095, 0.417586379,   1.755759353, 0.118724676}, ...
    {'clinical sample', 'PDAC-3642',  3.019298660, 0.172894599,  4.313284951, 0.479342311,   0.000000000, 1.202154944}, ...
    {'clinical sample', 'PDAC-3814',  7.686161273, 0.290909609, 71.070155160, 1.207925858, 111.114163500, 0.994596861}, ...
    {'clinical sample', 'PDAC-442',  26.384234920, 0.183763523, 19.592065310, 0.244041447,   0.201922416, 0.344316503}, ...
    {'clinical sample', 'PDAC-6944', 96.279196730, 0.273005609, 74.187496580, 0.509920516,  79.469600200, 0.235276706}, ...
    {'clinical sample', 'PDAC-7082', 93.943672650, 0.236363889, 77.437109880, 0.392010486,  75.562066790, 0.247958337}, ...
    {'clinical sample', 'PDAC-3438',  0.000000000, 2.414500000,  0.000000000, 4.196379026,   0.000000000, 6.350277008}, ...
    {'clinical sample', 'PDAC-4254',  0.000000000, 4.652500000,  0.000000000,14.018885700,   0.000000000,10.908137120}, ...
    {'clinical sample', 'PDAC-4486',  0.000000000, 3.207000000,  0.254570892, 1.771819324,   0.000000000, 8.118112188}, ...
    {'clinical sample', 'PDAC-478',   0.000000000, 3.971500000,  0.000000000,11.864909270,   0.000000000, 6.540785319}  ...
};
%% ==================== calibrated EFIRM cps/uL ===================================

clinData = cell2table(vertcat(clinicalTable{:}), ...
    'VariableNames', {'sampleCategory', 'SampleInfo', 'KRAS_Close_Avg', 'KRAS_Close_Stdev', ...
                      'KRAS_Moderate_Avg', 'KRAS_Moderate_Stdev', 'KRAS_Far_Avg', 'KRAS_Far_Stdev'});

isClinical = strcmp(clinData.sampleCategory, 'clinical sample');
sampleNames = clinData.SampleInfo(isClinical);
clinOnlyRaw = [clinData.KRAS_Close_Avg, clinData.KRAS_Moderate_Avg, clinData.KRAS_Far_Avg];
clinOnlyRaw = clinOnlyRaw(isClinical, :);

% KRAS_total options:
KRAS_total_ddPCR = [2.81, 161.5, 166, 91.95, 105.35, 101.85, 2.31, 1.05, 2.71, 1.69];  % Original ddPCR

%%========== Original EFIRM Current Readings ========================================%%
% Close_vals    = [57.497, 58.6465, 57.6275, 62.2925, 70.039, 69.728];
% Moderate_vals = [57.272, 58.769, 67.7485, 60.3155, 67.9255, 68.5095];
% Far_vals      = [65.688, 43.995, 82.525, 63.4955, 76.6985, 76.063];
%%========== Original EFIRM Current Readings End ====================================%%

%%========== Converted to cps/uL ========================================%%
% EFIRM signal (converted to ddPCR cps/uL) for each probe (Averages only)
Close_vals    = [0.911112802, 3.01929866, 7.686161273, 26.38423492, 96.27919673, 93.94367265, 0, 0, 0, 0];
Moderate_vals = [1.706672095, 4.313284951, 71.07015516, 19.59206531, 74.18749658, 77.43710988, 0, 0, 0.254570892, 0];
Far_vals      = [1.755759353, 0, 111.1141635, 0.201922416, 79.4696002, 75.56206679, 0, 0, 0, 0];
%%========== Converted to cps/uL End ====================================%%

%% ============= BEGIN PATCH 1: LOD + censored handling + weights =============
% Extract stdev from table (matching average order)
Close_sd    = clinData.KRAS_Close_Stdev(isClinical);
Moderate_sd = clinData.KRAS_Moderate_Stdev(isClinical);
Far_sd      = clinData.KRAS_Far_Stdev(isClinical);

Close_avg    = Close_vals(:);
Moderate_avg = Moderate_vals(:);
Far_avg      = Far_vals(:);

% --- Estimate per-probe LOD (can be replaced with fixed values) ---
nzClose    = Close_avg(Close_avg>0);    lodClose    = max(min(nzClose,[],'omitnan'), 1e-4);
nzModerate = Moderate_avg(Moderate_avg>0); lodModerate = max(min(nzModerate,[],'omitnan'), 1e-4);
nzFar      = Far_avg(Far_avg>0);        lodFar      = max(min(nzFar,[],'omitnan'), 1e-4);

% Replace 0 with half-LOD, retain a censoring mask
Close_is0    = (Close_avg==0);
Moderate_is0 = (Moderate_avg==0);
Far_is0      = (Far_avg==0);

Close_obs    = Close_avg;    Close_obs(Close_is0)    = 0.5*lodClose;
Moderate_obs = Moderate_avg; Moderate_obs(Moderate_is0) = 0.5*lodModerate;
Far_obs      = Far_avg;      Far_obs(Far_is0)        = 0.5*lodFar;

% Variance weighting (heteroscedasticity handling), down-weight censored points
W_close    = 1 ./ (Close_sd.^2    + 1e-6);
W_moderate = 1 ./ (Moderate_sd.^2 + 1e-6);
W_far      = 1 ./ (Far_sd.^2      + 1e-6);
down = 0.3;  % Weighting coefficient for censored points
W_close(Close_is0)         = W_close(Close_is0)         * down;
W_moderate(Moderate_is0)   = W_moderate(Moderate_is0)   * down;
W_far(Far_is0)             = W_far(Far_is0)             * down;

% Observation matrix (processed values) and weight matrix
X_obs = [Close_obs, Moderate_obs, Far_obs];             % N x 3 (Original space)
W_mat = [W_close,   W_moderate,  W_far  ];              % N x 3

% Fitting performed in log1p space
X_log = log1p(X_obs);          % N x 3
W_mat = W_mat ./ max(W_mat(:));% Normalize weights for numerical stability
%% ============= END PATCH 1 ==================================================

% EFIRM total signal (sum of 3 probes)
KRAS_total_sum3 = Close_vals + Moderate_vals + Far_vals;

% Clinical sample names (existing order)
sample_keys = {'PDAC-12', 'PDAC-3642', 'PDAC-3814', 'PDAC-442', ...
               'PDAC-6944', 'PDAC-7082', 'PDAC-3438', ...
               'PDAC-4254', 'PDAC-4486', 'PDAC-478'};

% Default selection:
KRAS_total = KRAS_total_sum3;  % <<< Toggleable

% Construct dictionary mapping
KRAS_total_dict = containers.Map(sample_keys, KRAS_total);
KRAS_total = zeros(numel(sampleNames), 1);
for i = 1:numel(sampleNames)
    key = sampleNames{i};
    if isKey(KRAS_total_dict, key)
        KRAS_total(i) = KRAS_total_dict(key);
    else
        KRAS_total(i) = NaN;
    end
end

%% ============= BEGIN PATCH 2: weighted log-space fitting ====================
% Use processed X_log as target; Y_abs initial values follow existing logic
X_raw = X_obs;      % Used only for final evaluation/export
X     = X_log;      % Fitting using log1p space
rng(0);
A_iter = rand(3,3) * 0.2 + 0.4;          % Initial response coefficients (gentle)
Y_prop_init = repmat([0.9 0.05 0.05], numel(KRAS_total), 1);
Y_abs = max(Y_prop_init .* KRAS_total, 1e-3);
max_iter = 50;
tol = 1e-8;
prev_loss = Inf;
loss_history = zeros(max_iter,1);
R2_history = zeros(max_iter,1);
proportion_history = zeros(max_iter, 3);
lambda_ridge = 1e-2;                     % Ridge regression regularization to stabilize condition numbers

% Pre-calculate weight diagonal matrices for each column (Weighted Least Squares)
W1 = diag(W_mat(:,1));
W2 = diag(W_mat(:,2));
W3 = diag(W_mat(:,3));

figure; hold on;
title('R^2 and Proportion Change per Iteration (log1p, WLS)');
xlabel('Iteration'); yyaxis left; ylabel('R^2'); ylim([0 1]);
r2_plot = plot(nan, nan, 'b-', 'DisplayName', 'R^2');
yyaxis right; ylabel('Mean Proportion'); ylim([0 1]);
p1 = plot(nan, nan, 'r-', 'DisplayName', 'Short');
p2 = plot(nan, nan, 'g-', 'DisplayName', 'Medium');
p3 = plot(nan, nan, 'm-', 'DisplayName', 'Long');
legend show;

%% ====== PATCH A: robust fit for A (3x3) and intercept b (3x1) in log1p space ======
% Uncensored mask
isC = [Close_is0, Moderate_is0, Far_is0];   % N x 3
U   = ~isC;
A_iter = eye(3);        % Initial values
b_iter = zeros(3,1);
delta = 0.5;            % Huber transition threshold
lambda_ridge = 1e-2;    % Ridge regularization

for j = 1:3
    yj = X_log(:,j);          
    Wj = W_mat(:,j);          
    mask = U(:,j);
    if nnz(mask) < 3, mask = true(size(mask)); end
    Z = [ones(sum(mask),1), Y_abs(mask,:)];   % (n_eff x 4) Intercept + Y_abs
    t = yj(mask);
    w = Wj(mask);
    beta = [0; A_iter(j,:)'];                 % Initial values for intercept and coefficients
    for ir = 1:10
        r = t - Z*beta;
        % Huber weighting
        s = median(abs(r)) + 1e-8;
        z = r / s;
        hub = ones(size(z)); idx = abs(z) > delta;
        hub(idx) = delta ./ abs(z(idx));
        W_eff = diag(w .* hub);
        beta = (Z' * W_eff * Z + lambda_ridge*eye(4)) \ (Z' * W_eff * t);
    end
    b_iter(j)   = beta(1);
    A_iter(j,:) = beta(2:4)';
end

for iter = 1:max_iter
    % === Update 3x3 A_iter; solve each of the three rows using weighted ridge regression (log space) ===
    Y = Y_abs;           
    % Row 1 (Close)
    At = (Y' * W1 * Y) + lambda_ridge * eye(3);
    bt = Y' * W1 * X(:,1);
    A_iter(1,:) = (At \ bt)';  
    % Row 2 (Moderate)
    At = (Y' * W2 * Y) + lambda_ridge * eye(3);
    bt = Y' * W2 * X(:,2);
    A_iter(2,:) = (At \ bt)';
    % Row 3 (Far)
    At = (Y' * W3 * Y) + lambda_ridge * eye(3);
    bt = Y' * W3 * X(:,3);
    A_iter(3,:) = (At \ bt)';

    % === Fixed A_iter, solve for Y_abs (per-sample weighted NNLS, non-negative) ===
    %% ====== PATCH B: per-sample solve on simplex (p>=eps, sum p=1) ======
    opts = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',200);
    epsilon = 1e-4;  Aeq = [1,1,1];  beq = 1;  lb = [epsilon,epsilon,epsilon];  ub = [1,1,1];
    Y_abs_new = zeros(size(Y_abs));
    for i = 1:size(X,1)
        xi_log = X_log(i,:)';                 % 3 x 1
        wi     = W_mat(i,:)';                 % 3 x 1
        ktotal = KRAS_total(i);
        fun = @(p) obj_huber_p(p, A_iter, b_iter, xi_log, wi, ktotal, delta);
        p0  = Y_abs(i,:); p0 = p0./max(sum(p0),eps);
        [p_opt, ~] = fmincon(fun, p0, [],[], Aeq, beq, lb, ub, [], opts);
        Y_abs_new(i,:) = ktotal * p_opt;
    end

    % Maintain total per sample = KRAS_total
    scale_factor = KRAS_total ./ max(sum(Y_abs_new,2), eps);
    Y_abs_new = Y_abs_new .* scale_factor;

    % === Evaluation: Compare in log space ===
    X_fit_log = Y_abs_new * A_iter';      % Linear approximation
    RMSE_iter = sqrt( mean( (X(:) - X_fit_log(:)).^2 .* repelem(W_mat(:),1,1) ) );  
    SSres = sum( ((X(:) - X_fit_log(:)).^2) .* W_mat(:) );
    SStot = sum( ((X(:) - mean(X(:))).^2) .* W_mat(:) );
    R2 = 1 - SSres / max(SStot, 1e-8);
    loss_history(iter) = RMSE_iter;
    R2_history(iter) = R2;
    Y_prop = Y_abs_new ./ max(sum(Y_abs_new, 2), eps);
    proportion_history(iter,:) = mean(Y_prop, 1);
    set(r2_plot, 'XData', 1:iter, 'YData', R2_history(1:iter));
    set(p1, 'XData', 1:iter, 'YData', proportion_history(1:iter,1));
    set(p2, 'XData', 1:iter, 'YData', proportion_history(1:iter,2));
    set(p3, 'XData', 1:iter, 'YData', proportion_history(1:iter,3));
    drawnow;

    if abs(prev_loss - RMSE_iter) < tol && iter > 1
        loss_history = loss_history(1:iter);
        R2_history = R2_history(1:iter);
        proportion_history = proportion_history(1:iter,:);
        Y_abs = Y_abs_new;
        break;
    end

    % Slight perturbation to avoid local minima
    Y_abs = max(Y_abs_new + 0.05*mean(Y_abs_new(:))*randn(size(Y_abs_new)), 1e-6);
    scale_factor = KRAS_total ./ max(sum(Y_abs,2), eps);
    Y_abs = Y_abs .* scale_factor;
    prev_loss = RMSE_iter;
end

Y_abs_final = Y_abs;
%% ============= END PATCH 2 ==================================================

%% ====== PATCH C: metrics & export (intercept-aware) ======
Y_prop_final = Y_abs_final ./ max(sum(Y_abs_final, 2), eps);
% Prediction (log space, including intercept)
X_fit_log = (Y_abs_final * A_iter') + b_iter';   % N x 3

% Weighted RMSE (log space), per sample
RMSE_w = sqrt( mean( ((X_log - X_fit_log).^2) .* W_mat, 2, 'omitnan') );
RMSE = RMSE_w;   % Backward compatibility

% Similarity: Bray-Curtis (more robust)
brayCurtisSim = 1 - sum(abs(Y_prop_final - repmat([0.9 0.05 0.05], size(Y_prop_final,1),1)),2) ...
                  ./ (sum(Y_prop_final,2) + 1 + eps);

% Pearson(log1p)
pearsonR = zeros(size(Y_prop_final,1),1);
for i = 1:size(Y_prop_final,1)
    a = log1p([0.9 0.05 0.05]);   
    b = log1p(Y_prop_final(i,:));
    R = corrcoef(a(:), b(:));
    pearsonR(i) = R(1,2);
end

% Backward compatibility: ensures subsequent code referencing cosineSim remains functional
cosineSim = brayCurtisSim;   

T_err = table(sampleNames, RMSE, brayCurtisSim, pearsonR, ...
    'VariableNames', {'Sample','RMSE_weighted_log','Similarity_BrayCurtis','Pearson_log1p'});
writetable(T_err, fullfile(outdir, 'error_metrics_per_sample.csv'));
%% ============= END PATCH C ==================================================

%% Step 4: Visualization and Export
% Set all subsequent figures to docked mode
set(0, 'DefaultFigureWindowStyle', 'docked')

fitY_abs = Y_abs_final;  % Use final absolute values after iteration
predictedY = Y_abs_final ./ max(sum(Y_abs_final, 2), eps);
fragNames = {'Short','Medium','Long'};

set(groot, 'defaultAxesFontSize', 24);
set(groot, 'defaultAxesLineWidth', 1.5);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultBarLineWidth', 1.5);
set(groot, 'defaultAxesFontWeight', 'bold');
set(groot, 'defaultLineMarkerSize', 10);

% Color settings
colorMap = lines(3);
cmap = lines(numel(sampleNames));

% === Data Export ===
T_prop = array2table(predictedY, 'VariableNames', fragNames);
T_prop.Sample = sampleNames;
writetable(T_prop, fullfile(outdir, 'final_predicted_fragment_proportion.csv'));

T_abs = array2table(Y_abs_final, 'VariableNames', fragNames);
T_abs.Sample = sampleNames;
writetable(T_abs, fullfile(outdir, 'final_absolute_fragment_concentration.csv'));

T_fitAbs = array2table(fitY_abs, 'VariableNames', fragNames);
T_fitAbs.Sample = sampleNames;
writetable(T_fitAbs, fullfile(outdir, 'predicted_absolute_concentration.csv'));

T_r2 = table((1:length(R2_history))', R2_history, 'VariableNames', {'Iteration','R_squared'});
writetable(T_r2, fullfile(outdir, 'R2_history.csv'));

% === Plot Export ===
% R2 Plot
figure; bar(R2_history, 0.4, 'FaceColor', colorMap(1,:));
ylabel('R^2', 'FontWeight', 'bold'); xlabel('Iteration', 'FontWeight', 'bold'); title('R^2 per Iteration');
grid on;
saveas(gcf, fullfile(outdir, 'iterative_R2_convergence.png'));

% Absolute Bar Chart
figure;
b = bar(fitY_abs, 1, 'grouped');
for k = 1:3
    b(k).FaceColor = colorMap(k,:);
end
title('Predicted Absolute Fragment Concentrations');
ylabel('KRAS EFIRM (cps/uL)', 'FontWeight', 'bold'); xlabel('Sample', 'FontWeight', 'bold');
set(gca, 'XTick', 1:numel(sampleNames), 'XTickLabel', sampleNames, 'FontWeight', 'bold'); xtickangle(45);
grid on;
legend(fragNames, 'Location', 'northeast');
saveas(gcf, fullfile(outdir, 'fragment_absolute_barplot.png'));

% Proportion Bar Chart
figure;
b = bar(predictedY, 1, 'grouped');
for k = 1:3
    b(k).FaceColor = colorMap(k,:);
end
title('Predicted Fragment Proportions after Iteration');
ylabel('Proportion', 'FontWeight', 'bold'); xlabel('Sample', 'FontWeight', 'bold');
set(gca, 'XTick', 1:numel(sampleNames), 'XTickLabel', sampleNames, 'FontWeight', 'bold'); xtickangle(45);
grid on;
legend(fragNames, 'Location', 'northeast');
saveas(gcf, fullfile(outdir, 'fragment_barplot_proportion.png'));

% Box Chart with points
Y_long = reshape(Y_abs_final', [], 1);
Group_frag = repmat(fragNames', numel(sampleNames), 1);
Group_sample = repelem(sampleNames, 3);
fragCats = categorical(Group_frag, fragNames, 'Ordinal', true);
sampleCats = categorical(Group_sample, sampleNames, 'Ordinal', true);

figure('Units','pixels','Position',[100, 100, 900, 1100]); 
hold on;
boxchart(fragCats, Y_long, 'BoxWidth', 0.4, 'BoxFaceColor', [0.6 0.6 0.85]);
legendHandles = []; legendLabels = {};
for i = 1:length(Y_long)
    x = double(fragCats(i)) + 0.12 * (rand - 0.5);
    y = Y_long(i);
    sName = string(sampleCats(i));
    sIdx = find(sampleNames == sName);
    C = cmap(sIdx, :);
    h = scatter(x, y, 60, 'filled', 'MarkerFaceColor', C);
    if ~ismember(char(sName), legendLabels)
        legendHandles(end+1) = h;
        legendLabels{end+1} = char(sName);
    end
end
lgd = legend(legendHandles, legendLabels, 'Location', 'northeastoutside');
lgd.Title.String = 'Sample';
title('Absolute Fragment Concentrations by Class');
ylabel('KRAS EFIRM (cps/uL)', 'FontWeight', 'bold'); xlabel('Fragment Class', 'FontWeight', 'bold');
grid on;
saveas(gcf, fullfile(outdir, 'boxchart_final_absolute.png'));

% RMSE Bar Chart
figure;
bar(RMSE, 0.4, 'FaceColor', colorMap(2,:));
title('RMSE per Sample (Absolute)'); ylabel('RMSE', 'FontWeight', 'bold'); xlabel('Sample', 'FontWeight', 'bold');
set(gca, 'XTick', 1:numel(sampleNames), 'XTickLabel', sampleNames, 'FontWeight', 'bold'); xtickangle(45);
grid on;
saveas(gcf, fullfile(outdir, 'bar_RMSE_absolute.png'));

% Cosine Similarity Bar Chart
figure;
bar(cosineSim, 0.4, 'FaceColor', colorMap(3,:));
title('Cosine Similarity per Sample (Absolute)'); ylabel('Cosine Similarity', 'FontWeight', 'bold'); xlabel('Sample', 'FontWeight', 'bold');
set(gca, 'XTick', 1:numel(sampleNames), 'XTickLabel', sampleNames, 'FontWeight', 'bold'); xtickangle(45);
grid on;
saveas(gcf, fullfile(outdir, 'bar_CosineSimilarity_absolute.png'));

%% ===== ADD: Residual diagnostics (histogram + residual vs predicted) =====
% Description:
% - Uses calculated X_fit_log (log-space prediction, including intercept b_iter)
% - Uses X_log (log-space observation) and W_mat (weights)
% - Plots weighted residual histogram and residual vs. predicted scatter plots
Res_raw = X_log - X_fit_log;              % N x 3 raw residuals (log space)
Wn      = W_mat ./ max(W_mat(:));         % Weights normalized to 0-1 for stability
Res_w   = Res_raw .* sqrt(Wn);            % Weighted residuals

% 1) Residual Histogram (stacked probes) + Normal Fitting Curve
figure('Units','pixels','Position',[100,100,900,600]);
h = histogram(Res_w(:), 30, 'FaceColor', [0.85 0.4 0.2], 'Normalization','pdf'); 
hold on;
% Normal distribution fit
mu = mean(Res_w(:));
sigma = std(Res_w(:));
x = linspace(min(Res_w(:)), max(Res_w(:)), 200);
y = normpdf(x, mu, sigma);
plot(x, y, 'b-', 'LineWidth', 2);
% Optional: KDE fit (smooth estimate)
[f, xi] = ksdensity(Res_w(:));
plot(xi, f, 'k--', 'LineWidth', 2);

title('Histogram of Weighted Residuals', 'FontSize', 22, 'FontWeight', 'bold');
xlabel('Weighted residuals in log1p space', 'FontWeight','bold');
ylabel('Density', 'FontWeight','bold');
legend({'Histogram (pdf)','Normal fit','KDE fit'}, 'Location','best');
grid on;
saveas(gcf, fullfile(outdir, 'hist_weighted_residuals_with_fit.png'));

% 2) Residual vs. Predicted (stacked probes)
figure;
scatter(X_fit_log(:), Res_w(:), 36, 'filled'); grid on;
xlabel('Predicted (log1p space)'); ylabel('Weighted residual');
title('Residuals vs Predicted (stacked probes)');
yline(0,'--');
saveas(gcf, fullfile(outdir, 'scatter_residual_vs_pred.png'));
%% ===== END ADD ============================================================

%% === Figure: usctDNA & mnctDNA by ddPCR status (box + scatter, separated) ===
% Dependencies: Y_abs_final, sampleNames, KRAS_total_ddPCR, outdir
% 1) Ensure ddPCR_aligned is aligned with sampleNames
if ~exist('ddPCR_aligned','var')
    ddPCR_keys = {'PDAC-12','PDAC-3642','PDAC-3814','PDAC-442','PDAC-6944','PDAC-7082','PDAC-3438','PDAC-4254','PDAC-4486','PDAC-478'};
    ddPCR_vals = KRAS_total_ddPCR(:);
    ddPCR_map  = containers.Map(ddPCR_keys, ddPCR_vals);
    ddPCR_aligned = nan(numel(sampleNames),1);
    for ii = 1:numel(sampleNames)
        k = sampleNames{ii};
        if isKey(ddPCR_map, k), ddPCR_aligned(ii) = ddPCR_map(k); end
    end
end

% 2) Data splitting: usctDNA=col 1, mnctDNA=col 3; ddPCR < 10 Negative, >= 10 Positive
us_all = Y_abs_final(:,1);     
mn_all = Y_abs_final(:,3);     
negMask = ddPCR_aligned < 10 & ~isnan(ddPCR_aligned);
posMask = ddPCR_aligned >= 10 & ~isnan(ddPCR_aligned);
us_neg = us_all(negMask);  us_pos = us_all(posMask);
mn_neg = mn_all(negMask);  mn_pos = mn_all(posMask);

% 4) Color mapping (consistent with full script)
us_fill  = [0.60 0.80 1.00];   us_edge = [0.10 0.30 0.80];  % Blue
mn_fill  = [1.00 0.60 0.80];   mn_edge = [0.80 0.00 0.50];  % Magenta

% 5) Plotting (Positions: 1=us-Neg, 2=us-Pos, 3=divider, 4=mn-Neg, 5=mn-Pos)
x_us_neg = 1; x_us_pos = 2; x_sep = 3; x_mn_neg = 4; x_mn_pos = 5;
figure('Units','pixels','Position',[120,120,900,650]); hold on;

% Boxcharts
boxchart(x_us_neg*ones(size(us_neg)), us_neg, 'BoxFaceColor', us_fill, 'BoxFaceAlpha', 0.55, 'LineWidth', 2, 'WhiskerLineColor', us_edge, 'BoxEdgeColor', us_edge);
boxchart(x_us_pos*ones(size(us_pos)), us_pos, 'BoxFaceColor', us_fill, 'BoxFaceAlpha', 0.55, 'LineWidth', 2, 'WhiskerLineColor', us_edge, 'BoxEdgeColor', us_edge);
boxchart(x_mn_neg*ones(size(mn_neg)), mn_neg, 'BoxFaceColor', mn_fill, 'BoxFaceAlpha', 0.55, 'LineWidth', 2, 'WhiskerLineColor', mn_edge, 'BoxEdgeColor', mn_edge);
boxchart(x_mn_pos*ones(size(mn_pos)), mn_pos, 'BoxFaceColor', mn_fill, 'BoxFaceAlpha', 0.55, 'LineWidth', 2, 'WhiskerLineColor', mn_edge, 'BoxEdgeColor', mn_edge);

% Scatter overlays
jit = 0.12;
scatter(x_us_neg + (rand(size(us_neg))-0.5)*jit, us_neg, 65, 'o', 'MarkerEdgeColor', us_edge, 'MarkerFaceColor','none', 'LineWidth',1.5);
scatter(x_us_pos + (rand(size(us_pos))-0.5)*jit, us_pos, 65, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', us_edge, 'LineWidth',0.5);
scatter(x_mn_neg + (rand(size(mn_neg))-0.5)*jit, mn_neg, 65, 'o', 'MarkerEdgeColor', mn_edge, 'MarkerFaceColor','none', 'LineWidth',1.5);
scatter(x_mn_pos + (rand(size(mn_pos))-0.5)*jit, mn_pos, 65, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', mn_edge, 'LineWidth',0.5);

% Divider
xline(x_sep,'--','Color',[0.5 0.5 0.5],'LineWidth',2);

set(gca,'XTick',[x_us_neg x_us_pos x_mn_neg x_mn_pos], 'XTickLabel',{'Neg','Pos','Neg','Pos'}, 'FontSize',24,'FontWeight','bold');
ylabel('KRAS EFIRM (cps/\muL)','FontSize',24,'FontWeight','bold');

% Top Labels
yl = ylim; yTop = yl(2);
text(mean([x_us_neg x_us_pos]), yTop*1.04, 'usctDNA', 'Color', us_edge, 'FontSize',22, 'FontWeight','bold', 'HorizontalAlignment','center');
text(mean([x_mn_neg x_mn_pos]), yTop*1.04, 'mnctDNA', 'Color', mn_edge, 'FontSize',22, 'FontWeight','bold', 'HorizontalAlignment','center');

grid on; box on;
ylim([0, yTop*1.22]);
saveas(gcf, fullfile(outdir, 'box_usctDNA_mnctDNA_by_ddPCR.png'));

%% Only display usctDNA and mnctDNA boxcharts
fragNames_SL = {'usctDNA','mnctDNA'};                  
Y_abs_SL = [Y_abs_final(:,1), Y_abs_final(:,3)];       
Y_long_SL = reshape(Y_abs_SL', [], 1);
Group_frag_SL = repmat(fragNames_SL', numel(sampleNames), 1);
Group_sample_SL = repelem(sampleNames, numel(fragNames_SL));
fragCats_SL = categorical(Group_frag_SL, fragNames_SL, 'Ordinal', true);
sampleCats_SL = categorical(Group_sample_SL, sampleNames, 'Ordinal', true);

figure('Units','pixels','Position',[100, 100, 900, 1100]);
hold on;
boxchart(fragCats_SL, Y_long_SL, 'BoxWidth', 0.4, 'BoxFaceColor', [0.6 0.6 0.85]);
legendHandles = []; legendLabels = {};
for i = 1:length(Y_long_SL)
    x = double(fragCats_SL(i)) + 0.12 * (rand - 0.5);
    y = Y_long_SL(i);
    sName = string(sampleCats_SL(i));
    sIdx = find(sampleNames == sName);
    C = cmap(sIdx, :);
    h = scatter(x, y, 60, 'filled', 'MarkerFaceColor', C);
    if ~ismember(char(sName), legendLabels)
        legendHandles(end+1) = h;
        legendLabels{end+1} = char(sName);
    end
end
lgd = legend(legendHandles, legendLabels, 'Location', 'northeast');
lgd.Title.String = 'Sample';
xlabel('Fragment Type','FontSize',22,'FontWeight','bold');
ylabel('KRAS EFIRM (cps/uL)','FontSize',22,'FontWeight','bold');
set(gca,'FontSize',20,'FontWeight','bold');
grid on;
saveas(gcf, fullfile(outdir, 'boxchart_final_absolute_usctDNA_mnctDNA.png'));

%% Histograms and QQ Plots
edges = 0:0.02:0.5;  
figure; histogram(RMSE, edges, 'FaceColor', colorMap(2,:));
title('Histogram of RMSE (Absolute)'); xlabel('RMSE', 'FontWeight', 'bold'); ylabel('Frequency'); grid on;
saveas(gcf, fullfile(outdir, 'histogram_RMSE_absolute.png'));

figure('Units','pixels','Position',[100, 100, 900, 900]); qqplot(RMSE); 
title('QQ Plot of RMSE (Absolute)'); grid on;
saveas(gcf, fullfile(outdir, 'qqplot_RMSE_absolute.png'));

edges = 0:0.02:1;
figure; histogram(cosineSim, edges, 'FaceColor', colorMap(3,:));
title('Histogram of Cosine Similarity (Absolute)'); xlabel('Cosine Similarity', 'FontWeight', 'bold'); ylabel('Frequency'); grid on;

figure('Units','pixels','Position',[100, 100, 900, 900]); qqplot(cosineSim); 
title('QQ Plot of Cosine Similarity (Absolute)'); grid on;
saveas(gcf, fullfile(outdir, 'qqplot_CosineSimilarity_absolute.png'));

%% ====== PATCH D: Stacked Short/Medium/Long with ddPCR overlay ======
% Align ddPCR array with sampleNames
ddPCR_keys = {'PDAC-12','PDAC-3642','PDAC-3814','PDAC-442','PDAC-6944','PDAC-7082','PDAC-3438','PDAC-4254','PDAC-4486','PDAC-478'};
ddPCR_vals = KRAS_total_ddPCR(:);                 
ddPCR_map  = containers.Map(ddPCR_keys, ddPCR_vals);
ddPCR_aligned = nan(numel(sampleNames),1);
for i = 1:numel(sampleNames)
    k = sampleNames{i};
    if isKey(ddPCR_map, k), ddPCR_aligned(i) = ddPCR_map(k); end
end

idx_show = 1:numel(sampleNames);
Xshort = Y_abs_final(idx_show,1);
Xmed   = Y_abs_final(idx_show,2);
Xlong  = Y_abs_final(idx_show,3);
ddpcr  = ddPCR_aligned(idx_show);
labs   = sampleNames(idx_show);

figure('Units','pixels','Position',[100,100,1100,700]); hold on; box on; grid on;
barH = bar([Xshort Xmed Xlong], 'stacked', 'BarWidth', 0.75);
barH(1).FaceColor = [0.30 0.50 0.85];
barH(2).FaceColor = [0.95 0.55 0.35];
barH(3).FaceColor = [0.35 0.75 0.55];
for k = 1:3, barH(k).EdgeColor = 'none'; end

plot(1:numel(ddpcr), ddpcr, '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'k');
set(gca, 'XTick', 1:numel(labs), 'XTickLabel', labs, 'XTickLabelRotation', 30, 'FontWeight', 'bold');
ylabel('Concentration (cps/\muL)', 'FontWeight','bold');
legend({'Short','Medium','Long','ddPCR'}, 'Location','northwest', 'Box','on');
saveas(gcf, fullfile(outdir, 'stacked_fragments_with_ddPCR.png'));

%% ====== PATCH E: EFIRM vs ddPCR correlation ======
EFIRM_total = sum(Y_abs_final,2);
validIdx = ~isnan(ddPCR_aligned) & ddPCR_aligned>0;   
x = ddPCR_aligned(validIdx);
y = EFIRM_total(validIdx);

[Rp, pval_p] = corr(x, y, 'Type','Pearson');
fprintf('Pearson R=%.3f (p=%.3g)\n', Rp, pval_p);

figure('Units','pixels','Position',[100,100,900,700]);
scatter(x, y, 80, 'filled'); hold on; lsline;
xlabel('ddPCR (cps/\muL)', 'FontWeight','bold');
ylabel('EFIRM total (cps/\muL)', 'FontWeight','bold');
title(sprintf('EFIRM total vs ddPCR (Pearson r=%.2f)', Rp), 'FontSize',22, 'FontWeight','bold');
grid on;
saveas(gcf, fullfile(outdir, 'scatter_corr_totalEFIRM_ddPCR.png'));

%% === New Plot: Stacked bars for usctDNA + mnctDNA with ddPCR line overlay ===
us_vals = Y_abs_final(:,1);
mn_vals = Y_abs_final(:,3);
ddPCR_map = containers.Map(sample_keys, KRAS_total_ddPCR);
ddPCR_vec = nan(numel(sampleNames),1);
for i = 1:numel(sampleNames)
    key = string(sampleNames{i});
    if isKey(ddPCR_map, key); ddPCR_vec(i) = ddPCR_map(key); end
end

fig = figure('Units','pixels','Position',[100,100,1200,700]); hold on; box on; grid on;
barH = bar([us_vals mn_vals], 'stacked', 'BarWidth', 0.75);
barH(1).FaceColor = [0.30 0.50 0.85]; barH(1).EdgeColor = 'none';
barH(2).FaceColor = [0.85 0.10 0.55]; barH(2).EdgeColor = 'none';

if any(~isnan(ddPCR_vec))
    plot(1:numel(sampleNames), ddPCR_vec, '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    legend({'usctDNA','mnctDNA','ddPCR'}, 'Location', 'northwest');
else
    legend({'usctDNA','mnctDNA'}, 'Location', 'northwest');
end
set(gca,'XTick',1:numel(sampleNames),'XTickLabel',sampleNames,'XTickLabelRotation',30,'FontWeight','bold');
ylabel('Concentration (cps/\muL)','FontWeight','bold');
saveas(fig, fullfile(outdir,'stacked_usctDNA_mnctDNA_with_ddPCR.png'));

%% ====== PATCH F: ROC-like analysis ======
labels = ddPCR_aligned > 5;   % Ground truth: ddPCR > 5 cps/uL considered Positive
scores = EFIRM_total;         
[Xroc,Yroc,T,AUC] = perfcurve(labels, scores, true);

figure('Units','pixels','Position',[100,100,900,700]);
plot(Xroc,Yroc,'-o','LineWidth',2); hold on; plot([0 1],[0 1],'--k');
xlabel('False Positive Rate (1-Specificity)','FontWeight','bold');
ylabel('True Positive Rate (Sensitivity)','FontWeight','bold');
title(sprintf('ROC of EFIRM vs ddPCR (AUC=%.2f)', AUC), 'FontSize',22,'FontWeight','bold');
grid on;
saveas(gcf, fullfile(outdir, 'roc_totalEFIRM_ddPCR.png'));

fprintf('Analysis complete. Results saved in: %s\n', outdir);
disp('All computations complete. Workspace cleared.');
clearvars -except outdir

function f = obj_huber_p(p, Aiter, biter, xlog, w, ktotal, delta)
    pred = biter + Aiter * (ktotal * p(:));   
    r = xlog - pred;
    idx = abs(r) <= delta;
    hub = zeros(size(r));
    hub(idx)  = 0.5 * r(idx).^2;
    hub(~idx) = delta*(abs(r(~idx)) - 0.5*delta);
    f = sum( w(:) .* hub(:) );
end