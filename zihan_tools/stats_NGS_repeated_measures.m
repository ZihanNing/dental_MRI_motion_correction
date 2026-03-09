%% stats_NGS_repeated_measures.m
% Load NGS_all_cases.mat, run repeated-measures ANOVA per region, Tukey post-hoc,
% print results, and plot diagrams.
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-03

clear; clc; close all;

% -----------------------------
% 1) Load pre-saved table
% -----------------------------
matPath = '/home/zn23/Data/ddMRI/NGS_all_cases.mat';

if ~exist(matPath, 'file')
    error('Cannot find MAT file: %s', matPath);
end

S = load(matPath);

% Expect variable "T" saved in MAT
if ~isfield(S, 'T')
    error('MAT file does not contain variable "T". Please check %s', matPath);
end
T = S.T;

if ~istable(T)
    error('Loaded variable "T" is not a MATLAB table.');
end

fprintf('Loaded %d cases from: %s\n', height(T), matPath);

% -----------------------------
% 2) Define regions and columns
% -----------------------------
% IMPORTANT: These names must match your column names in T.
% (They should, if you used the previous collector script.)
reconLevels = {'Aq','MoCo_Full','MoCo_Upper','MoCo_Lower'};

regionSpecs = struct( ...
    'name', {'FullFOV','UpperJaw','LowerJaw'}, ...
    'cols', { ...
        {'NoMoCo_Aq_FullFOV','MoCo_FullFOV_FullFOV','MoCo_UpperJaw_FullFOV','MoCo_LowerJaw_FullFOV'}, ...
        {'NoMoCo_Aq_UpperJaw','MoCo_FullFOV_UpperJaw','MoCo_UpperJaw_UpperJaw','MoCo_LowerJaw_UpperJaw'}, ...
        {'NoMoCo_Aq_LowerJaw','MoCo_FullFOV_LowerJaw','MoCo_UpperJaw_LowerJaw','MoCo_LowerJaw_LowerJaw'} ...
    } ...
);

% Sanity check: verify columns exist
for r = 1:numel(regionSpecs)
    cols = regionSpecs(r).cols;
    missing = cols(~ismember(cols, T.Properties.VariableNames));
    if ~isempty(missing)
        fprintf('\n[ERROR] Missing columns for region %s:\n', regionSpecs(r).name);
        disp(missing(:));
        error('Fix the column names above to match T.Properties.VariableNames.');
    end
end

% -----------------------------
% 3) Run stats + plots per region
% -----------------------------
results = struct();

for r = 1:numel(regionSpecs)
    regionName = regionSpecs(r).name;
    cols = regionSpecs(r).cols;

    fprintf('\n============================================================\n');
    fprintf('Region: %s\n', regionName);
    fprintf('============================================================\n');

    % Extract data table for fitrm (only the 4 measure columns)
    data = T(:, cols);

    % Remove rows with any missing values in these columns
    data = rmmissing(data);
    n = height(data);
    fprintf('N used (after removing missing): %d\n', n);

    % Within-subject design table (factor levels)
    within = table(categorical(reconLevels(:)), 'VariableNames', {'Reconstruction'});

    % Repeated measures model
    % Use "first-last ~ 1" syntax
    rm = fitrm(data, sprintf('%s-%s ~ 1', cols{1}, cols{end}), 'WithinDesign', within);

    % Repeated measures ANOVA table
    ranovatbl = ranova(rm);
    disp('Repeated-measures ANOVA (ranova):');
    disp(ranovatbl);

    % Tukey pairwise comparisons
    % Note: multcompare supports several corrections; Tukey is appropriate for all pairwise.
    posthoc = multcompare(rm, 'Reconstruction', 'ComparisonType', 'tukey-kramer');
    disp('Post-hoc pairwise comparisons (Tukey-Kramer):');
    disp(posthoc);

    % Store
    results(r).region = regionName;
    results(r).n = n;
    results(r).ranova = ranovatbl;
    results(r).posthoc = posthoc;

    % -----------------------------
    % 4) Plot diagrams
    % -----------------------------
    % Convert to numeric matrix for plotting (n x 4)
    Y = table2array(data);

    % (A) Boxplot
    figure('Name', ['NGS Boxplot - ' regionName], 'Color', 'w');
    boxplot(Y, 'Labels', reconLevels);
    ylabel('NGS');
    title(['NGS across reconstructions - ' regionName], 'Interpreter', 'none');
    grid on;

    % -----------------------------
    % (B) Mean ± SD BAR PLOT WITH SIGNIFICANCE
    % -----------------------------
    Y = table2array(data);

    mu = mean(Y, 1, 'omitnan');
    sd = std(Y, 0, 1, 'omitnan');

    figure('Name', ['NGS Mean±SD - ' regionName], 'Color', 'w');
    hold on;

    x = 1:numel(reconLevels);

    % Bar plot
    bh = bar(x, mu);
    bh.FaceColor = [0.6 0.6 0.6];
    bh.EdgeColor = 'k';

    % Error bars
    errorbar(x, mu, sd, 'k', 'LineStyle', 'none', 'LineWidth', 1);

    xticks(x);
    xticklabels(reconLevels);
    ylabel('NGS');
    title(['Mean \pm SD NGS - ' regionName], 'Interpreter', 'tex');
    grid on;

    % -----------------------------
    % Add significance annotations
    % -----------------------------
    sigPairs = posthoc(posthoc.pValue < 0.05, :);

    if ~isempty(sigPairs)

        yMax = max(mu + sd);
        yOffset = 0.02 * range(mu + sd);   % spacing between bars
        currentHeight = yMax + yOffset;

        for s = 1:height(sigPairs)

            g1 = sigPairs.Reconstruction_1(s);
            g2 = sigPairs.Reconstruction_2(s);
            p  = sigPairs.pValue(s);

            idx1 = find(strcmp(reconLevels, char(g1)));
            idx2 = find(strcmp(reconLevels, char(g2)));

            if isempty(idx1) || isempty(idx2)
                continue
            end

            % Draw horizontal line
            plot([idx1 idx1 idx2 idx2], ...
                 [currentHeight currentHeight+yOffset currentHeight+yOffset currentHeight], ...
                 'k', 'LineWidth', 1.2);

            % Determine stars
            if p < 0.001
                stars = '***';
            elseif p < 0.01
                stars = '**';
            else
                stars = '*';
            end

            text(mean([idx1 idx2]), currentHeight + yOffset*1.2, ...
                 stars, 'HorizontalAlignment', 'center', ...
                 'FontSize', 12, 'FontWeight', 'bold');

            currentHeight = currentHeight + 3*yOffset;
        end
    end

    hold off;

    % (C) Spaghetti plot (each subject)
    figure('Name', ['NGS Subject Lines - ' regionName], 'Color', 'w');
    plot(x, Y', '-o');  % each column becomes one x-point, transpose to plot per subject
    xlim([0.5, numel(reconLevels)+0.5]);
    xticks(x);
    xticklabels(reconLevels);
    ylabel('NGS');
    title(['Per-case NGS (paired) - ' regionName], 'Interpreter', 'none');
    grid on;

end

% -----------------------------
% 5) Save stats results
% -----------------------------
outStatsMat = '/home/zn23/Data/ddMRI/NGS_stats_repeated_measures.mat';
save(outStatsMat, 'results', 'regionSpecs', 'reconLevels', 'matPath');
fprintf('\nSaved stats struct to: %s\n', outStatsMat);