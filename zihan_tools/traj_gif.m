%% create_sampling_pattern_gif.m
% Create animated GIFs for MRI trajectory sampling patterns:
% 1) Linear Cartesian
% 2) DISORDER-like Cartesian (uniformly distributed within k-space)
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 08-Apr-2026

clear; close all; clc;

%% Parameters
N = 16;                         % 16 x 16 grid
gifDelay = 0.15;               % delay time between frames (s)
markerSize = 220;              % size of sampled points
lineWidthGrid = 0.8;
outDir = pwd;

% Colors
bgColor      = [1 1 1];
gridColor    = [0.85 0.85 0.85];
unsampledCol = [0.88 0.88 0.88];
sampledCol   = [0.10 0.45 0.90];
currentCol   = [0.95 0.35 0.20];
textCol      = [0.15 0.15 0.15];

%% Generate trajectories
% Grid coordinates:
% x: column index, y: row index
[xg, yg] = meshgrid(1:N, 1:N);

% --- 1. Linear Cartesian: one row from left to right
linear_order = [];
targetRow = 1; % "most edge row"; you may change to N if preferred
for c = 1:N
    linear_order = [linear_order; c, targetRow];
end

% --- 2. DISORDER Cartesian: one shot, uniformly distributed across k-space
% Here we place one point per row, spread across columns with a permutation.
disorder_order = [1,1;4,2;9,3;15,1;...
    3,5;7,5;10,7;13,8;...
    1,12;7,9;12,11;16,11;...
    4,13;5,16;10,15;14,13];

%% Create GIFs
create_sampling_gif( ...
    linear_order, N, ...
    fullfile(outDir, 'linear_cartesian.gif'), ...
    'Linear Cartesian Sampling', ...
    bgColor, gridColor, unsampledCol, sampledCol, currentCol, textCol, ...
    markerSize, lineWidthGrid, gifDelay);

create_sampling_gif( ...
    disorder_order, N, ...
    fullfile(outDir, 'disorder_cartesian.gif'), ...
    'DISORDER Cartesian Sampling', ...
    bgColor, gridColor, unsampledCol, sampledCol, currentCol, textCol, ...
    markerSize, lineWidthGrid, gifDelay);

disp('GIFs created successfully.');

function create_sampling_gif(order, N, gifName, figTitle, ...
    bgColor, gridColor, unsampledCol, sampledCol, currentCol, textCol, ...
    markerSize, lineWidthGrid, gifDelay)

    fig = figure('Color', bgColor, 'Position', [100 100 700 700]);
    ax = axes(fig);
    hold(ax, 'on');
    axis(ax, 'equal');
    xlim(ax, [0.5 N+0.5]);
    ylim(ax, [0.5 N+0.5]);
    set(ax, 'YDir', 'reverse'); % makes row 1 appear at top
    ax.XTick = 1:N;
    ax.YTick = 1:N;
    ax.FontSize = 12;
    ax.LineWidth = 1.0;
    ax.Box = 'on';
    ax.XColor = textCol;
    ax.YColor = textCol;
    xlabel(ax, 'k_x', 'FontSize', 24, 'Color', textCol);
    ylabel(ax, 'k_y', 'FontSize', 24, 'Color', textCol);
    title(ax, figTitle, 'FontSize', 24, 'FontWeight', 'bold', 'Color', textCol);

    % Draw grid
    for k = 0.5:1:(N+0.5)
        plot(ax, [0.5 N+0.5], [k k], '-', 'Color', gridColor, 'LineWidth', lineWidthGrid);
        plot(ax, [k k], [0.5 N+0.5], '-', 'Color', gridColor, 'LineWidth', lineWidthGrid);
    end

    % Draw all grid points as unsampled background
    [xg, yg] = meshgrid(1:N, 1:N);
    scatter(ax, xg(:), yg(:), markerSize*0.35, ...
        'o', 'MarkerFaceColor', unsampledCol, ...
        'MarkerEdgeColor', 'none');

    nFrames = size(order, 1);
    
    pause(1)

    for f = 1:nFrames
        cla(ax);
        hold(ax, 'on');
        axis(ax, 'equal');
        xlim(ax, [0.5 N+0.5]);
        ylim(ax, [0.5 N+0.5]);
        set(ax, 'YDir', 'reverse');
        ax.XTick = 1:N;
        ax.YTick = 1:N;
        ax.FontSize = 12;
        ax.LineWidth = 1.0;
        ax.Box = 'on';
        ax.XColor = textCol;
        ax.YColor = textCol;
        xlabel(ax, 'k_x', 'FontSize', 24, 'Color', textCol);
        ylabel(ax, 'k_y', 'FontSize', 24, 'Color', textCol);
        title(ax, figTitle, 'FontSize', 24, 'FontWeight', 'bold', 'Color', textCol);

        % Grid
        for k = 0.5:1:(N+0.5)
            plot(ax, [0.5 N+0.5], [k k], '-', 'Color', gridColor, 'LineWidth', lineWidthGrid);
            plot(ax, [k k], [0.5 N+0.5], '-', 'Color', gridColor, 'LineWidth', lineWidthGrid);
        end

        % All points in background
        scatter(ax, xg(:), yg(:), markerSize*0.35, ...
            'o', 'MarkerFaceColor', unsampledCol, ...
            'MarkerEdgeColor', 'none');

        % Sampled points so far
        sampled = order(1:f, :);
        scatter(ax, sampled(:,1), sampled(:,2), markerSize, ...
            'o', 'MarkerFaceColor', sampledCol, ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.2);

        % Current point highlight
        scatter(ax, order(f,1), order(f,2), markerSize*1.15, ...
            'o', 'MarkerFaceColor', currentCol, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2);

        % Optional trajectory line
        if f > 1
            plot(ax, sampled(:,1), sampled(:,2), '-', ...
                'Color', sampledCol*0.75, 'LineWidth', 2);
        end

        % Progress text
        text(0.7, N+0.15, sprintf('Shot progress: %d / %d', f, nFrames), ...
            'FontSize', 13, 'Color', textCol, 'FontWeight', 'bold');

        drawnow;

        % Capture frame
        frame = getframe(fig);
        img = frame2im(frame);
        [A, map] = rgb2ind(img, 256);

        if f == 1
            imwrite(A, map, gifName, 'gif', 'LoopCount', Inf, 'DelayTime', gifDelay);
        else
            imwrite(A, map, gifName, 'gif', 'WriteMode', 'append', 'DelayTime', gifDelay);
        end
    end
    
    pause(2)

    close(fig);
end