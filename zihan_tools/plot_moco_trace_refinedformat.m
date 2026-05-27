% This is in use 24-Mar-2026 by Zihan

%%
clear
clc
close all

%%
close all

%% load the T (T1 & T2)
T = squeeze(T); 
T_rmDiv = 0.5*(T(1:2:end-1,:) + T(2:2:end,:));

t = 1:size(T_rmDiv,1);

%% load the T (PD)
T = squeeze(T); 
T_rmDiv = T;

t = 1:size(T_rmDiv,1);

%% load the T (T2)
T = squeeze(T); 
T_rmDiv = T(1:round(size(T,1))/2,:);

t = 1:size(T_rmDiv,1);

   

%% Define 6 distinct colors (MATLAB default palette)
figure('Color','w','Units','centimeters','Position',[10 10 25 5]);
colors = lines(6);

% Left y-axis: translation (mm)
yyaxis left
hold on
for i = 1:3
    plot(t, T_rmDiv(:,i), '-', 'LineWidth', 1.5, 'Color', colors(i,:));
end
ylabel('Translation (mm)');
grid on; box on;

% Right y-axis: rotation (deg)
yyaxis right
hold on
for i = 4:6
    plot(t, T_rmDiv(:,i), '-', 'LineWidth', 1.5, 'Color', colors(i,:));
end
ylabel('Rotation (deg)');

xlabel('Shot index');
% title('Motion trace (translation + rotation)');

legend({'Tra-LR','Tra-AP','Tra-FH','Rot-FH','Rot-AP','Rot-LR'}, ...
       'Location','bestoutside');
