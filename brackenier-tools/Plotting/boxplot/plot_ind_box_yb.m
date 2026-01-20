
% y = rand([1 1000]);
% 
% figure
% 
% xCord = 2;
% width = 1;
% co = [0 1 0]

function [pBox] = plot_ind_box_yb(xCord,y,width, co)

y(isnan(y))=[];

% get percentiles
pt = prctile(y,[9 25 50 75 91]); 
means = mean(y);
med = pt(3);
IQR = (pt(4)-pt(2));

% box filled with color 
xBox = [xCord-width/2 xCord-width/2 xCord+width/2 xCord+width/2];
yBox = [pt(2)  pt(4)  pt(4)  pt(2)];

pBox = fill(xBox,yBox,co);            
set(pBox,'FaceAlpha',1,'LineWidth',1.5); 
hold on

% draw the median
line([xCord-width/2 xCord+width/2], [med med],'color','k', 'LineWidth', 2);

% draw whiskers
plot([xCord xCord],[pt(4) pt(5)],'k--','LineWidth',2)
plot([xCord xCord],[pt(1) pt(2)],'k--','LineWidth',2)
plot([xCord xCord],[pt(1) pt(2)],'k--','LineWidth',2)

plot([xCord-width/4 xCord+width/4],[pt(1) pt(1)],'k--','LineWidth',2)
plot([xCord-width/4 xCord+width/4],[pt(5) pt(5)],'k--','LineWidth',2)


