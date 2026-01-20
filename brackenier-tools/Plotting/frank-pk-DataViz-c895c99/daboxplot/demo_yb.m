condition_names = {'Water', 'Land', 'Moon', 'Hyperspace'};

% an alternative color scheme for some plots
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30];  
   
figure();

h = daboxplot(data1,'scatter',2,'whiskers',0,'boxalpha',0.7,...
    'xtlabels', condition_names); 
ylabel('Performance');
xl = xlim; xlim([xl(1), xl(2)+0.75]);       % make space for the legend
legend([h.bx(1,:)],group_names);            % add the legend manually
set(gca,'FontSize',9);

