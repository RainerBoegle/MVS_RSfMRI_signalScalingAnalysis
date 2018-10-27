%% select figure windows to plot
FigNumsForPlotting = 1:4; %[2,3,4];
NamesForPlotting = {'dLog2_AggregateData';'dLog2_MedSL_Data';'dLog2_MedS_Data';'OverlapBase2Encode_dlog2_green1MedS_blue2Agg_red4MedSL'};

%% try to plot
%pdf
for IndPlot= 1:length(FigNumsForPlotting)
    print(FigNumsForPlotting(IndPlot),[NamesForPlotting{IndPlot},'.pdf'],'-dpdf','-noui','-opengl','-r600','-loose');
end

% %tiff
% for IndPlot= 1:length(FigNumsForPlotting)
%     print(FigNumsForPlotting(IndPlot),[NamesForPlotting{IndPlot},'.tiff'],'-dtiff','-noui','-painters','-r600');
% end

%-depsc2 
for IndPlot= 1:length(FigNumsForPlotting)
    print(FigNumsForPlotting(IndPlot),[NamesForPlotting{IndPlot},'.ps'],'-depsc2 ','-noui','-painters','-r600');
end

% %-depsc2 
% for IndPlot= 1:length(FigNumsForPlotting)
%     print(FigNumsForPlotting(IndPlot),[NamesForPlotting{IndPlot},'_2.ps'],'-depsc2 ','-noui','-opengl','-r600');
% end
