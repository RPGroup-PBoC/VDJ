close all
hmgb=[5 10 25 50 112 225]
f=[539	2.03	0.66		4.81	0.91		6.17	0.8		16.27	1.16		23.27	1.26		40.67	1.16;
736	0.3	0.49		2.28	0.57		6.42	0.91		12.44	0.83		23.11	0.64		40.98	0.56;
946	0.18	0.85		1.29	1.25		8.05	2.01		13.18	0.62		24.54	0.62		36.26	0.88;
1124	4.57	0.54		7.01	0.56		10.09	0.65		18.47	1.04		29.57	0.96		43.96	0.96;
1316	0.99	1.31		5.61	1.03		8.15	0.97		14.28	1.17		24.23	1.1		35.56	0.74;
1521	2.73	0.64		4.42	0.98		8.94	1.11		14.63	0.82		24.77	0.94		34.49	0.91;
1717	2.72	0.45		6.55	0.69		11.88	0.48		20.51	0.63		26.88	0.57		41.69	0.47;
1910	1.35	1.34		6.24	1.11		10.05	1.05		16.94	1.16		26.4	0.7		43.52	0.77;
2077	2.24	0.79		5.78	0.81		12.58	0.57		20.13	0.5		33.6	0.42		43.75	0.44;
2280	3.14	0.44		8.82	1.02		15.63	0.92		26.02	0.98		37.5	0.87		47.64	1.12];
figure(120)
hold
k=9
cc=hsv(k);
box on
for i=1:9
    plotHandle = ploterr(hmgb,f(i,2:2:end),[],f(i,3:2:end),'o-');
        set(plotHandle(1),'LineWidth',1,'Color',cc(mod(i,k)+1,:),'MarkerSize',6,'MarkerFaceColor',cc(mod(i,k)+1,:),'MarkerEdgeColor',cc(mod(i,k)+1,:));
        set(plotHandle(2),'LineWidth',1,'Color',cc(mod(i,k)+1,:));
        set(get(get(plotHandle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%errorbar(hmgb,f(i,2:2:end),f(i,3:2:end),'ko-')
end
plotHandle=ploterr(hmgb,f(10,2:2:end),[],f(10,3:2:end),'o-')
set(plotHandle(1),'LineWidth',1,'Color','k','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k');
        set(plotHandle(2),'LineWidth',1,'Color','k');
        set(get(get(plotHandle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
figureSize = [600 350];
%set(gca, 'xscale','log', 'yscale','log','ylim',[3e-2 300],'xlim',[3e-2 300], 'fontsize', 16);
%set(gca, 'xscale','log', 'yscale','linear','ylim',[0.8 4],'xlim',[3e-2 11], 'fontsize', 28,'Layer','top');
%set(gca, 'XTick',[10.^-2,10.^-1,10.^0,10.^1,10.^2]);
%set(gca, 'YTick',[0,1,2,3,4,5,6,7,8,9,10,11,12]);)

%set(120,'Markersize',1)
xlabel('HMGB1 concentration','fontsize',17);
ylabel('Percent DNA compaction','fontsize',17);
set(gca, 'xscale','linear', 'yscale','linear','ylim',[-5 50],'fontsize', 15,'Layer','top');
set(120, 'Position', [100 50 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
set(gca, 'XTick',[0 50 100 150 200 250]);
set(gca, 'YTick',[0 10 20 30 40 50]);
leg = legend('539bp substrate','736bp substrate','946bp substrate','1124bp substrate','1316bp substrate','1521bp substrate','1717bp substrate','1910bp substrate','2077bp substrate','2280bp substrate',4);
set(leg,'Box','on','fontsize', 14);
saveas(120, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics\Figs','hmgb1compaction.pdf'), 'pdf');