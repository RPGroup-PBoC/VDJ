close all
clear all
c=0;
mon12=[];
moff12=[];
moff12_fin=[];
mon12_fin=[];
mon23=[];
moff23=[];
moff23_fin=[];
mon23_fin=[];
meanon12_std_all=[];
meanon12_std_fin_all=[];
meanoff12_std_all=[];
meanoff12_std_fin_all=[];
meanon23_std_all=[];
meanon23_std_fin_all=[];
meanoff23_std_all=[];
meanoff23_std_fin_all=[];
numoffs12=[];
numons12=[];
numoffs12_fin=[];
numons12_fin=[];
numoffs23=[];
numons23=[];
numoffs23_fin=[];
numons23_fin=[];

sumcalc = @(x) sum(x);
meanoncalc = @(x) mean(x);
for k=157.5:0.25:157.5
c=c+1;
conc=[1, 2.5, 5, 10, 25, 50];
names={{'12RSS-1nM'},{'12RSS-2.5nM'},{'12RSS-5nM'},{'12RSS-10nM'},{'12RSS-25nM'},{'12RSS-50nM'}};
offtime_nolac_comp_c=[];
ontime_nolac_comp_c=[];
for ff=1:length(names)
currexp=names{ff};
currexp=currexp{1};
%load data file
load( fullfile('C:\RAGDynamics',currexp,['\analysis_',num2str(k),'.mat']))
offtime_nolac_comp_c=cat(1,offtime_nolac_comp_c,offtime_nolac_comp);
ontime_nolac_comp_c=cat(1,ontime_nolac_comp_c,ontime_nolac_comp);
end
timeoff_c=sum(offtime_nolac_comp_c);
%so transition density is number of events divided by timeoff
transdens=length(ontime_nolac_comp_c)./timeoff_c;
%mean length of a fake state
mean_fake=mean(ontime_nolac_comp_c);
meanon_corr=[];
meanon_corr_fin=[];
meanoff_corr=[];
meanoff_corr_fin=[];
meanon=[];
meanoff=[];
meanon_fin=[];
meanoff_fin=[];
wall=[];
wall_fin=[];
ontime_comp_fin=[];
offtime_comp_fin=[];
meanon_std=[];
meanoff_std=[];
meanon_fin_std=[];
meanoff_fin_std=[];
numoffs=[];
for ff=1:length(names)
    currexp=names{ff};
    currexp=currexp{1};
    %load data file
load( fullfile('C:\RAGDynamics',currexp,['\analysis_',num2str(k),'.mat']))
    meanon(ff)=mean(ontime_comp);
    meanoff(ff)=mean(offtime_comp);
    
    sumoff=sum(offtime_comp);
    sizeoff=length(offtime_comp);
    ontime_comp_fin=ontime_comp(ontime_accept==0 | ontime_accept==2);
    offtime_comp_fin=offtime_comp(offtime_accept==0 | offtime_accept==2);
    meanon_fin(ff)=mean(ontime_comp_fin);
    meanoff_fin(ff)=mean(offtime_comp_fin);
    sumoff_fin=sum(offtime_comp_fin);
    sizeoff_fin=length(offtime_comp_fin);
    
    numfakes=transdens*sum(offtime_comp);
    w=(length(ontime_comp)-numfakes)./length(ontime_comp);
    meanon_corr(ff)=(mean(ontime_comp)-((1-w)*mean_fake))./w;
    meanoff_corr(ff)=(sumoff+mean_fake*numfakes)./(sizeoff-numfakes);
    wall(ff)=w;
    meanon_std(ff)=std(((bootstrp(1000,meanoncalc,ontime_comp./30)-(1-w).*bootstrp(1000,meanoncalc,ontime_nolac_comp_c./30))./w));
    meanoff_std(ff)=std((bootstrp(1000,sumcalc,offtime_comp./30))./(sizeoff-numfakes));
    numons(ff)=length(ontime_comp)-numfakes;
    numoffs(ff)=sizeoff-numfakes;
    
    numfakes_fin=transdens*sum(offtime_comp);
    w=(length(ontime_comp(ontime_accept==0 | ontime_accept==2))-numfakes)./length(ontime_comp(ontime_accept==0 | ontime_accept==2));
    meanon_corr_fin(ff)=(mean(ontime_comp(ontime_accept==0 | ontime_accept==2))-((1-w)*mean_fake))./w;
    meanoff_corr_fin(ff)=(sumoff_fin+mean_fake*numfakes)./(sizeoff_fin-numfakes);
    wall_fin(ff)=w;
    numons_fin(ff)=length(ontime_comp(ontime_accept==0 | ontime_accept==2))-numfakes;
    numoffs_fin(ff)=sizeoff_fin-numfakes;
    meanon_fin_std(ff)=std(((bootstrp(1000,meanoncalc,ontime_comp./30)-(1-w).*bootstrp(1000,meanoncalc,ontime_nolac_comp_c./30))./w));
    meanoff_fin_std(ff)=std((bootstrp(1000,sumcalc,offtime_comp_fin./30))./(sizeoff-numfakes));
end
meanon_corr12=meanon_corr./30;
meanon_corr_fin12=meanon_corr_fin./30;
meanoff12=meanoff./30;
meanoff_fin12=meanoff_fin./30;
meanoff_corr_fin12=meanoff_corr_fin./30;
meanoff_corr12=meanoff_corr./30;

meanon12_std_all(c,:)=meanon_std;
meanon12_std_fin_all(c,:)=meanon_fin_std;
meanoff12_std_all(c,:)=meanoff_std;
meanoff12_std_fin_all(c,:)=meanoff_fin_std;
numoffs12(c,:)=numoffs;
numons12(c,:)=numons;
numoffs12_fin(c,:)=numoffs_fin;
numons12_fin(c,:)=numons_fin;

mon12(c,:)=meanon_corr12;
mon12_fin(c,:)=meanon_corr_fin12;
moff12(c,:)=meanoff_corr12;
moff12_fin(c,:)=meanoff_corr_fin12;
wtot12(c,:)=wall;
wtot12_fin(c,:)=wall_fin;

names={{'23RSS-1nM'},{'23RSS-2.5nM'},{'23RSS-5nM'},{'23RSS-10nM'},{'23RSS-25nM'},{'23RSS-50nM'}};
offtime_nolac_comp_c=[];
ontime_nolac_comp_c=[];
for ff=1:length(names)
currexp=names{ff};
currexp=currexp{1};
%load data file
load( fullfile('C:\RAGDynamics',currexp,['\analysis_',num2str(k),'.mat']))
offtime_nolac_comp_c=cat(1,offtime_nolac_comp_c,offtime_nolac_comp);
ontime_nolac_comp_c=cat(1,ontime_nolac_comp_c,ontime_nolac_comp);
end
timeoff_c=sum(offtime_nolac_comp_c);
%so transition density is number of events divided by timeoff
transdens=length(ontime_nolac_comp_c)./timeoff_c;
%mean length of a fake state
mean_fake=mean(ontime_nolac_comp_c);
meanon_corr=[];
meanon_corr_fin=[];
meanoff_corr=[];
meanoff_corr_fin=[];
meanon=[];
meanoff=[];
meanon_fin=[];
meanoff_fin=[];
wall=[];
wall_fin=[];
ontime_comp_fin=[];
offtime_comp_fin=[];
meanon_std=[];
meanoff_std=[];
meanon_fin_std=[];
meanoff_fin_std=[];
for ff=1:length(names)
    currexp=names{ff};
    currexp=currexp{1};
    %load data file
load( fullfile('C:\RAGDynamics',currexp,['\analysis_',num2str(k),'.mat']))
    meanon(ff)=mean(ontime_comp);
    meanoff(ff)=mean(offtime_comp);
    
    sumoff=sum(offtime_comp);
    sizeoff=length(offtime_comp);
    ontime_comp_fin=ontime_comp(ontime_accept==0 | ontime_accept==2);
    offtime_comp_fin=offtime_comp(offtime_accept==0 | offtime_accept==2);
    meanon_fin(ff)=mean(ontime_comp_fin);
    meanoff_fin(ff)=mean(offtime_comp_fin);
    sumoff_fin=sum(offtime_comp_fin);
    sizeoff_fin=length(offtime_comp_fin);
    
    numfakes=transdens*sum(offtime_comp);
    w=(length(ontime_comp)-numfakes)./length(ontime_comp);
    meanon_corr(ff)=(mean(ontime_comp)-((1-w)*mean_fake))./w;
    meanoff_corr(ff)=(sumoff+mean_fake*numfakes)./(sizeoff-numfakes);
    wall(ff)=w;
    meanon_std(ff)=std(((bootstrp(1000,meanoncalc,ontime_comp./30)-(1-w).*bootstrp(1000,meanoncalc,ontime_nolac_comp_c./30))./w));
    meanoff_std(ff)=std((bootstrp(1000,sumcalc,offtime_comp./30))./(sizeoff-numfakes));
    numons(ff)=length(ontime_comp)-numfakes;
    numoffs(ff)=sizeoff-numfakes;
    
    numfakes_fin=transdens*sum(offtime_comp);
    w=(length(ontime_comp(ontime_accept==0 | ontime_accept==2))-numfakes)./length(ontime_comp(ontime_accept==0 | ontime_accept==2));
    meanon_corr_fin(ff)=(mean(ontime_comp(ontime_accept==0 | ontime_accept==2))-((1-w)*mean_fake))./w;
    meanoff_corr_fin(ff)=(sumoff_fin+mean_fake*numfakes)./(sizeoff_fin-numfakes);
    wall_fin(ff)=w;
    numons_fin(ff)=length(ontime_comp(ontime_accept==0 | ontime_accept==2))-numfakes;
    numoffs_fin(ff)=sizeoff_fin-numfakes;
    meanon_fin_std(ff)=std(((bootstrp(1000,meanoncalc,ontime_comp./30)-(1-w).*bootstrp(1000,meanoncalc,ontime_nolac_comp_c./30))./w));
    meanoff_fin_std(ff)=std((bootstrp(1000,sumcalc,offtime_comp_fin./30))./(sizeoff-numfakes));
end
meanon_corr23=meanon_corr./30;
meanon_corr_fin23=meanon_corr_fin./30;
meanoff23=meanoff./30;
meanoff_fin23=meanoff_fin./30;
meanoff_corr_fin23=meanoff_corr_fin./30;
meanoff_corr23=meanoff_corr./30;

meanon23_std_all(c,:)=meanon_std;
meanon23_std_fin_all(c,:)=meanon_fin_std;
meanoff23_std_all(c,:)=meanoff_std;
meanoff23_std_fin_all(c,:)=meanoff_fin_std;
numoffs23(c,:)=numoffs;
numons23(c,:)=numons;
numoffs23_fin(c,:)=numoffs_fin;
numons23_fin(c,:)=numons_fin;

mon23(c,:)=meanon_corr23;
mon23_fin(c,:)=meanon_corr_fin23;
moff23(c,:)=meanoff_corr23;
moff23_fin(c,:)=meanoff_corr_fin23;
wtot23(c,:)=wall;
wtot23_fin(c,:)=wall_fin;
end

%%
figure;hold;
for k=1:5
mm23=mon23_fin(k,:);
mm12=mon12_fin(k,:);
ww23=wtot23_fin(k,:);
ww12=wtot12_fin(k,:);

mmf23=moff23_fin(k,:);
mmf12=moff12_fin(k,:);
ww23=wtot23_fin(k,:);
ww12=wtot12_fin(k,:);
plot(conc(ww23>0.66), mm23(ww23>0.66),'o','color','k')
plot(conc(ww12>0.66), mm12(ww12>0.66),'o','color','r')
plot(conc(ww23>0.66), mmf23(ww23>0.66),'.','color','k')
plot(conc(ww12>0.66), mmf12(ww12>0.66),'.','color','r')
end
%%
figure;hold;
for k=1:5
mm23=mon23(k,:);
mm12=mon12(k,:);
ww23=wtot23(k,:);
ww12=wtot12(k,:);

mmf23=moff23(k,:);
mmf12=moff12(k,:);
ww23=wtot23(k,:);
ww12=wtot12(k,:);
plot(conc(ww23>0.66), mm23(ww23>0.66),'o','color','k')
plot(conc(ww12>0.66), mm12(ww12>0.66),'o','color','r')
plot(conc(ww23>0.66), mmf23(ww23>0.66),'.','color','k')
plot(conc(ww12>0.66), mmf12(ww12>0.66),'.','color','r')
end
%%
close all
lim=0.5
figure(190);hold;
for k=3:3
mm23=mon23_fin(k,:);
mm12=mon12_fin(k,:);
ww23=wtot23_fin(k,:);
ww12=wtot12_fin(k,:);
mm12error=meanon12_std_fin_all(k,:)
mm23error=meanon23_std_fin_all(k,:)
errorbar(conc(ww12>lim), mm12(ww12>lim),mm12error(ww12>lim),'o','color','k','MarkerSize',8,'MarkerFaceColor','k')
errorbar(conc(ww23>lim), mm23(ww23>lim),mm23error(ww23>lim),'o','color','r','MarkerSize',8,'MarkerFaceColor','r')

mmf23=moff23_fin(k,:);
mmf12=moff12_fin(k,:);
wwf23=wtot23_fin(k,:);
wwf12=wtot12_fin(k,:);
mmf12error=meanoff12_std_fin_all(k,:)
mmf23error=meanoff23_std_fin_all(k,:)
errorbar(conc(ww12>lim), mmf12(ww12>lim),mmf12error(ww12>lim),'s','color','k','MarkerSize',8,'MarkerFaceColor','k')
errorbar(conc(ww23>lim), mmf23(ww23>lim),mmf23error(ww23>lim),'s','color','r','MarkerSize',8,'MarkerFaceColor','r')
end
ca= conc(ww12>lim);
ma=mmf12(ww12>lim);
me=mmf12error(ww12>lim);
box on
figureSize = [700 700];
% set(gca, 'fontsize', 24,'Layer','top','ylim',[0 1500],'xlim',[0 55]);
xlabel('[RAG1/2c](nM)');
 ylabel('Mean dwell time (s)');

 set(190, 'Position', [100 50 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
 leg = legend('<\tau_{bound}> : 12 RSS','<\tau_{bound}> : 23 RSS','<\tau_{unbound}> : 12 RSS','<\tau_{unbound}> : 23 RSS',1);
 set(leg,'Box','on','fontsize', 24);
 
 %%
close all
lim=7
figure(190);hold;
for k=3:3
mm23=mon23_fin(k,:);
mm12=mon12_fin(k,:);
ww23=numons23_fin(k,:);
ww12=numons12_fin(k,:);
mm12error=meanon12_std_fin_all(k,:)
mm23error=meanon23_std_fin_all(k,:)
errorbar(conc(ww12>lim), mm12(ww12>lim),mm12error(ww12>lim),'o','color','k','MarkerSize',8)
errorbar(conc(ww23>lim), mm23(ww23>lim),mm23error(ww23>lim),'o','color','r','MarkerSize',8)

mmf23=moff23_fin(k,:);
mmf12=moff12_fin(k,:);
ww23=numoffs23_fin(k,:);
ww12=numoffs12_fin(k,:);
mmf12error=meanoff12_std_fin_all(k,:)
mmf23error=meanoff23_std_fin_all(k,:)
errorbar(conc(ww12>lim), mmf12(ww12>lim),mmf12error(ww12>lim),'s','color','k','MarkerSize',8)
errorbar(conc(ww23>lim), mmf23(ww23>lim),mmf23error(ww23>lim),'s','color','r','MarkerSize',8)
end
ca= conc(ww12>lim);
ma=mmf12(ww12>lim);
me=mmf12error(ww12>lim);
box on
figureSize = [700 700];
% set(gca, 'fontsize', 24,'Layer','top','ylim',[0 1500],'xlim',[0 55]);
xlabel('[RAG1/2c](nM)');
 ylabel('Mean dwell time (s)');

 set(190, 'Position', [100 50 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
 leg = legend('<\tau_{bound}> : 12 RSS','<\tau_{bound}> : 23 RSS','<\tau_{unbound}> : 12 RSS','<\tau_{unbound}> : 23 RSS',1);
 set(leg,'Box','on','fontsize', 24);
 
 %%
 close all
lim=0.1
figure(190);hold;
for k=3:3
mm23=mon23(k,:);
mm12=mon12(k,:);
ww23=wtot23(k,:);
ww12=wtot12(k,:);
mm12error=meanon12_std_all(k,:)
mm23error=meanon23_std_all(k,:)
errorbar(conc(ww12>lim), mm12(ww12>lim),mm12error(ww12>lim),'o','color','k','MarkerSize',8,'MarkerFaceColor','k')
errorbar(conc(ww23>lim), mm23(ww23>lim),mm23error(ww23>lim),'o','color','r','MarkerSize',8,'MarkerFaceColor','r')

mmf23=moff23(k,:);
mmf12=moff12(k,:);
wwf23=wtot23(k,:);
wwf12=wtot12(k,:);
mmf12error=meanoff12_std_all(k,:)
mmf23error=meanoff23_std_all(k,:)
errorbar(conc(ww12>lim), mmf12(ww12>lim),mmf12error(ww12>lim),'s','color','k','MarkerSize',8,'MarkerFaceColor','k')
errorbar(conc(ww23>lim), mmf23(ww23>lim),mmf23error(ww23>lim),'s','color','r','MarkerSize',8,'MarkerFaceColor','r')
end
box on
figureSize = [700 700];
 %set(gca, 'fontsize', 24,'Layer','top','ylim',[0 1500],'xlim',[0 55]);
xlabel('[RAG1/2c](nM)');
 ylabel('Mean dwell time (s)');

 set(190, 'Position', [100 50 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
 leg = legend('<\tau_{bound}> : 12 RSS','<\tau_{bound}> : 23 RSS','<\tau_{unbound}> : 12 RSS','<\tau_{unbound}> : 23 RSS',1);
 set(leg,'Box','on','fontsize', 24);
 
 %%
close all
lim=1
figure(190);hold;
for k=1:1
mm23=mon23(k,:);
mm12=mon12(k,:);
ww23=numons23(k,:);
ww12=numons12(k,:);
mm12error=meanon12_std_all(k,:)
mm23error=meanon23_std_all(k,:)
errorbar(conc(ww12>lim), mm12(ww12>lim),mm12error(ww12>lim),'o','color','k','MarkerSize',10,'MarkerFaceColor','k', 'LineWidth',1.7)
errorbar(conc(ww23>lim), mm23(ww23>lim),mm23error(ww23>lim),'o','color','r','MarkerSize',10,'MarkerFaceColor','r', 'LineWidth',1.7)

mmf23=moff23(k,:);
mmf12=moff12(k,:);
ww23=numoffs23(k,:);
ww12=numoffs12(k,:);
mmf12error=meanoff12_std_all(k,:)
mmf23error=meanoff23_std_all(k,:)
errorbar(conc(ww12>lim), mmf12(ww12>lim),mmf12error(ww12>lim),'s','color','k','MarkerSize',10, 'LineWidth',1.7)
errorbar(conc(ww23>lim), mmf23(ww23>lim),mmf23error(ww23>lim),'s','color','r','MarkerSize',10, 'LineWidth',1.7)
end
ca= conc(ww12>lim);
ma=mmf12(ww12>lim);
me=mmf12error(ww12>lim);
box on
figureSize = [600 400];
set(gca, 'fontsize', 26,'Layer','top','ylim',[0 4750],'xlim',[0 55]);
set(gca, 'YTick',[0,1000, 2000, 3000, 4000],'YTickLabel',{'0','1','2','3','4'})
xlabel('[RAG1/2c](nM)');
 ylabel('Mean dwell time (x10^3 s)');

 set(190, 'Position', [100 50 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
 leg = legend('<\tau_{bound}> : 12 RSS','<\tau_{bound}> : 23 RSS','<\tau_{unbound}> : 12 RSS','<\tau_{unbound}> : 23 RSS',1);
 set(leg,'Box','on','fontsize', 26);
 
 saveas(190, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics','meandwelltimeTEST.pdf'), 'pdf');
saveas(190, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics','meandwelltimeTEST.fig'), 'fig');


%% Hists (ACCEPT ALL transitions) 12 RSS
close all
names={{'12RSS-1nM'},{'12RSS-2.5nM'},{'12RSS-5nM'},{'12RSS-10nM'},{'12RSS-25nM'},{'12RSS-50nM'}}
figure(120)
for ff=1:length(names)
currexp=names{ff};
currexp=currexp{1};
%load data file
load( fullfile('C:\RAGDynamics',currexp,'\analysis_157.5.mat'))
subplot(length(names),2,2*ff-1)
[a,b]=hist(offtime_comp./30,[0:400:20000],'FaceColor','r','EdgeColor','k');
bar(b,a./sum(a(:)))
hf = findobj(gca,'Type','patch');
set(hf,'FaceColor','r','EdgeColor','k')
set(gca, 'xscale','linear', 'yscale','linear','xlim',[0 7000],'ylim',[0 1],'xscale','linear');
%leg = legend([names{ff}],'Location','NorthEast');
%set(leg,'Box','on','fontsize', 8);
text(0.975, 0.8,names{ff},'units','normalized','HorizontalAlignment', 'Right','FontSize',12,'EdgeColor','k');

if ff==1
    title('Dwell time in unbound state','fontsize',12)
end
subplot(length(names),2,2*ff)
[a,b]=hist(ontime_comp./30,[0:200:20000]);
bar(b,a./sum(a(:)))
hf = findobj(gca,'Type','patch');
set(hf,'FaceColor','r','EdgeColor','k')
set(gca, 'xscale','linear', 'yscale','linear','xlim',[0 5000],'ylim',[0 1],'xscale','linear');
if ff==1
    title('Dwell time in bound state','fontsize',12)
end
end
figureSize = [900 700];
set(120, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
saveas(120, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics','timehist_all12.pdf'), 'pdf');
saveas(120, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics','timehist_all12.fig'), 'fig');
%% Hists (ACCEPT ALL transitions) 23RSS
close all
names={{'23RSS-1nM'},{'23RSS-2.5nM'},{'23RSS-5nM'},{'23RSS-10nM'},{'23RSS-25nM'},{'23RSS-50nM'}}
figure(120)
for ff=1:length(names)
currexp=names{ff};
currexp=currexp{1};
%load data file
load( fullfile('C:\RAGDynamics',currexp,'\analysis_157.5.mat'))
subplot(length(names),2,2*ff-1)
[a,b]=hist(offtime_comp./30,[0:400:20000],'FaceColor','r','EdgeColor','k');
bar(b,a./sum(a(:)))
hf = findobj(gca,'Type','patch');
set(hf,'FaceColor','r','EdgeColor','k')
set(gca, 'xscale','linear', 'yscale','linear','xlim',[0 7000],'ylim',[0 1],'xscale','linear');
%leg = legend([names{ff}],'Location','NorthEast');
%set(leg,'Box','on','fontsize', 8);
text(0.975, 0.8,names{ff},'units','normalized','HorizontalAlignment', 'Right','FontSize',12,'EdgeColor','k');

if ff==1
    title('Dwell time in unbound state','fontsize',12)
end
subplot(length(names),2,2*ff)
max(ontime_comp./30)
[a,b]=hist(ontime_comp./30,[0:200:20000]);
bar(b,a./sum(a(:)))
hf = findobj(gca,'Type','patch');
set(hf,'FaceColor','r','EdgeColor','k')
set(gca, 'xscale','linear', 'yscale','linear','xlim',[0 5000],'ylim',[0 1],'xscale','linear');
if ff==1
    title('Dwell time in bound state','fontsize',12)
end
end
figureSize = [900 700];
set(120, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
saveas(120, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics','timehist_all23.pdf'), 'pdf');
saveas(120, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics','timehist_all23.fig'), 'fig');

%% Fraction binding times (ACCEPT All data)

names={{'23RSS-1nM'},{'23RSS-2.5nM'},{'23RSS-5nM'},{'23RSS-10nM'},{'23RSS-25nM'},{'23RSS-50nM'}};
conc=[1, 2.5, 5, 10, 25, 50];
for ff=1:length(names)
currexp=names{ff};
currexp=currexp{1};
%load data file
load( fullfile('C:\RAGDynamics',currexp,'\analysis_157.5.mat'))
toton(ff)=sum(ontime_comp);
totoff(ff)=sum(offtime_comp);
end
fracoff=totoff./(toton+totoff);
fracon=toton./(toton+totoff);
close all
c=35.07
xx=[00:.01:100];
fx=(xx/c)./(1+(xx/c));
figure(150)
plot(conc,fracoff,'ro')
hold
plot(conc,fracon,'ko')
plot(xx,fx,'k-')
plot(xx,1-fx,'r-')
%h=text(16,0.505,'\leftarrow 16.5nM','fontsize',16,'VerticalAlignment','Middle','HorizontalAlignment','Left');
%set(h,'rotation',90)
leg = legend('fraction of time unbound','fraction of time bound',1);
set(leg,'Box','on','fontsize', 18);
xlabel('RAG concentration (nM)','fontsize',18);
ylabel('Fractional occupancy','fontsize',18);
set(gca, 'xscale','linear', 'yscale','linear','xlim',[0 51],'ylim',[0 1],'xscale','linear','fontsize',18);
figureSize = [600 600];
set(150, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
saveas(150, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics','fracbinding23RSS.pdf'), 'pdf');
saveas(150, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics','fracbinding23RSS.fig'), 'fig');

%% Print all traces for each sample
close all
lowstateaverage=153;
highstateaverage=163;
names={{'23RSS-1nM'},{'23RSS-2.5nM'},{'23RSS-5nM'},{'23RSS-10nM'},{'23RSS-25nM'},{'23RSS-50nM'}};
for ff=1:length(names)
    close all
    currexp=names{ff};
    currexp=currexp{1};
    load( fullfile('C:\RAGDynamics',currexp,'\analysis_157.5.mat'))
    for j=1:length(statetrace_comp)
        seg=statetrace_comp{j};
        tra=trace_comp{j};
        seg_no=statetrace_nolac_comp{j};
        tra_no=trace_nolac_comp{j};
        seg(seg==3)=highstateaverage;
        seg(seg==2)=lowstateaverage;
        seg_no(seg_no==3)=highstateaverage;
        seg_no(seg_no==2)=lowstateaverage;
        figure(j)
          subplot(2,1,2);
          plot(seg,'k--')
          hold
          plot(tra,'r-')
          plot([1 length(tra)+1],[158 158],'b--')
          set(gca, 'xscale','linear', 'yscale','linear','ylim',[135 195],'xscale','linear','xlim',[1 length(tra)]);
          subplot(2,1,1);
          plot(seg_no,'k--')
          hold
          plot(tra_no,'r-')
          plot([1 length(tra_no)+1],[158 158],'b--')
          set(gca, 'xscale','linear', 'yscale','linear','ylim',[135 195],'xscale','linear','xlim',[1 length(tra_no)]);
          if accepted(j)==0
          leg = legend('This trace is accepted for analysis');
          elseif accepted(j)==1
          leg = legend('This trace is NOT accepted (no transitions during data)');
          elseif accepted(j)==2
          leg = legend('This trace is NOT accepted (transitions during calibration)');
          elseif accepted(j)==3
          leg = legend('This trace is NOT accepted (failed both tests)');
          end
          set(leg,'Box','on','fontsize', 16,'Location','NorthOutside');
          saveas(j, fullfile('C:\Users\Brewster\Documents\My Dropbox\RAGDynamics',currexp,['trace' num2str(j,'%2.3d') '.pdf']))
    end
end


