names={{'1200L-5nM'},{'1800L-5nM'},{'12RSS-5nM'},{'NoRSS-5nM'}};
loops_comp=[];
pccut_comp=[];
totallength_comp=[];
for j=1:4
    load(fullfile('Y:\VDJ\hairpins\analyzed\all',names{j}{:},'analysis_270.mat'))
    loops=[];
    pccut=[];
    cutframes=[];
    totallength=0;
    bb=0;
    bbb=0;
    meds=[];
    meds_cut=[];
    for i=1:length(statetrace_comp)
        currtrace=trace_comp{i};
        if sum(currtrace(1:1.5*10^4)>300)>0
            totallength=totallength+length(currtrace);
            %i
            %find all loops
            regs=[];
            if sum(statetrace_comp{i}==2)>0
                regs=regionprops(onstate_comp{i},'PixelList','Area','Centroid','BoundingBox');
                loops=cat(2,loops,[regs.Area]);
                for k=1:length(regs)
                    bb=bb+1;
                    meds(bb)=median(currtrace(regs(k).PixelList(:,2)));
                    %j
                    %i
                    median(currtrace(regs(k).PixelList(:,2)));
                end
                %    plot(trace_comp{i},'r')
                %    waitforbuttonpress
                %    close
                %    i
            end
            %find ending loops
            if sum(statetrace_comp{i}==2)>0
                laststep=max(find(onstate_comp{i}==1));
                thetrace=trace_comp{i};
                as=allstates_comp{i};
                regs2=[];
                % i
                if length(thetrace)>laststep
                    %i
                    %laststep+1
                    %length(as)
                    laststep_nextstate=max(find(as==as(laststep+1)));
                    if thetrace(laststep_nextstate)==0
                        regs2=regionprops(onstate_comp{i},'PixelList','Area','Centroid','BoundingBox');
                        pccut=cat(2,pccut,regs2(end).Area);
                        cutframes=cat(2,cutframes, i);
                        i
                        for k=1:length(regs2)
                            bbb=bbb+1;
                            meds_cut(bbb)=median(currtrace(regs(k).PixelList(:,2)));
                        end
                        %i
                        %plot(trace_comp{i},'k')
                        %median(currtrace(regs(k).PixelList(:,2)))
                        %waitforbuttonpress
                        %close
                        
                    end
                end
            end
        end
    end
    loops_comp{j}=loops;
    pccut_comp{j}=pccut;
    totallength_comp{j}=totallength;
    meds_comp{j}=meds;
    meds_cut_comp{j}=meds_cut;
end
%%
close all
figure(150)
for j=1:4
subplot(4,1,j)
currloop=loops_comp{j};
currloopcut=pccut_comp{j};
currmeds=meds_comp{j};
currmedscut=meds_cut_comp{j};

[loopa,loopb]=hist(currloop(currmeds>150)./30,[50:50:2000]);
[pca,pcb]=hist(currloopcut(currmedscut>150)./30,[50:50:2000]);
bar(loopb,loopa,'r')
hold on
bar(pcb,pca,'k')
xlabel('Paired complex dwell time','fontsize',18);
ylabel('Occurences','fontsize',18);
set(gca, 'xscale','linear', 'yscale','linear','ylim',[0 15],'xlim',[0 1550],'xtick',[0 250 500 750 1000 1250 1500],'xscale','linear','fontsize',18);
%leg = legend('No bead loss','Bead loss',1);
%set(leg,'Box','on','fontsize', 18);
figureSize = [600 600];
set(150, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
end
%saveas(150, fullfile('C:\hairpins\analyzed\1200L','pc.pdf'), 'pdf');
%saveas(150, fullfile('C:\hairpins\analyzed\1200L','pc.fig'), 'fig');

%%
figure(170)
mc=meds_comp{2};
mcc=meds_cut_comp{2};
[ma,mb]=hist(mc(mc>150),[150:4:300]);
[mca,mcb]=hist(mcc(mcc>150),[150:4:300]);
plot([245 245], [0 25],'b--','LineWidth',3)
hold on
bar(mb,ma,'r')
bar(mcb,mca,'k')

xlabel('Median dwell length','fontsize',30);
ylabel('Occurences','fontsize',30);
set(gca, 'xscale','linear','xlim',[100 300], 'ylim',[0 20],'ytick',[0 5 10 15 20],'yscale','linear','xscale','linear','fontsize',24);
%leg = legend('No bead loss','Bead loss',1);
%set(leg,'Box','on','fontsize', 18);
figureSize = [600 200];
set(170, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
%saveas(170, fullfile('C:\hairpins\analyzed\1200L','pct.pdf'), 'pdf');
%saveas(170, fullfile('C:\hairpins\analyzed\1200L','pct.fig'), 'fig');
%%
load('Y:\VDJ\hairpins\analyzed\all\1200L-5nM\analysis_270.mat')
si=length(statetrace_comp);
figue
for i=1:length(statetrace_comp)
    subplot(si,1,i);
    currtrace=trace_comp{i};
    median(currtrace(statetrace_comp{i}==3))
    if sum(currtrace(1:1.5*10^4)>300)>0
    plot(100.*statetrace_comp{i})
    set(gca,'ylim',[0 350])
    hold;
    plot(trace_comp{i},'r')
    plot(trace_nolac_comp{i},'g')
    waitforbuttonpress
    close all
    end
end
%% Trace figure
load('Y:\VDJ\hairpins\analyzed\all\1200L-5nM\analysis_270.mat')
currtrace=trace_comp{19};
currstate=statetrace_comp{19};
currstate(currstate==3)=300;
currstate(currstate==2)=250;
currstate(currstate==1)=nan;
load('Y:\VDJ\hairpins\analyzed\all\1800L-5nM\analysis_270.mat')
currtrace2=trace_comp{52};
currstate2=statetrace_comp{52};
currstate2(currstate2==3)=300;
currstate2(currstate2==2)=200;
currstate2(currstate2==1)=nan;
%%
close all
figure(150)
hold
plot([1:1:length(currstate)]./30,currtrace,'g','LineWidth',1.5)
plot([1:1:length(currstate2(75000:end))]./30,currtrace2(75000:end),'b','LineWidth',1.5)
plot([1:1:length(currstate2)]./30,245.*ones(length(currstate2),1),'g--','LineWidth',1.5)
plot([1:1:length(currstate2)]./30,196.*ones(length(currstate2),1),'b--','LineWidth',1.5)
set(gca,'ylim',[190 330],'xlim',[0 2500])
box on
figureSize = [400 600];
xlabel('Time (s)','fontsize',36);
ylabel('<RMS> (nm)','fontsize',36);
set(gca, 'fontsize',36);
saveas(150, fullfile('C:\hairpins\analyzed\all\traces.pdf'),'pdf');


%% print traces
load('Y:\VDJ\hairpins\analyzed\all\1800L-5nM\analysis_270.mat')
si=length(statetrace_comp);
k=0;
j=0;
figure(150)
hold
for i=1:2:si
    if j<5
        j=j+1;
        subplot(5,2,2*j-1);
        currtrace=trace_comp{i};
 %   median(currtrace(statetrace_comp{i}==3))
 %   if sum(currtrace(1:1.5*10^4)>300)>0
        plot(100.*statetrace_comp{i},'k','LineWidth',1.5)
        set(gca,'ylim',[0 350])
        hold on;
        plot(trace_comp{i},'r','LineWidth',1.5)
        plot(trace_nolac_comp{i},'g','LineWidth',1.5)
        
        if i+1<si+1
        subplot(5,2,2*j);
        currtrace=trace_comp{i};
 %   median(currtrace(statetrace_comp{i}==3))
 %   if sum(currtrace(1:1.5*10^4)>300)>0
        plot(100.*statetrace_comp{i+1},'k','LineWidth',1.5)
        set(gca,'ylim',[0 350])
        hold on;
        plot(trace_comp{i+1},'r','LineWidth',1.5)
        plot(trace_nolac_comp{i+1},'g','LineWidth',1.5)
        end
    end
    if j==5
        k=k+1;
        figureSize = [1200 1600];
        set(150, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
        saveas(150, fullfile('C:\hairpins\analyzed\all',['1800trajectories',num2str(k),'.pdf']), 'pdf');
        j=0
        close all
        figure(150)
        hold on
    end
end
k=k+1;
figureSize = [1200 1600];
set(150, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
saveas(150, fullfile('C:\hairpins\analyzed\all',['1800trajectories',num2str(k),'.pdf']), 'pdf');
j=0
close all
