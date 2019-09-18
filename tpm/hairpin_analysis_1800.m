names={{'1800L-5nM'},{'12RSS-5nM'},{'NoRSS-5nM'}};
loops_comp=[];
pccut_comp=[];
totallength_comp=[];
for j=1:3
    load(fullfile('C:\hairpins\analyzed\1800L',names{j}{:},'analysis_270.mat'))
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
for j=1:3
subplot(3,1,j)
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
set(gca, 'xscale','linear', 'yscale','linear','xlim',[0 1550],'xtick',[0 250 500 750 1000 1250 1500],'ylim',[0 6],'ytick',[0 2 4 6],'xscale','linear','fontsize',18);
%leg = legend('No bead loss','Bead loss',1);
%set(leg,'Box','on','fontsize', 18);
figureSize = [600 600];
set(150, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
end
saveas(150, fullfile('C:\hairpins\analyzed\1800L','pc.pdf'), 'pdf');
saveas(150, fullfile('C:\hairpins\analyzed\1800L','pc.fig'), 'fig');

%%
figure(170)
mc=meds_comp{1};
mcc=meds_cut_comp{1};
[ma,mb]=hist(mc(mc>150),[150:5:300]);
[mca,mcb]=hist(mcc(mcc>150),[150:5:300]);
plot([196 196], [0 20],'b--','LineWidth',3)
hold on
bar(mb,ma,'r')
bar(mcb,mca,'k')

xlabel('Median dwell length','fontsize',30);
ylabel('Occurences','fontsize',30);
set(gca, 'xscale','linear','xlim',[100 300], 'ylim',[0 9],'ytick',[0 4 8],'yscale','linear','xscale','linear','fontsize',24);
%leg = legend('No bead loss','Bead loss',1);
%set(leg,'Box','on','fontsize', 18);
figureSize = [600 200];
set(170, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
saveas(170, fullfile('C:\hairpins\analyzed\1800L','pct.pdf'), 'pdf');
saveas(170, fullfile('C:\hairpins\analyzed\1800L','pct.fig'), 'fig');