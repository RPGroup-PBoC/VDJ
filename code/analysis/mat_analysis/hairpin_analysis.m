path = '/Users/soichihirokawa/Documents/vdj_recombination/analysis';
specific_path = 'Analyzed Data/new_analysis';
names={{'fix'}};
loops_comp=[];
pccut_comp=[];
totallength_comp=[];
for j=1:1
    data = load(...
        fullfile(path, specific_path, names{j}{:},'analysis_275.mat'));
    loops=[];
    pccut=[];
    cutframes=[];
    totallength=0;
    bb=0;
    bbb=0;
    meds=[];
    meds_cut=[];
    for i=1:length(data.statetrace_comp)
        currtrace=data.trace_comp{i};
        if sum(currtrace(1:1.5*10^4)>290)>0
            totallength=totallength+length(currtrace);
            %i
            %find all loops
            regs=[];
            if sum(data.statetrace_comp{i}==2)>0
                
                onstate_label=bwlabel(data.onstate_comp{i});
                % PixelList gives every pixel in every looped state
                regs=regionprops(onstate_label,'PixelList','Area','Centroid','BoundingBox');
                [regs.Area]./30
                loops=cat(2,loops,[regs.Area]);
                for k=1:length(regs)
                    bb=bb+1;
                    
                    %Determine the median value of RMS in first loop state
                    meds(bb)=median(currtrace(regs(k).PixelList(:,2)));
                    %j
                    %i
                    %median(currtrace(regs(k).PixelList(:,2)));
                end
                    %plot(trace_comp{i},'r')
                    %waitforbuttonpress
                    %close
                %    i
            end
            %find ending loops
            if sum(data.statetrace_comp{i}==2)>0
                laststep=max(find(data.onstate_comp{i}==1));
                thetrace=data.trace_comp{i};
                as=data.allstates_comp{i};
                regs2=[];
                % i
                if length(thetrace)>laststep
                    %i
                    %laststep+1
                    %length(as)
                    laststep_nextstate=max(find(as==as(laststep+1)));
                    if thetrace(laststep_nextstate)==0
                        onstate_label2=bwlabel(data.onstate_comp{i});
                        regs2 = ...
                            regionprops(onstate_label2,...
                            'PixelList','Area','Centroid','BoundingBox');
                        pccut=cat(2,pccut,regs2(end).Area);
                        cutframes=cat(2,cutframes, i);
                        %for k=1:length(regs2)
                            bbb=bbb+1;
                            meds_cut(bbb) = ...
                                median(currtrace(regs2(end).PixelList(:,2)));
                        %end
                        %i
                        %plot(trace_comp{i},'k')
                        %median(currtrace(regs2(end).PixelList(:,2)))
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
for j=1:1
subplot(3,1,j)
currloop=loops_comp{j};
currloopcut=pccut_comp{j};
currmeds=meds_comp{j};
currmedscut=meds_cut_comp{j};

[loopa,loopb]=hist(currloop(currmeds>150)./30,[50:50:2000]);
[pca,pcb]=hist(currloopcut(currmedscut>150)./30,[50:50:2000]);
%probloop=[loopa./(sum(loopa));pca./(sum(loopa))]';
%bar(loopb',probloop)
bar(loopb,loopa./(sum(loopa)),'r')
hold on
bar(pcb,pca./(sum(loopa)),'k')
xlabel('Paired complex dwell time (s)','FontSize',18);
ylabel('Probabilities','FontSize',18);
set(gca, 'xscale','linear', 'yscale','linear','ylim',[0 0.35],...
    'xlim',[0 2050],'xtick',[0 250 500 750 1000 1250 1500 1750 2000],...
    'xscale','linear','fontsize',20);
leg = legend('No bead loss','Bead loss');
set(leg,'Box','on','fontsize', 12);

annotation('textbox',[0.72 0.815 0.14 0.0345],...
    'String', {['N = ' num2str(data.n) ' beads']},...
    'FontSize',12,...
    'LineWidth',1)
figureSize = [600 600];
set(150, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
end
saveas(150, fullfile(path, specific_path,names{1}{:},'pc1200.pdf'), 'pdf');
saveas(150, fullfile(path, specific_path,names{1}{:},'pc1200.fig'), 'fig');

%%
figure(170)
mc=meds_comp{1};
mcc=meds_cut_comp{1};
% [ma,mb]=hist(mc(mc>150),[150:9:300]);
% [mca,mcb]=hist(mcc(mcc>150),[150:9:300]);
for i=1:length(mc)
    mbbp(i) = real((-0.14 + sqrt(0.14^2 - 4*(-1.92*10^-5)*(102-mc(i))))./(2*-1.92*10^-5));
end
for i=1:length(mcc)
    mcbbp(i) = real((-0.14 + sqrt(0.14^2 - 4*(-1.92*10^-5)*(102-mcc(i))))./(2*-1.92*10^-5));
end
[ma,mb]=hist(mbbp(mbbp>800),[800:100:2000]);
[mca,mcb]=hist(mcbbp(mcbbp>800),[800:100:2000]);
dlengthprob=[ma./(sum(ma)+sum(mca));mca./(sum(ma)+sum(mca))]';

hold on
%bar(mbbp',dlengthprob,'stacked')
bar(mb,ma./(sum(ma)),'r')
bar(mcb,mca./(sum(ma)),'k')
plot([1419 1419], [0 1],'b--','LineWidth',3)
xlabel('Median dwell length (bp)','FontSize',18);
ylabel('Probabilities','FontSize',18);
set(gca, 'xscale','linear','xlim',[800 1800], 'ylim',[0 0.6],...
    %{
'ytick',[0 5 10 15 20],...
%}
'yscale','linear','xscale','linear','FontSize',20);
leg = legend('No bead loss','Bead loss','Predicted Loop Length');
%set(leg,'FontSize',14);
set(leg,'Box','on','fontsize', 12);
figureSize = [600 300];
set(170, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
saveas(170, fullfile(path, specific_path,names{1}{:},'pct1200.pdf'), 'pdf');
saveas(170, fullfile(path, specific_path,names{1}{:},'pct1200.fig'), 'fig');

%%
% Find fraction of paired complex states leading to bead loss from all
% paired complex states.
frac_beadloss = sum(pca)./sum(loopa)