path = '/Users/soichihirokawa/Documents/vdj_recombination/analysis';
specific_path = '/Analyzed Data/new_analysis/coding_flank/';

all_mutations = dir(fullfile(path,specific_path));

for n=1:length(all_mutations)
    
    if all_mutations(n).name(1) ~= '1' | all_mutations(n).name(end-3:end) == '.pdf'
        continue
    end
    
    underscore_index = findstr(all_mutations(n).name, '_');
    mutation_name = all_mutations(n).name(underscore_index(1)+1:underscore_index(2)-1)
    
    names={{all_mutations(n).name}};

    % For naming some of our figures later on, we are going to extract the
    % RSS mutation name.
    
    mutation_files = dir(fullfile(path, specific_path, all_mutations(n).name));
        
    for m=1:length(mutation_files)
        if contains(mutation_files(m).name,'_275_analyzed.mat')
            analysis_file = mutation_files(m).name
        end
    end

    loops_comp=[];
    pccut_comp=[];
    totallength_comp=[];
    meds_comp = [];
    meds_cut_comp = [];
    

    % Adding the length of time in the paired complex state. This will be
    % summed and divided by the total time of the experiment to get the
    % fraction of time that the collection of selected beads remain in the
    % paired complex state.

    for j=1:1
        data = load(...
            fullfile(path, specific_path, names{j}{:},analysis_file));
        loops=[];
        pccut=[];
        cutframes=[];
        totallength=[];
        meds=[];
        meds_cut=[];

        loops_comp_tmp = [];
        pccut_comp_tmp = [];
        meds_comp_tmp = [];
        meds_cut_comp_tmp = [];
        totallength_comp_tmp = [];
        
        for i=1:length(data.statetrace_comp)
            bb=0;
            bbb=0;
            
            totallength{i} = 0;
            loops{i} = 0;
            pccut{i} = 0;
            cutframes{i} = 0;
            meds{i} = 0;
            meds_cut{i} = 0;
            for l=1:size(data.statetrace_comp{i},1)
                currtrace=data.trace_comp{i}(l,:);
                
                if sum(currtrace(1:1.5*10^4)>290)>0
                    totallength{i}=totallength{i}+length(currtrace);
                    %i
                    %find all loops
                    regs=[];
                    if sum(data.statetrace_comp{i}(l,:)==2)>0

                        onstate_label=bwlabel(transpose(data.onstate_comp{i}(l,:)));
                        % PixelList gives every pixel in every looped state
                        regs=regionprops(onstate_label,'PixelList',...
                            'Area','Centroid','BoundingBox');
                        [regs.Area]./30;
                        loops{i}=cat(2,loops{i},[regs.Area]);
                        for k=1:length(regs)
                            bb=bb+1;

                            %Determine the median value of RMS in first loop state
                            meds{i}(bb)=median(currtrace(regs(k).PixelList(:,2)));
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
                    if sum(data.statetrace_comp{i}(l,:)==2)>0
                        laststep=max(find(transpose(data.onstate_comp{i}(l,:)==1)));
                        thetrace=data.trace_comp{i}(l,:);
                        as=data.allstates_comp{i}(l,:);
                        regs2=[];
                        % i
                        if length(thetrace)>laststep
                            %i
                            %laststep+1
                            %length(as)
                            laststep_nextstate=max(find(as==as(laststep+1)));
                            if thetrace(laststep_nextstate)==0
                                onstate_label2=bwlabel(transpose(data.onstate_comp{i}(l,:)));
                                regs2 = ...
                                    regionprops(onstate_label2,...
                                    'PixelList','Area','Centroid','BoundingBox');
                                pccut{i}=cat(2,pccut{i},regs2(end).Area);
                                cutframes{i}=cat(2,cutframes{i}, i);
                                %for k=1:length(regs2)
                                    bbb=bbb+1;
                                    meds_cut{i}(bbb) = ...
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

        end
        for k=1:length(loops)
            loops_comp_tmp=cat(2,loops_comp_tmp,loops{k});
            pccut_comp_tmp=cat(2,pccut_comp_tmp,pccut{k});
            totallength_comp_tmp=cat(2,totallength_comp_tmp,totallength{k});
            meds_comp_tmp=cat(2,meds_comp_tmp,meds{k});
            meds_cut_comp_tmp=cat(2,meds_cut_comp_tmp,meds_cut{k});
        end
        
        loops_comp{j} = loops_comp_tmp(loops_comp_tmp~=0);
        pccut_comp{j} = pccut_comp_tmp(pccut_comp_tmp~=0);
        totallength_comp{j} = totallength_comp_tmp(totallength_comp_tmp~=0);
        meds_comp{j} = meds_comp_tmp(meds_comp_tmp~=0);
        meds_cut_comp{j} = meds_cut_comp_tmp(meds_cut_comp_tmp~=0);
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

    annotation('textbox',[0.74 0.80 0.15 0.05],...
        'String', {['N = ' num2str(data.n) ' beads'],...
        ['n_{PC} = ' num2str(sum(loopa)) ' paired complexes']},...
        'FontSize',12,'LineWidth',1,'FitBoxToText','on',...
        'HorizontalAlignment','right','VerticalAlignment','middle')
    figureSize = [600 600];
    set(150, 'Position', [0 0 figureSize+[100 50]],'PaperUnits','points','PaperSize',figureSize,'PaperPosition',[0 0 figureSize]);
    end
    saveas(150, fullfile(path, specific_path,names{1}{:},['pc_' mutation_name '.pdf']), 'pdf');
    saveas(150, fullfile(path, specific_path,names{1}{:},['pc_' mutation_name '.fig']), 'fig');

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
    if ~isempty(mcc)
        [mca,mcb]=hist(mcbbp(mcbbp>800),[800:100:2000]);
    else
        mca=0;
    end

    hold on
    %bar(mbbp',dlengthprob,'stacked')
    bar(mb,ma./(sum(ma)),'r')
    if ~isempty(mcc)
        bar(mcb,mca./(sum(ma)),'k')
    end
    
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
    saveas(170, fullfile(path, specific_path,names{1}{:},['pct_' mutation_name '.pdf']), 'pdf');
    saveas(170, fullfile(path, specific_path,names{1}{:},['pct_' mutation_name '.fig']), 'fig');

    %%
    filename = fullfile(path, specific_path, names{1}{:}, analysis_file);
    save(filename, 'loops','pccut','meds','meds_cut','-append')

end