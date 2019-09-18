path = '/Users/soichihirokawa/Documents/V(D)J Recombination/Code';
specific_path = '/Offline analysis--Internal/Analyzed Data';
names={{'consensus_comp'}};
loops_comp=[];
pccut_comp=[];
totallength_comp=[];
for j=1:1
    data = load(...
        fullfile(path, specific_path, names{j}{:},'analysis_285.mat'));
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
                    
                    meds(bb)=median(currtrace(regs(k).PixelList(:,2)));
                end
            end
            %find ending loops
            if sum(data.statetrace_comp{i}==2)>0
                laststep=max(find(data.onstate_comp{i}==1));
                thetrace=data.trace_comp{i};
                as=data.allstates_comp{i};
                regs2=[];
                if length(thetrace)>laststep
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
                    end
                end
            end
        end
    end
    % Compiles the dwell time in the paired complex state
    loops_comp{j}=loops;
    % Compiles the dwell time in the paired complex state prior to cleavage
    pccut_comp{j}=pccut;
    totallength_comp{j}=totallength;
    meds_comp{j}=meds;
    meds_cut_comp{j}=meds_cut;
end

% Place values into posterior
total_time = sum(pccut./30);
t_max = 3600;
n = length(pccut);

tau = linspace(0.01, 1000, 10000);
log_post = -total_time./ tau - n * log(tau) - n * log(1 - exp(-t_max./tau));
post = exp(log_post - max(log_post));
post = post / trapz(tau, post);

plot(tau, post)
xlabel('\tau (sec)', 'FontSize', 16)
ylabel('P(\tau | \bf t)', 'FontSize', 16)

% Find maximum likelihood value of tau
[M, I] = max(post);
readout = ['The maximum likelihood value of tau is ', num2str(tau(I))];
disp(readout)