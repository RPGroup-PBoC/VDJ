path = '/Users/soichihirokawa/Documents/V(D)J Recombination/Code';
specific_path = '/Offline analysis--Internal/Analyzed Data/new_analysis';

% Create tables of mutations and resultant rates and probabilities
mutation_name = {};

mean_bead_pcf = [];
std_bead_pcf = [];
pcf_count_comp = [];

mean_frac_pc = [];
std_frac_pc = [];

pop_cleaved = [];
num_loops = [];

mean_cleaved = [];
std_cleaved = [];

cut_rate = [];
err_cut_rate = [];

unloop_rate = [];
err_unloop_rate = [];

loop_rate = [];
err_loop_rate = [];

true_frac_pc = [];
std_true_frac_pc = [];

all_mutations = dir(fullfile(path,specific_path));

for n=1:length(all_mutations)
    
    if all_mutations(n).name(1) ~= '1'
        continue
    end
    
    names={{all_mutations(n).name}};

    % For naming some of our figures later on, we are going to extract the
    % RSS mutation name.
    
    mutation_filenames = dir(fullfile(path,specific_path,all_mutations(n).name));
    
    for findmatfile=1:length(mutation_filenames)
        if length(mutation_filenames(findmatfile).name) < 12
            continue            
        elseif strcmp(mutation_filenames(findmatfile).name(end-11:end),'analyzed.mat')
            analysis_file = mutation_filenames(findmatfile).name;
        end
    end

    % Adding the length of time in the paired complex state. This will be
    % summed and divided by the total time of the experiment to get the
    % fraction of time that the collection of selected beads remain in the
    % paired complex state.
    k = strfind(names{1}{1},'_');
    if length(k)==1
        mutation_name{end+1,1} = names{1}{1}(1:k-1);
    else
        mutation_name{end+1,1} = names{1}{1}(k(1)+1:k(2)-1);
    end

    data = load(fullfile(path, specific_path, names{1}{:},analysis_file));
    
    % Determine the probability of paired complex formation, the unlooping 
    % rate, the on rate of paired complex formation and the cutting rate.
    fate_rate = 0;
    std_fate_rate = 0;
    
    % Find the fraction of beads that go into the paired complex state
    pcf_frac = [];
    pcf_count = 0;
    for l = 1:length(data.ontime_comp)
        num_pcf = 0;
        slice_beads = data.ontime_comp{l}(:,1);
        num_pcf = slice_beads(slice_beads~=0);
        pcf_count = pcf_count + size(num_pcf,1);
        pcf_frac(l) = size(num_pcf,1)./size(data.ontime_comp{l},1);
    end
    
    mean_bead_pcf(end+1,1) = mean(pcf_frac);
    std_bead_pcf(end+1,1) = nanstd(pcf_frac);
    pcf_count_comp(end+1,1) = pcf_count;
    
    % Find fraction of time spent in paired complex state for each 
    % experiment
    tot_time_pc_bead = [];
    frac_time_pc_exp = [];
    count = 0;
    for l = 1:length(data.ontime_comp)
        if strcmp(mutation_name{end,1},'WT1223RSS') && l==6
            continue
        end
        count = count + 1;
        tot_time_pc_bead{count} = sum(data.ontime_comp{l},2)./30;
        frac_time_pc_exp(count) = sum(tot_time_pc_bead{end})./sum(data.time_comp{l});
    end
    
    % Determine the mean and standard deviation for the fraction of time
    % spent in the paired complex state
    mean_frac_pc(end+1,1) = round(nanmean(frac_time_pc_exp),4);
    std_frac_pc(end+1,1) = round(nanstd(frac_time_pc_exp,1),4);
    
    % Determine rates of looping and paired complex fate, which considers
    % both unlooping and cleavage
    
    % ontime_exp is the mean time in the paired complex state within a
    % replicate
    ontime_rep = [];
    ontime_tot = [];
    for l=1:length(data.ontime_comp)
        ontime_rep = cat(1,ontime_rep,data.ontime_comp{l}(:));
        ontime_tot = ontime_rep(ontime_rep~=0);
    end
    fate_rate = 30./nanmean(ontime_tot);
    std_fate_rate = 30.*nanstd(ontime_tot)/(nanmean(ontime_tot)^2);
    
    % Find fraction of paired complex states leading to bead loss from all
    % paired complex states, the fraction of time that the collection of beads
    % are in the paired complex state and the fraction of time of paired
    % complexes that led to cleavage over the total time of paired complex
    % states.
    loops = [];
    meds = [];
    pccut = [];
    meds_cut = [];
    loops_comp = [];
    meds_comp = [];
    pccut_comp = [];
    meds_cut_comp = [];
    
    num_loop = [];
    num_cleave = [];
    frac_cleaved = [];
    
    for l=1:length(data.loops)
        loops{l} = data.loops{l}(data.loops{l} ~= 0);
        loops_comp = cat(2, loops_comp, loops{l});
        
        meds{l} = data.meds{l}(data.meds{l} ~= 0);
        meds_comp = cat(2, meds_comp, meds{l});
        
        pccut{l} = data.pccut{l}(data.pccut{l} ~= 0);
        pccut_comp = cat(2, pccut_comp, pccut{l});
        
        meds_cut{l} = data.meds_cut{l}(data.meds_cut{l} ~= 0);
        meds_cut_comp = cat(2, meds_cut_comp, meds_cut{l});
    end
    [num_loop,~] = hist(loops_comp(meds_comp>150)./30,...
        [50:50:2000]);
    [num_cleave,~] = hist(pccut_comp(meds_cut_comp>150)./30,...
        [50:50:2000]);
    pop_cleaved(end+1,1) = sum(num_cleave,2)./sum(num_loop,2);
    num_loops(end+1,1) = sum(num_loop,2);

    frac_cleaved_replicate = [];
    
    for l=1:length(data.loops)
        
        loops_rep = [];
        meds_rep = [];
        pccut_rep = [];
        meds_cut_rep = [];
    
        num_loop_replicate = [];
        num_cleave_replicate = [];
        
        loops_rep = data.loops{l}(data.loops{l} ~= 0);
        meds_rep = data.meds{l}(data.meds{l} ~= 0);
        pccut_rep = data.pccut{l}(data.pccut{l} ~= 0);
        meds_cut_rep = data.meds_cut{l}(data.meds_cut{l} ~= 0);
        
        [num_loop_replicate,~] = hist(loops_rep(meds_rep>150)./30,...
            [50:50:2000]);
        [num_cleave_replicate,~] = hist(pccut_rep(meds_cut_rep>150)./30,...
            [50:50:2000]);
        frac_cleaved_replicate(l) = sum(num_cleave_replicate)./sum(num_loop_replicate);
    end
    mean_cleaved(end+1,1) = nanmean(frac_cleaved_replicate);
    std_cleaved(end+1,1) = nanstd(frac_cleaved_replicate);
        
    
    % Find rates of unlooping and cutting using the fate rate and the
    % fraction of paired complexes cleaved   
    unloop_rate(end+1,1) = fate_rate * (1 - pop_cleaved(end,1));
    err_unloop_rate(end+1,1) = sqrt(((1 - pop_cleaved(end,1))^2*std_fate_rate)^2);
    cut_rate(end+1,1) = fate_rate * pop_cleaved(end,1);
    err_cut_rate(end+1,1) = sqrt((pop_cleaved(end,1) * std_fate_rate)^2);
    
    % Find rate looping using the rate of cutting, unlooping and fraction
    % of time spent in the paired complex state
    loop_rate(end+1,1) = (cut_rate(end,1)+unloop_rate(end,1)) * ...
        (mean_frac_pc(end,1))/(1 - mean_frac_pc(end,1));
    err_loop_rate(end+1,1) = sqrt(((1 - mean_frac_pc(end,1))/mean_frac_pc(end,1))^2*...
        (err_unloop_rate(end,1)^2 + err_cut_rate(end,1)^2)+...
        ((unloop_rate(end,1) + cut_rate(end,1))*std_frac_pc(end,1)/(mean_frac_pc(end,1))^2)^2);

    if pop_cleaved(end,1) == 1
        true_frac_pc(end+1,1) = mean_frac_pc(end,1);
        std_true_frac_pc(end+1,1) = std_frac_pc(end,1);
    else
        true_frac_pc(end+1,1) = loop_rate(end,1)/(unloop_rate(end,1) + loop_rate(end,1));
        std_true_frac_pc(end+1,1) = sqrt((std_frac_pc(end,1)*(1 - pop_cleaved(end,1)))^2/...
            (1 - pop_cleaved(end,1) + pop_cleaved(end,1)*mean_frac_pc(end,1))^4);
    end
    
    filename = fullfile(path, specific_path, names{1}{:}, analysis_file);
    save(filename,'tot_time_pc_bead','frac_time_pc_exp','-append');  
end

T = table(mutation_name,mean_bead_pcf,std_bead_pcf,pcf_count_comp,...
    mean_frac_pc,std_frac_pc,pop_cleaved,num_loops,mean_cleaved,std_cleaved,cut_rate,...
    err_cut_rate,unloop_rate,err_unloop_rate,loop_rate,err_loop_rate,...
    true_frac_pc,std_true_frac_pc,'RowNames',mutation_name);

if ~exist(fullfile(path, specific_path,'all_results.txt'),'file')
    writetable(T,fullfile(path,specific_path,'all_results.txt'));
else
    Date = 10000*(year(today)-2000) + 100*month(today) + day(today);
    writetable(T,fullfile(path,specific_path,...
        ['all_results_', num2str(Date),'.txt']));
end
%}