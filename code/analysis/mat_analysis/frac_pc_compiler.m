root_dir = '/Volumes/Soichi''s Backup/V(D)J Recombination/Code/Offline analysis--Internal/new_analysis'

all_mutations = dir(fullfile(root_dir));

mutation_name = {};
replicate = [];
frac_time_pc = [];

rowname={};

for n=1:length(all_mutations)
    
    if all_mutations(n).name(1) ~= '1'
        continue
    end
    
    names={{all_mutations(n).name}};
        
    dir_struct = dir(fullfile(root_dir, names{1}{1}));
    
    num_270s = 0;
    p = 0;
    while num_270s == 0
        p = p+1;
        num_270s = num_270s + strfind(dir_struct(p).name,'270');
    end
    
    if num_270s > 0
        analysis_file = 'analysis_270_analyzed.mat';
    else
        analysis_file = 'analysis_280_analyzed.mat';
    end

    % Adding the length of time in the paired complex state. This will be
    % summed and divided by the total time of the experiment to get the
    % fraction of time that the collection of selected beads remain in the
    % paired complex state.
    k = strfind(names{1}{1},'_');
    if length(k)==1
        rss_mutation = names{1}{1}(1:k-1);
    else
        rss_mutation = names{1}{1}(k(1)+1:k(2)-1);
    end
    
    for j=1:1
        data = load(fullfile(root_dir, names{j}{:},analysis_file),'frac_time_pc_exp');

        for m=1:length(data.frac_time_pc_exp)
            mutation_name{end+1,1} = rss_mutation;
            replicate(end+1,1)=m;
            frac_time_pc(end+1,1)=data.frac_time_pc_exp(m);
            rowname{end+1,1} = [rss_mutation num2str(m)];
        end
    
    end
    
end

T = table(mutation_name,replicate,frac_time_pc,'RowNames',rowname);

if ~exist(fullfile(root_dir,'frac_pc_data.txt'),'file')
    writetable(T,fullfile(root_dir,'frac_pc_data.txt'));
else
    Date = 10000*(year(today)-2000) + 100*month(today) + day(today);
    writetable(T,fullfile(root_dir,...
        ['frac_pc_data_', num2str(Date),'.txt']));
end
    