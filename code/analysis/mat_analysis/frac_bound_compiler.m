root_dir = '/Volumes/Soichi''s Backup/V(D)J Recombination/Code/Offline analysis--Internal/Analyzed Data/single_binding/';
mutation_name = 'REF12RSS';
data_struct = 'GaussFit.mat';

% Find concentration of RAG for each experiment
all_concentrations = dir(fullfile(root_dir,mutation_name));

concentration = [];
frac_bound = [];

for n = 1:length(all_concentrations)
    concentration_val = 0;
    
    if all_concentrations(n).name(1) ~= '1'
        continue
    end
    
    concentration_name = {{all_concentrations(n).name}};
    
    underscore_counts = strfind(concentration_name{1}{1},'_');
    concentration_val = str2double(concentration_name{1}{1}(underscore_counts(end)+1:end-10));
    
    prob_values = load(fullfile(root_dir,mutation_name,concentration_name{1}{1},...
        'Single bead analysis',data_struct));
    
    for m=1:length(prob_values.GaussFit)
        concentration(end+1,1) = concentration_val;
        frac_bound(end+1,1) = prob_values.GaussFit(m).p(1);
    end
    
end

T = table(concentration,frac_bound);

if ~exist(fullfile(root_dir,mutation_name,[mutation_name '_' 'individual_frac_bound.txt']),'file')
    writetable(T,fullfile(root_dir,mutation_name,[mutation_name '_' 'individual_frac_bound.txt']));
else
    Date = 10000*(year(today)-2000) + 100*month(today) + day(today);
    writetable(T,fullfile(root_dir,mutation_name,...
        [mutation_name '_' 'individual_frac_bound_', num2str(Date),'.txt']));
end
    
