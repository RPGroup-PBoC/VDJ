function data_comp
clear all

files = {{'Ca2+/190702_12HeptA4T_comp/190701_12HeptA4T_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca'},
    {'Ca2+/190702_12HeptA4T_comp/190702_12HeptA4T_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca'}};
     
save_path = 'Ca2+/190702_12HeptA4T_comp/';

path = '/Users/soichihirokawa/Documents/vdj_recombination/analysis/';
code_path = 'Analyzed Data/new_analysis';

mkdir(fullfile(path, code_path, save_path));

rec = [];
lacd = [];
nolacd = [];

load(fullfile(path,code_path,char(files{1}),...
    '/lacdata.mat'),'alllacnames','lacANAL_RMS','record')
load(fullfile(path,code_path,char(files{1}),...
    '/nolacdata.mat'),'nolacANAL_RMS')

names = alllacnames;
rec = record;
lacd = lacANAL_RMS;
nolacd = nolacANAL_RMS;


if length(files) > 1
    for j = 2:length(files)
        load(fullfile(path,code_path,char(files{j}),...
            '/lacdata.mat'),'alllacnames','lacANAL_RMS','record');
        load(fullfile(path,code_path,char(files{j}),...
            '/nolacdata.mat'),'nolacANAL_RMS');
        
        for i = 1:length(lacANAL_RMS);
            %if ~(j==2 & i==11);
                alllacnames(i)
                names{1,end+1} = alllacnames{1,i};
                rec{1,end+1} = record{1,i};
                lacd{1,end+1} = lacANAL_RMS{1,i};
                nolacd{1,end+1} = nolacANAL_RMS{1,i};
            
        end
    end
    
end

clear alllacnames lacANAL_RMS nolacANAL_RMS record

alllacnames = names;
lacANAL_RMS = lacd;
nolacANAL_RMS = nolacd;
record = rec;

save(fullfile(path,code_path,save_path,'/lacdata.mat'),'alllacnames','lacANAL_RMS','record');
save(fullfile(path,code_path,save_path,'/nolacdata.mat'),'nolacANAL_RMS','record');

end
