%function [path, nolacnames, lacnames, lacnames2,corres] = masterscriptV2params()
%
%List of the input parameters to masterscriptV2:
%path=''; %Path must end with a \ or /
%norssnames(1)=cellstr(''); %Nolac and lac sets must be in the same order
        %in the two arrays; if any  set doesn't have a nolac set, enter '' as
        %the name; nolacnames and lacnames must have the same length
%lacnames(1)=cellstr('');
%lacnames2(1)=cellstr('');%If more than one set is taken on a given area; 
        %if there are some areas with multiple sets and some without, leave
        %the string name empty for areas with only one set, i.e. lacnames
        %and lacnames2 must be the same length OR if there are no second sets
        %on any areas, leave this as lacnames2(1)=cellstr('');
%corres{1} = []; %Bead number correspondences: each element of each matrix 
        %is the number of the corresponding bead in the nolac set; must
        %have the same length as lacnames; no elements can be empty!
        
%EXAMPLE: nolac data was taken on one area with 5 beads; beads 1 and 4 were elminated when 
%lac data was taken; two lac sets were taken on this area; then lac data
%was taken on a new area for which there is no nolac data:
%nolacnames(1)=cellstr('ex_area1'); 
%nolacnames(2)=cellstr(''); %Placeholder to indicate the second lac set
        %doesn't have a corresponding nolac set
%lacnames(1)=cellstr('ex_lac_area1');
%lacnames(2)=cellstr('ex_lac_area2');
%lacnames2(1)=cellstr('ex_lac_area1_2');
%lacnames2(2)=cellstr(''); %Placeholder to make sure the second set on
        %area1 gets matched with the correct file in lacnames
%corres{1} = [2 3 5]; {[1 2 5 10 11 12 14 16 17 22]  [1]  [3 10 13]  [2 5 11 15]}
%corres{2}= [1 2 3 4];%Lacnames(1) and lacnames2(1) will use the same 
        %corres matrix; lacnames(2) has a 1-to-1 correspondance
%
%S. Johnson

function [path, nolacnames, lacnames, lacnames2,corres] = masterscriptV2params()

path = '/Volumes/Soichi Backup/Volumes/Seagate Expansion Drive/tpm_data/12CodC6A/';
%%%% CG %%%%%
%W prom
%65

 nolacnames(1)=cellstr('170810_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_NoProtein1');
 nolacnames(2)=cellstr('170810_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_NoProtein2');
 nolacnames(3)=cellstr('170810_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_NoProtein3');
 nolacnames(4)=cellstr('170815_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_NoProtein1');
 nolacnames(5)=cellstr('170815_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_NoProtein2');
 nolacnames(6)=cellstr('170816_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Protein1');
 nolacnames(7)=cellstr('170816_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_NoProtein2');
 %nolacnames(8)=cellstr('190627_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein3');
 %nolacnames(9)=cellstr('190627_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein4');
 %nolacnames(10)=cellstr('190628_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein1');
 %nolacnames(11)=cellstr('190628_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein2');
 %nolacnames(12)=cellstr('190628_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein3');
 %nolacnames(13)=cellstr('190628_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein4');
 %nolacnames(14)=cellstr('190606_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein1');
 %nolacnames(15)=cellstr('190606_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein2');
 %nolacnames(16)=cellstr('190606_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein3');
 %nolacnames(17)=cellstr('190606_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein4');
 %nolacnames(18)=cellstr('190627_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein1');
 %nolacnames(19)=cellstr('190627_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein2');
 %nolacnames(20)=cellstr('190627_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein3');
 %nolacnames(21)=cellstr('190627_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein4');
 %nolacnames(22)=cellstr('190628_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein1');
 %nolacnames(23)=cellstr('190628_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein2');
 %nolacnames(24)=cellstr('190628_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein3');
 %nolacnames(25)=cellstr('190628_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_NoProtein4');
 
 lacnames(1)=cellstr('170810_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Protein1');
 lacnames(2)=cellstr('170810_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Protein2');
 lacnames(3)=cellstr('170810_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Protein3');
 lacnames(4)=cellstr('170815_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Protein1');
 lacnames(5)=cellstr('170815_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Protein2');
 lacnames(6)=cellstr('170816_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Protein1');
 lacnames(7)=cellstr('170816_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Protein2');
 %lacnames(8)=cellstr('190627_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein3');
 %lacnames(9)=cellstr('190627_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein4');
 %lacnames(10)=cellstr('190628_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein1');
 %lacnames(11)=cellstr('190628_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein2');
 %lacnames(12)=cellstr('190628_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein3');
 %lacnames(13)=cellstr('190628_fluxus_cons1223rss_2900bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein4');
 %lacnames(14)=cellstr('190606_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein1');
 %lacnames(15)=cellstr('190606_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein2');
 %lacnames(16)=cellstr('190606_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein3');
 %lacnames(17)=cellstr('190606_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein4');
 %lacnames(18)=cellstr('190627_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein1');
 %lacnames(19)=cellstr('190627_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein2');
 %lacnames(20)=cellstr('190627_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein3');
 %lacnames(21)=cellstr('190627_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein4');
 %lacnames(22)=cellstr('190628_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein1');
 %lacnames(23)=cellstr('190628_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein2');
 %lacnames(24)=cellstr('190628_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein3');
 %lacnames(25)=cellstr('190628_cons1223rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_Ca2_Protein4');

 lacnames2(1)=cellstr('');  
 
 for i = 1:length(nolacnames)
     dataname = [nolacnames{i} '_1_POS.mat'];
     nroi = load(fullfile(path,nolacnames{i},dataname),'nroi');
     corres{i} = [1:nroi.nroi];
     nroi.nroi
 end
 
