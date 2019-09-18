%[conc,data] = LoadDataAnalysis(DNAset,kind,varargin)
%
%Loads the most current data I have on:
%
%"DNAset": this is the name of the subfolder in SJLacI or SJLac2; for
%   example, "E894 conc curve".  The exceptions are: for WTLac or a 
%   no-operator control ("noOidE872O1" or "noOidE872noO1"), need the name
%   of the subfolder in WTLac conc curves or in Operator KO controls; 
%   for linked channels, need "linked channels"; for small beads, 
%   "Indicia Beads"; for HG seqs wprom, 'HGseqsWProm'; 
%   for E898 conc curve with small beads, "E898conccurve".
%"Kind": this can be "stats","GaussFit","lacdataconcat", "HistoAnal" (MUST
%   be followed by one, optionally two, additional inputs that select only part
%   of that file; see below), "names" which returns a cell array of paths
%   to the current data analysis folders for a DNAset.
%Optional inputs:
%--"SJLacI","SJLacI2","Old","New" (indicates the repressor
%   batch).  Default is "SJLacI".
%--analmethod: if "Kind" is "HistoAnal", pass a repressor batch then
%   "Means"/"Meanssub0"/"GaussMeans"/"WeightMeans"/"WeightMeanssub0"/
%   "alldistribs"/"alldistribssub0".  THIS MUST BE PASSED IF KIND = HISTOANAL
%--numsecs: if "kind" is "HistoAnal", pass a repressor batch, then an analmethod,
%   then 1000/2000/3000/4000/5000 to select data that only lasted 1000, ..., 5000
%   seconds.
%
%Outputs are "conc" (which is either the appropriate set of concentrations
%or loop lengths) and "data' which is a cell of GaussFit's, stat's, etc,
%one for each concentration or length.  (If "kind" is HistoAnal and
%analmethod is anything but alldistribs/_sub0 then data is a 10x2 cell,
%with the second set of elements being the SEs.)
%
%Steph 9/10

function [conc,data] = LoadDataAnalysis(DNAset,kind,varargin)

%%%%%%% LIST OF GOOD DATA %%%%%%
%Concentration curves
%O1-E894-O1 conc curve
concE8O1O1 = 10^-12.*[1 2.5 5 10 50 100 500 1000];
E8O1O1conc{1} = 'O1E894O1_1pMSJLacI';
E8O1O1conc{2} = 'O1E894O1_2500fMSJLacI_moredata';
E8O1O1conc{3} = 'O1E894O1_5pMSJLacI';
E8O1O1conc{4} = 'O1E894O1_10pMSJLacI';
%E8O1O1conc{5} = 'O1E894O1_50pMSJLacI';
%E8O1O1conc{5} = 'O1E894O1_50pMSJLacI_redo';
E8O1O1conc{5} = 'O1E894O1_50pMSJLacI_redo2';
E8O1O1conc{6} = 'O1E894O1_100pMSJLacI';
E8O1O1conc{7} = 'O1E894O1_500pMSJLacI';
E8O1O1conc{8} = 'O1E894O1_1nMSJLacI';

%O2-E894-O1 conc curve
concE8O2O1 = 10^-12.*[1 5 10 50 100 500];
E8O2O1conc{1} = 'O2E894O1_1pMSJLacI';
E8O2O1conc{2} = 'O2E894O1_5pMSJLacI';
E8O2O1conc{3} = 'O2E894O1_10pMSJLacI';
E8O2O1conc{4} = 'O2E894O1_50pMSJLacI';
%E8O2O1conc{5} = 'O2E894O1_100pMSJLacI';
E8O2O1conc{5} = 'O2E894O1_100pMSJLacI_moredata';
%E8O2O1conc{6} = 'O2E894O1_500pMSJLacI';
E8O2O1conc{6} = 'O2E894O1_500pMSJLacI_moredata';

%Oid-E894-O1 conc curve
concE8 = 10^-12.*[0.5 1 2.5 5 10 25 50 100 500 10^3 10^4];
E8conc{1} = 'E894_500fMSJLacI_moredata';
E8conc{2} = 'E894_1pMSJLacI_moredata';
E8conc{3} = 'E894_2500fMSJLacI';
%E8conc{4} = 'E894_5pMSJLacI';
E8conc{4} = 'E894_5pMSJLacI_moredata2';
E8conc{5} = 'E894_10pMSJLacI';
E8conc{6} = 'E894_25pMSJLacI_redo_moredata';
E8conc{7} = 'E894_50pMSJLacI_doublebds';
E8conc{8} = 'E894_100pMSJLacI_moredata';
E8conc{9} = 'E894_500pMSJLacI';
E8conc{10} = 'E894_1nMSJLacI_moredata';
E8conc{11} = 'E894_10nMSJLacI';

% concE8 = 10^-12.*[1 2.5 5 10 25 50 100 500 10^3 10^4];
% E8conc{1} = 'E894_1pMSJLacI_moredata';
% E8conc{2} = 'E894_2500fMSJLacI';
% %E8conc{3} = 'E894_5pMSJLacI';
% E8conc{3} = 'E894_5pMSJLacI_moredata2';
% E8conc{4} = 'E894_10pMSJLacI';
% E8conc{5} = 'E894_25pMSJLacI_redo_moredata';
% E8conc{6} = 'E894_50pMSJLacI_doublebds';
% E8conc{7} = 'E894_100pMSJLacI_moredata';
% E8conc{8} = 'E894_500pMSJLacI';
% E8conc{9} = 'E894_1nMSJLacI_moredata';
% E8conc{10} = 'E894_10nMSJLacI';

%Oid-TA94-O1 conc curve
concTA = 10^-12.*[0.1 0.25 0.5 1 2.5 5 10 50 100 500 10^3 10^4 10^5];
TAconc{1} = 'TA94_100fMSJLacI';
TAconc{2} = 'TA94_250fMSJLacI';
TAconc{3} = 'TA94_500fMSJLacI';
TAconc{4} = 'TA94_1pMSJLacI';
TAconc{5} = 'TA94_2500fMSJLacI';
%TAconc{} = 'TA94_5pMSJLacI_moredata';
TAconc{6} = 'TA94_5pMSJLacI_moredata2';
TAconc{7} = 'TA94_10pMSJLacI_moredata';
TAconc{8} = 'TA94_50pMSJLacI_moredata';
%TAconc{} = 'TA94_100pMSJLacI_redo';
TAconc{9} = 'TA94_100pMSJLacI_redo_moredata';
%TAconc{} = 'TA94_500pMSJLacI';
TAconc{10} = 'TA94_500pMSJLacI_moredata';
%TAconc{} = 'TA94_1nMSJLacI_moredata';
TAconc{11} = 'TA94_1nMSJLacI_moredata2';
%TAconc{} = 'TA94_10nMSJLacI';
%TAconc{} = 'TA94_10nMSJLacI redo'; %Which one of these did I want to use???
TAconc{12} = 'TA94_10nMSJLacI_redo_moredata'; %decided to go with the redo ...
TAconc{13} = 'TA94_100nMSJLacI';

% concTA = 10^-12.*[1 2.5 5 10 50 100 500 10^3 10^4 10^5];
% TAconc{1} = 'TA94_1pMSJLacI';
% TAconc{2} = 'TA94_2500fMSJLacI';
% %TAconc{} = 'TA94_5pMSJLacI_moredata';
% TAconc{3} = 'TA94_5pMSJLacI_moredata2';
% TAconc{4} = 'TA94_10pMSJLacI_moredata';
% TAconc{5} = 'TA94_50pMSJLacI_moredata';
% %TAconc{} = 'TA94_100pMSJLacI_redo';
% TAconc{6} = 'TA94_100pMSJLacI_redo_moredata';
% %TAconc{} = 'TA94_500pMSJLacI';
% TAconc{7} = 'TA94_500pMSJLacI_moredata';
% %TAconc{} = 'TA94_1nMSJLacI_moredata';
% TAconc{8} = 'TA94_1nMSJLacI_moredata2';
% %TAconc{} = 'TA94_10nMSJLacI';
% %TAconc{} = 'TA94_10nMSJLacI redo'; %Which one of these did I want to use???
% TAconc{9} = 'TA94_10nMSJLacI_redo_moredata'; %decided to go with the redo ...
% TAconc{10} = 'TA94_100nMSJLacI';

%O1-PUC306-Oid conc curve
concPUC = 10^-12.*[0.01 0.1 0.5 1 2.5 5 10 50 100 500 10^3 10^4];
PUCconc{1} = 'PUC306_10fMSJLacI';
PUCconc{2} = 'PUC306_100fMSJLacI';
%PUCconc{3} = 'PUC306_500fMSJLacI';
PUCconc{3} = 'PUC306_500fMSJLacI_moredata';
%PUCconc{4} = 'PUC306_1pMSJLacI';
PUCconc{4} = 'PUC306_1pMSJLacI_moredata';
PUCconc{5} = 'PUC306_2500fMSJLacI';
%PUCconc{6} = 'PUC306_5pMSJLacI';
PUCconc{6} = 'PUC306_5pMSJLacI_moredata';
PUCconc{7} = 'PUC306_10pMSJLacI';
%PUCconc{8} = 'PUC306_50pMSJLacI';
PUCconc{8} = 'PUC306_50pMSJLacI_moredata';
PUCconc{9} = 'PUC306_100pMSJLacI';
PUCconc{10} = 'PUC306_500pMSJLacI';
PUCconc{11} = 'PUC306_1nMSJLacI';
%PUCconc{12} = 'PUC306_10nMSJLacI';
PUCconc{12} = 'PUC306_10nMSJLacI_moredata';

% concPUC = 10^-12.*[1 2.5 5 10 50 100 500 10^3 10^4];
% PUCconc{1} = 'PUC306_1pMSJLacI_moredata';
% PUCconc{2} = 'PUC306_2500fMSJLacI';
% %PUCconc{} = 'PUC306_5pMSJLacI';
% PUCconc{3} = 'PUC306_5pMSJLacI_moredata';
% PUCconc{4} = 'PUC306_10pMSJLacI';
% %PUCconc{} = 'PUC306_50pMSJLacI';
% PUCconc{5} = 'PUC306_50pMSJLacI_moredata';
% PUCconc{6} = 'PUC306_100pMSJLacI';
% PUCconc{7} = 'PUC306_500pMSJLacI';
% PUCconc{8} = 'PUC306_1nMSJLacI';
% %PUCconc{} = 'PUC306_10nMSJLacI';
% PUCconc{9} = 'PUC306_10nMSJLacI_moredata';

%O1-PUC306-Oid DNA check
concPUCcheck = 10^-12.*[1 10 100];
PUCcheck{1} = 'PUC306_newDNA_1pMSJLacI';
PUCcheck{2} = 'PUC306_newDNA_10pMSJLacI';
PUCcheck{3} = 'PUC306_newDNA_100pMSJLacI';

%Oid-E898-O1 conc curve
concE898 = 10^-12.*[1 5 25 100 1000];
E898conc{1} = 'OidE898O1_1pMSJLacI_smallbds_fBpt07fGpt0461';
E898conc{2} = 'OidE898O1_5pMSJLacI_smallbds_fBpt07fGpt0461';
E898conc{3} = 'OidE898O1_25pMSJLacI_smallbds_fBpt07fGpt0461';
E898conc{4} = 'OidE898O1_100pMSJLacI_smallbds_fBpt07fGpt0461';
E898conc{5} = 'OidE898O1_1nMSJLacI_smallbds_fBpt07fGpt0461';

%Oid-E8107-O1 conc curve
concE8107 = 10^-12.*[1 5 25 100 1000];
E8107conc{1} = 'OidE8107O1_1pMSJLacI';
E8107conc{2} = 'OidE8107O1_5pMSJLacI';
E8107conc{3} = 'OidE8107O1_25pMSJLacI';
E8107conc{4} = 'OidE8107O1_100pMSJLacI';
E8107conc{5} = 'OidE8107O1_1nMSJLacI';

%Oid-E8108-O1 conc curve
concE8108 = 10^-12.*[1 5 25 100 1000];
E8108conc{1} = 'OidE8108O1_1pMSJLacI';
E8108conc{2} = 'OidE8108O1_5pMSJLacI';
E8108conc{3} = 'OidE8108O1_25pMSJLacI';
E8108conc{4} = 'OidE8108O1_100pMSJLacI';
E8108conc{5} = 'OidE8108O1_1nMSJLacI';

%Single- and no-operator:
%O1 only, with E72
concnoOidE872O1 = 10^-12.*[1000 10000];
noOidE872O1conc{1} = 'noOidE872O1_1nMSJLacI';
noOidE872O1conc{2} = 'noOidE872O1_10nMSJLacI';

%no operators, E72
concnoOidE872noO1 = 10^-12.*[1000 10000];
noOidE872noO1conc{1} = 'noOidE872noO1_1nMSJLacI';
noOidE872noO1conc{2} = 'noOidE872noO1_10nMSJLacI';

%%%% LENGTH SERIES: default is to load 100 pM %%%%
%lengthsE8 = 83:100;
%lengthsTA = 83:100;
%lengthsE8 = [83 84 89 90 91 92 93 94 95 96 97 98 99 100];
%lengthsTA = [83 84 89 90 91 92 93 94 95 96 97 98 99 100];
%lengthsE8 = 89:100;
%lengthsTA = 89:100;
lengthsE8 = [89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 116];
lengthsTA = [89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 116];

%This is 50 pM
% E8length{1} = fullfile('Oid-E,T 89 to 100-O1','OidE889O1_50pMSJLacI');
% E8length{2} = fullfile('Oid-E,T 89 to 100-O1','OidE892O1_50pMSJLacI');
% E8length{3} = fullfile('E894 conc curve',E8conc{logical(concE8==50*10^-12)}); %This is from the Oid-O1 series, above
% E8length{4} = fullfile('Oid-E,T 89 to 100-O1','OidE895O1_50pMSJLacI');
% E8length{5} = fullfile('Oid-E,T 89 to 100-O1','OidE896O1_50pMSJLacI');
% E8length{6} = fullfile('Oid-E,T 89 to 100-O1','OidE898O1_50pMSJLacI');
% E8length{7} = fullfile('Oid-E,T 89 to 100-O1','OidE8100O1_50pMSJLacI');

%This is 100 pM
%E8length{1} = fullfile('Oid-E,T 89 to 100-O1','OidE883O1_100pMSJLacI');
%E8length{2} = fullfile('Oid-E,T 89 to 100-O1','OidE884O1_100pMSJLacI');
E8length{1} = fullfile('Oid-E,T 89 to 100-O1','OidE889O1_100pMSJLacI');
E8length{2} = fullfile('Oid-E,T 89 to 100-O1','OidE890O1_100pMSJLacI');
E8length{3} = fullfile('Oid-E,T 89 to 100-O1','OidE891O1_100pMSJLacI');
E8length{4} = fullfile('Oid-E,T 89 to 100-O1','OidE892O1_100pMSJLacI');
E8length{5} = fullfile('Oid-E,T 89 to 100-O1','OidE893O1_100pMSJLacI');
E8length{6} = fullfile('E894 conc curve',E8conc{logical(concE8==100*10^-12)}); %This is from the Oid-O1 series, above
E8length{7} = fullfile('Oid-E,T 89 to 100-O1','OidE895O1_100pMSJLacI');
E8length{8} = fullfile('Oid-E,T 89 to 100-O1','OidE896O1_100pMSJLacI');
E8length{9} = fullfile('Oid-E,T 89 to 100-O1','OidE897O1_100pMSJLacI');
E8length{10} = fullfile('Oid-E,T 89 to 100-O1','OidE898O1_100pMSJLacI');
E8length{11} = fullfile('Oid-E,T 89 to 100-O1','OidE899O1_100pMSJLacI');
E8length{12} = fullfile('Oid-E,T 89 to 100-O1','OidE8100O1_100pMSJLacI');
E8length{13} = fullfile('Oid-E,T 89 to 100-O1','OidE8101O1_100pMSJLacI');
E8length{14} = fullfile('Oid-E,T 89 to 100-O1','OidE8102O1_100pMSJLacI');
E8length{15} = fullfile('Oid-E,T 89 to 100-O1','OidE8103O1_100pMSJLacI');
E8length{16} = fullfile('Oid-E,T 89 to 100-O1','OidE8104O1_100pMSJLacI');
E8length{17} = fullfile('Oid-E,T 89 to 100-O1','OidE8105O1_100pMSJLacI');
E8length{18} = fullfile('Oid-E,T 89 to 100-O1','OidE8106O1_100pMSJLacI');
E8length{19} = fullfile('Oid-E,T 89 to 100-O1','OidE8107O1_100pMSJLacI');
E8length{20} = fullfile('Oid-E,T 89 to 100-O1','OidE8108O1_100pMSJLacI');
E8length{21} = fullfile('Oid-E,T 89 to 100-O1','OidE8109O1_100pMSJLacI');
E8length{22} = fullfile('Oid-E,T 89 to 100-O1','OidE8116O1_100pMSJLacI');

%TAlength{1} = fullfile('Oid-E,T 89 to 100-O1','OidTA83O1_100pMSJLacI');
%TAlength{2} = fullfile('Oid-E,T 89 to 100-O1','OidTA84O1_100pMSJLacI');
TAlength{1} = fullfile('Oid-E,T 89 to 100-O1','OidTA89O1_100pMSJLacI');
TAlength{2} = fullfile('Oid-E,T 89 to 100-O1','OidTA90O1_100pMSJLacI');
TAlength{3} = fullfile('Oid-E,T 89 to 100-O1','OidTA91O1_100pMSJLacI');
TAlength{4} = fullfile('Oid-E,T 89 to 100-O1','OidTA92O1_100pMSJLacI');
TAlength{5} = fullfile('Oid-E,T 89 to 100-O1','OidTA93O1_100pMSJLacI');
TAlength{6} = fullfile('TA94 conc curve',TAconc{logical(concTA==100*10^-12)});
TAlength{7} = fullfile('Oid-E,T 89 to 100-O1','OidTA95O1_100pMSJLacI');
TAlength{8} = fullfile('Oid-E,T 89 to 100-O1','OidTA96O1_100pMSJLacI');
TAlength{9} = fullfile('Oid-E,T 89 to 100-O1','OidTA97O1_100pMSJLacI');
TAlength{10} = fullfile('Oid-E,T 89 to 100-O1','OidTA98O1_100pMSJLacI');
TAlength{11} = fullfile('Oid-E,T 89 to 100-O1','OidTA99O1_100pMSJLacI');
TAlength{12} = fullfile('Oid-E,T 89 to 100-O1','OidTA100O1_100pMSJLacI');
TAlength{13} = fullfile('Oid-E,T 89 to 100-O1','OidTA101O1_100pMSJLacI');
TAlength{14} = fullfile('Oid-E,T 89 to 100-O1','OidTA102O1_100pMSJLacI');
TAlength{15} = fullfile('Oid-E,T 89 to 100-O1','OidTA103O1_100pMSJLacI');
TAlength{16} = fullfile('Oid-E,T 89 to 100-O1','OidTA104O1_100pMSJLacI');
TAlength{17} = fullfile('Oid-E,T 89 to 100-O1','OidTA105O1_100pMSJLacI');
TAlength{18} = fullfile('Oid-E,T 89 to 100-O1','OidTA106O1_100pMSJLacI');
TAlength{19} = fullfile('Oid-E,T 89 to 100-O1','OidTA107O1_100pMSJLacI');
TAlength{20} = fullfile('Oid-E,T 89 to 100-O1','OidTA108O1_100pMSJLacI');
TAlength{21} = fullfile('Oid-E,T 89 to 100-O1','OidTA109O1_100pMSJLacI');
TAlength{22} = fullfile('Oid-E,T 89 to 100-O1','OidTA116O1_100pMSJLacI');

%HGseqs

lengthsEHG = [56 57 58 59 60 61 65 67 69 70 71 72 78 79 80 81 82 83 84 85 86 87 88];
lengthsTHG = [56 57 58 59 60 61 65 67 69 70 71 72 78 79 80 81 82 83 84 85 86 87 88];

ElengthHG{1} = 'OidE856O2wprom_100pMSJLacI';
ElengthHG{2} = 'OidE857O2wprom_100pMSJLacI';
ElengthHG{3} = 'OidE858O2wprom_100pMSJLacI';
ElengthHG{4} = 'OidE859O2wprom_100pMSJLacI';
ElengthHG{5} = 'OidE860O2wprom_100pMSJLacI';
ElengthHG{6} = 'OidE861O2wprom_100pMSJLacI';
ElengthHG{7} = 'OidE865O2wprom_100pMSJLacI';
ElengthHG{8} = 'OidE867O2wprom_100pMSJLacI';
ElengthHG{9} = 'OidE869O2wprom_100pMSJLacI';
ElengthHG{10} = 'OidE870O2wprom_100pMSJLacI';
ElengthHG{11} = 'OidE871O2wprom_100pMSJLacI';
ElengthHG{12} = 'OidE872O2wprom_100pMSJLacI';
%ElengthHG{} = 'OidE878O2wprom_smallbds_100pMSJLacI_thresh_pt07fB_pt0461fG';
ElengthHG{13} = 'OidE878O2wprom_100pMSJLacI';
%ElengthHG{} = 'OidE879O2wprom_smallbds_100pMSJLacI_thresh_pt07fB_pt0461fG';
ElengthHG{14} = 'OidE879O2wprom_100pMSJLacI';
ElengthHG{15} = 'OidE880O2wprom_100pMSJLacI';
%ElengthHG{} = 'OidE881O2wprom_100pMSJLacI';
ElengthHG{16} = 'OidE881O2wprom_100pMSJLacI_RIGHTDNA';
%ElengthHG{} = 'OidE882O2wprom_100pMSJLacI';
ElengthHG{17} = 'OidE882O2wprom_100pMSJLacI_RIGHTDNA';
ElengthHG{18} = 'OidE883O2wprom_100pMSJLacI';
ElengthHG{19} = 'OidE884O2wprom_100pMSJLacI';
ElengthHG{20} = 'OidE885O2wprom_100pMSJLacI';
ElengthHG{21} = 'OidE886O2wprom_100pMSJLacI';
ElengthHG{22} = 'OidE887O2wprom_100pMSJLacI';
ElengthHG{23} = 'OidE888O2wprom_100pMSJLacI';

TlengthHG{1} = 'OidTA56O2wprom_100pMSJLacI';
TlengthHG{2} = 'OidTA57O2wprom_100pMSJLacI';
TlengthHG{3} = 'OidTA58O2wprom_100pMSJLacI';
TlengthHG{4} = 'OidTA59O2wprom_100pMSJLacI';
TlengthHG{5} = 'OidTA60O2wprom_100pMSJLacI';
TlengthHG{6} = 'OidTA61O2wprom_100pMSJLacI';
TlengthHG{7} = 'OidTA65O2wprom_100pMSJLacI';
TlengthHG{8} = 'OidTA67O2wprom_100pMSJLacI';
TlengthHG{9} = 'OidTA69O2wprom_100pMSJLacI';
TlengthHG{10} = 'OidTA70O2wprom_100pMSJLacI';
TlengthHG{11} = 'OidTA71O2wprom_100pMSJLacI';
TlengthHG{12} = 'OidTA72O2wprom_100pMSJLacI';
%TlengthHG{} = 'OidTA78O2wprom_smallbds_100pMSJLacI_thresh_pt07fB_pt0461fG';
TlengthHG{13} = 'OidTA78O2wprom_100pMSJLacI';
%TlengthHG{} = 'OidTA79O2wprom_smallbds_100pMSJLacI';
TlengthHG{14} = 'OidTA79O2wprom_100pMSJLacI';
TlengthHG{15} = 'OidTA80O2wprom_100pMSJLacI';
TlengthHG{16} = 'OidTA81O2wprom_100pMSJLacI';
TlengthHG{17} = 'OidTA82O2wprom_100pMSJLacI';
TlengthHG{18} = 'OidTA83O2wprom_100pMSJLacI';
TlengthHG{19} = 'OidTA84O2wprom_100pMSJLacI';
TlengthHG{20} = 'OidTA85O2wprom_100pMSJLacI';
TlengthHG{21} = 'OidTA86O2wprom_100pMSJLacI';
TlengthHG{22} = 'OidTA87O2wprom_100pMSJLacI';
TlengthHG{23} = 'OidTA88O2wprom_100pMSJLacI';

%%%% CONTROLS %%%%
%Linked channels
concLinkedChannels = 10^-12.*[5 100 500];
linkedchannels{1} = 'E894_5pMSJLacI_daisychain';
linkedchannels{2} = 'E894_100pMSJLacI_daisychain';
linkedchannels{3} = 'E894_500pMSJLacI_daisychain';

%Dry ice
concDI = 10^-12.*[1 5 100 10^3];
DI{1} = 'E894_1pMSJLacI_dryice_moredata';
DI{2} = 'E894_5pMSJLacI_dryice';
DI{3} = 'E894_100pMSJLacI_dryice';
DI{4} = 'E894_1nMSJLacI_dryice';

%No heparin
conchep = 10^-12.*[10 10^3];
hep{1} = 'E894_10pMSJLacI_nohep';
hep{2} = 'E894_1nMSJLacI_nohep';

%Small beads
concsmallbds = 10^-12.*[0.5 1 10 100 10^3];
smallbds{1} = 'E894_smallbds_500fMSJLacI';
smallbds{2} = 'E894_smallbds_1pMSJLacI';
smallbds{3} = 'E894_smallbds_10pMSJLacI';
smallbds{4} = 'E894_smallbds_100pMSJLacI';
smallbds{5} = 'E894_smallbds_1nMSJLacI';

%Second repressor batch
concSJLac2E8 = 10^-12.*[1 5 100 10^3];
E8SJLac2{1} = 'E894_1pMSJLacI2';
E8SJLac2{2} = 'E894_5pMSJLacI2';
E8SJLac2{3} = 'E894_100pMSJLacI2_moredata';
E8SJLac2{4} = 'E894_1nMSJLacI2';

concSJLac2TA = 10^-12.*[1];
TASJLac2{1} = 'TA94_1pMSJLacI2';

%%%% OLD BATCHES %%%%
%"New" repressor (KMLacI2)
concNew=10^-12.*[0.5 1 5 10 100 10^3 10^4];
New{1} = 'E8_94, 500 fM NEW lac, no BSA, 23 deg, BF';
New{2} = 'E8_94, 1 pM NEW lac, no BSA, 23 deg, BF';
New{3} = 'E8_94, 5 pM NEW lac, no BSA, 23 deg, BF';
New{4} = 'E8_94, 10 pM NEW lac, no BSA, 23 deg, BF';
New{5} = 'E8_94, 100 pM NEW lac, no BSA, 23 deg, BF';
New{6} = 'E8_94, 1 nM NEW lac, no BSA, 23 deg, BF';
New{7} = 'E8_94, 10 nM NEW lac, no BSA, 23 deg, BF';


%"Old" repressor (KMLacI), Lin's repressor
concOld=10^(-12).*[0.1 0.5 1 5 10 100 1000 10^4 10^5];
Old{1} = 'E8_94, 100 fM lac, no BSA, 23 deg, BF';
Old{2} = 'E8_94, 500 fM lac, no BSA, 23 deg, BF';
Old{3} = 'E8_94, 1 pM lac, no BSA, 23 deg, BF';
Old{4} = 'E8_94, 5 pM lac, no BSA, 23 deg, BF';
Old{5} = 'E8_94, 10 pM lac, no BSA, 23 deg, BF';
Old{6} = 'E8_94, 100 pM lac, no BSA, 23 deg, BF';
Old{7} = 'E8_94, 1 nM lac, no BSA, 23 deg, BF';
Old{8} = 'E8_94, 10 nM lac, no BSA, 23 deg, BF';
Old{9} = 'E8_94, 100 nM lac, no BSA, 23 deg, BF';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find the correct mother directory: organized by repressor batch
if ~ismac
    pathroot = 'Z:\';
end

if nargin<3 || (nargin>=3 && strcmpi(varargin{1},'SJLacI'))
    repbatch = 'SJLacI';
    if(ismac)
        pathroot='/Volumes/dumbo/stephj/TPM data analysis/SJLacI';
        %pathroot='/Volumes/RANDIR2';
    else
        pathroot=strcat(pathroot,'stephj\TPM data analysis\SJLacI');
    end
else
    repbatch = varargin{1};
    if strcmpi(repbatch,'SJLacI2')
        if(ismac)
            pathroot='/Volumes/dumbo-4/stephj/TPM data analysis/SJLacI2';
        else
            pathroot=strcat(pathroot,'stephj\TPM data analysis\SJLacI2');
        end
    elseif strcmpi(repbatch,'Old') || strcmpi(repbatch,'New')
        if(ismac)
            pathroot='/Volumes/dumbo-4/stephj/TPM data analysis/Shortloops titrations BEST';
        else
            pathroot=strcat(pathroot,'stephj\TPM data analysis\Shortloops titrations BEST');
        end
    else
        disp('Invalid repressor batch entered.')
        return
    end
end

%Create the path for the requested DNA type
if strcmpi(DNAset,'WTLac_noO1_O3toO1') || strcmpi(DNAset,'WTLac_noO3')
    path = fullfile(pathroot,'WTLac conc curves',DNAset);
elseif strcmpi(DNAset,'noOidE872O1') || strcmpi(DNAset,'noOidE872noO1')
    path = fullfile(pathroot,'Operator KO controls',DNAset);
elseif strcmpi(DNAset,'Oid-E,T 89 to 100-O1') || strcmpi(repbatch,'Old') || strcmpi(repbatch,'New')
    path = pathroot; %The length series contains data in different folders, so the names 
        %in the lists of good above contain subfolder information.  "Old"
        %and "New" don't have subfolders.
elseif strcmpi(DNAset,'E894 conc curve') || strcmpi(DNAset,'TA94 conc curve') || ...
        strcmpi(DNAset,'PUC306 conc curve') || strcmpi(DNAset,'OidE94Oid conc curve') ||...
        strcmpi(DNAset,'Indicia Beads') || ...
        strcmpi(DNAset,'dry ice') || strcmpi(DNAset,'O1E894O1 conc curve') || ...
        strcmpi(DNAset,'O2E894O1 conc curve')
    path = fullfile(pathroot,DNAset);
elseif strcmpi(DNAset,'linked channels') || strcmpi(DNAset,'heparin')
    path = fullfile(pathroot,'E894, cause of f shift');
elseif strcmpi(DNAset,'PUC check')
    path = fullfile(pathroot,'PUC306 conc curve');
elseif strcmpi(DNAset,'HGseqsWProm')
    path = fullfile(pathroot,'HGseqs (Oid-E,T 80 to 90-O2 w prom)');
elseif strcmpi(DNAset,'E898conccurve') || strcmpi(DNAset,'E8107conccurve') || ...
        strcmpi(DNAset,'E8108conccurve')
    path = fullfile(pathroot,'Oid-E,T 89 to 100-O1');
else
    disp('Invalid DNAset entered.')
    return
end

%Load the requested data
if strcmpi(DNAset,'noOidE872noO1')
    conc = concnoOidE872noO1;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,noOidE872noO1conc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,noOidE872noO1conc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,noOidE872noO1conc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'Loop Lengths')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,noOidE872noO1conc{i},'Loop Lengths','lengths.mat'));
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'noOidE872O1')
    conc = concnoOidE872O1;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,noOidE872O1conc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,noOidE872O1conc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,noOidE872O1conc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'Loop Lengths')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,noOidE872O1conc{i},'Loop Lengths','lengths.mat'));
        end
    else
        disp('Invalid data kind entered.')
        return
    end

elseif strcmpi(DNAset,'E894 conc curve') && strcmpi(repbatch,'SJLacI') && ...
        ~strcmpi(repbatch,'SJLacI2') && ~strcmpi(repbatch,'Old') && ~strcmpi(repbatch,'New')
    conc = concE8;
    
    if strcmpi(kind,'names')
        data = E8conc;
    elseif strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8conc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8conc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8conc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'E8HistoAnal','E8HistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        elseif strcmpi(varargin{2},'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'alldistribssub02')
            data = datatemp.alldistribssub02;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'Meanssub02')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'TA94 conc curve') && strcmpi(repbatch,'SJLacI') && ...
        ~strcmpi(repbatch,'SJLacI2') && ~strcmpi(repbatch,'Old') && ~strcmpi(repbatch,'New')
    conc = concTA;

    if strcmpi(kind,'names')
        data = TAconc;
    elseif strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,TAconc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,TAconc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,TAconc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'TAHistoAnal','TAHistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        elseif strcmpi(varargin{2},'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'alldistribssub02')
            data = datatemp.alldistribssub02;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'Meanssub02')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
        
elseif strcmpi(DNAset,'PUC306 conc curve')
    conc = concPUC;

    if strcmpi(kind,'GaussFit')
        data = cell(length(concPUC),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,PUCconc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(concPUC),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,PUCconc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(concPUC),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,PUCconc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        if ~strcmpi(varargin{2},'alldistribsB') && ~strcmpi(varargin{2},'alldistribsM')
            datatemp = load(fullfile(path,'PUCHistoAnal','PUCHistoAnal.mat'));
            if strcmpi(varargin{2},'Means')
                datatemp2 = datatemp.Means;
                datatemp3 = datatemp.SEs;
                data = {datatemp2 datatemp3};
            elseif strcmpi(varargin{2},'Meanssub0')
                datatemp2 = datatemp.Meanssub0;
                datatemp3 = datatemp.SEssub0;
                data = {datatemp2 datatemp3};
            elseif strcmpi(varargin{2},'Gaussmeans')
                datatemp2 = datatemp.Gaussmeans;
                datatemp3 = datatemp.GaussSEs;
                data = {datatemp2 datatemp3};
            elseif strcmpi(varargin{2},'Weightmeans')
                data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
            elseif strcmpi(varargin{2},'Weightmeanssub0')
                data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
            elseif strcmpi(varargin{2},'alldistribs')
                data = datatemp.alldistribs;
            elseif strcmpi(varargin{2},'alldistribssub0')
                data = datatemp.alldistribssub0;
            else
                disp('Invalid analmethod.')
                return
            end
        else
            if strcmpi(varargin{2},'alldistribsB')
                data = cell(length(concPUC),1);
                for i=1:length(conc)
                    tempGauss = load(fullfile(path,PUCconc{i},'Single bead analysis','GaussFit.mat'));
                    [notneeded, data{i}, notneeded] = pLoopsSortByTime(tempGauss);
                end
            elseif strcmpi(varargin{2},'alldistribsM')
                data = cell(length(concPUC),1);
                for i=1:length(conc)
                    tempGauss = load(fullfile(path,PUCconc{i},'Single bead analysis','GaussFit.mat'));
                    [notneeded, notneeded, data{i}] = pLoopsSortByTime(tempGauss);
                end
            end
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Gaussmeans')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'O1E894O1 conc curve')
    conc = concE8O1O1;

    if strcmpi(kind,'names')
        data = E8O1O1conc;
    elseif strcmpi(kind,'GaussFit')
        data = cell(length(concE8O1O1),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8O1O1conc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(concE8O1O1),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8O1O1conc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(concE8O1O1),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8O1O1conc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'O1E8O1HistoAnal','O1E8O1HistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        elseif strcmpi(varargin{2},'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'alldistribssub02')
            data = datatemp.alldistribssub02;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'Meanssub02')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'O2E894O1 conc curve')
    conc = concE8O2O1;

    if strcmpi(kind,'names')
        data = E8O2O1conc;
    elseif strcmpi(kind,'GaussFit')
        data = cell(length(concE8O2O1),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8O2O1conc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(concE8O2O1),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8O2O1conc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(concE8O2O1),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8O2O1conc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'O2E8O1HistoAnal','O2E8O1HistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        elseif strcmpi(varargin{2},'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'alldistribssub02')
            data = datatemp.alldistribssub02;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'Meanssub02')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
elseif strcmpi(DNAset,'Oid-E,T 89 to 100-O1')
    conc{1} = lengthsE8;
    conc{2} = lengthsTA;
    
    if strcmpi(kind,'names')
        data{1} = E8length;
        data{2} = TAlength;
    elseif strcmpi(kind,'GaussFit')
        data{1} = cell(length(conc{1}),1);
        data{2} = cell(length(conc{2}),1);
        for i=1:length(conc{1})
            data{1}{i} = load(fullfile(path,E8length{i},'Single bead analysis','GaussFit.mat'));
        end
        for i=1:length(conc{2})
            data{2}{i} = load(fullfile(path,TAlength{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'ThreshFit')
        data{1} = cell(length(conc{1}),1);
        data{2} = cell(length(conc{2}),1);
        %As of 7/2011, loading all of the ThreshFits for this data set
        %exceeds my memory limits, probably because the traces are
        %included.  So returning structures without the trace field.
        disp('Cannot load all ThreshFit files with current memory limits.')
        disp('Use paths in data.') 
        for i=1:length(conc{1})
            %data{1}{i} = load(fullfile(path,E8length{i},'Thresholding','ThreshFit.mat'));
            %temp = load(fullfile(path,E8length{i},'Thresholding','ThreshFit.mat'));
%             temp = temp.ThreshFit;
%             for dd = 1:length(temp)
%                 ThreshFit(dd).name = temp(dd).name;
%                 ThreshFit(dd).Approved = temp(dd).Approved;
%                 ThreshFit(dd).Identity = temp(dd).Identity;
%                 ThreshFit(dd).topline = temp(dd).topline;
%                 ThreshFit(dd).midline = temp(dd).midline;
%                 ThreshFit(dd).bottomline = temp(dd).bottomline;
%                 ThreshFit(dd).p = temp(dd).p;
%             end
%             data{1}{i} = ThreshFit;
%             clear temp ThreshFit
            %Better way!
%             ThreshFit = load(fullfile(path,E8length{i},'Thresholding','ThreshFit.mat'));
%             ThreshFit = ThreshFit.ThreshFit;
%             ThreshFit = rmfield(ThreshFit,'trace');
%             data{1}{i} = ThreshFit;
%             clear ThreshFit;
            %Drat but I do need "trace" for almost everything, including
            %MakepLoopsHistos.
            data{1}{i} = fullfile(path,E8length{i},'Thresholding','ThreshFit.mat');
        end
        for i=1:length(conc{2})
            %data{2}{i} = load(fullfile(path,TAlength{i},'Thresholding','ThreshFit.mat'));
%             ThreshFit = load(fullfile(path,TAlength{i},'Thresholding','ThreshFit.mat'));
%             ThreshFit = ThreshFit.ThreshFit;
%             ThreshFit = rmfield(ThreshFit,'trace');
%             data{2}{i} = ThreshFit;
%             clear ThreshFit;
            data{2}{i} = fullfile(path,TAlength{i},'Thresholding','ThreshFit.mat');
        end    
        
    elseif strcmpi(kind,'stats')
        data{1} = cell(length(conc{1}),1);
        data{2} = cell(length(conc{2}),1);
        for i=1:length(conc{1})
            data{1}{i} = load(fullfile(path,E8length{i},'Single bead analysis','stats.mat'));
        end
        for i=1:length(conc{2})
            data{2}{i} = load(fullfile(path,TAlength{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'Threshstats')
        data{1} = cell(length(conc{1}),1);
        data{2} = cell(length(conc{2}),1);
        for i=1:length(conc{1})
            data{1}{i} = load(fullfile(path,E8length{i},'Thresholding','stats.mat'));
        end
        for i=1:length(conc{2})
            data{2}{i} = load(fullfile(path,TAlength{i},'Thresholding','stats.mat'));
        end
        disp('Thresh stats INCORRECT as of 3-2011!')
    elseif strcmpi(kind,'lacdataconcat')
        data{1} = cell(length(conc{1}),1);
        data{2} = cell(length(conc{2}),1);
        for i=1:length(conc{1})
            data{1}{i} = load(fullfile(path,E8length{i},'lacdataconcat.mat'));
        end
        for i=1:length(conc{2})
            data{2}{i} = load(fullfile(path,TAlength{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal') || strcmpi(kind,'ThreshAnal')
        if strcmpi(kind,'HistoAnal')
            datatempE8 = load(fullfile(path,'Oid-E,T 89 to 100-O1','E8lengthsHistoAnal','E8lengthsHistoAnal.mat'));
            datatempTA = load(fullfile(path,'Oid-E,T 89 to 100-O1','TAlengthsHistoAnal','TAlengthsHistoAnal.mat'));
        else
            datatempE8 = load(fullfile(path,'Oid-E,T 89 to 100-O1','E8lengthsThreshAnal','E8lengthsThreshAnal.mat'));
            datatempTA = load(fullfile(path,'Oid-E,T 89 to 100-O1','TAlengthsThreshAnal','TAlengthsThreshAnal.mat'));
        end
        if strcmpi(varargin{2},'Means')
            datatemp2E8 = datatempE8.Means;
            datatemp3E8 = datatempE8.SEs;
            datatemp2TA = datatempTA.Means;
            datatemp3TA = datatempTA.SEs;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'MeansB') && isfield(datatempE8,'MeansB')
            datatemp2E8 = datatempE8.MeansB;
            datatemp3E8 = datatempE8.SEsB;
            datatemp2TA = datatempTA.MeansB;
            datatemp3TA = datatempTA.SEsB;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'MeansM') && isfield(datatempE8,'MeansM')
            datatemp2E8 = datatempE8.MeansM;
            datatemp3E8 = datatempE8.SEsM;
            datatemp2TA = datatempTA.MeansM;
            datatemp3TA = datatempTA.SEsM;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2E8 = datatempE8.Meanssub0;
            datatemp3E8 = datatempE8.SEssub0;
            datatemp2TA = datatempTA.Meanssub0;
            datatemp3TA = datatempTA.SEssub0;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Meanssub02') && isfield(datatempE8,'Meanssub02')
            datatemp2E8 = datatempE8.Meanssub02;
            datatemp3E8 = datatempE8.SEssub02;
            datatemp2TA = datatempTA.Meanssub02;
            datatemp3TA = datatempTA.SEssub02;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Meanssub02B') && isfield(datatempE8,'Meanssub02B')
            datatemp2E8 = datatempE8.Meanssub02B;
            datatemp3E8 = datatempE8.SEssub02B;
            datatemp2TA = datatempTA.Meanssub02B;
            datatemp3TA = datatempTA.SEssub02B;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Meanssub02M') && isfield(datatempE8,'Meanssub02M')
            datatemp2E8 = datatempE8.Meanssub02M;
            datatemp3E8 = datatempE8.SEssub02M;
            datatemp2TA = datatempTA.Meanssub02M;
            datatemp3TA = datatempTA.SEssub02M;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2E8 = datatempE8.Gaussmeans;
            datatemp3E8 = datatempE8.GaussSEs;
            datatemp2TA = datatempTA.Gaussmeans;
            datatemp3TA = datatempTA.GaussSEs;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Weightmeans')
            data{1} = [datatempE8.pLoopWeightMeans,datatempE8.SEWeightMean];
            data{2} = [datatempTA.pLoopWeightMeans,datatempTA.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data{1} = [datatempE8.pLoopWeightMeanssub0,datatempE8.SEWeightMeansub0];
            data{2} = [datatempTA.pLoopWeightMeanssub0,datatempTA.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data{1} = datatempE8.alldistribs;
            data{2} = datatempTA.alldistribs;
        elseif strcmpi(varargin{2},'alldistribsB') && isfield(datatempE8,'alldistribsB')
            data{1} = datatempE8.alldistribsB;
            data{2} = datatempTA.alldistribsB;
        elseif strcmpi(varargin{2},'alldistribsM') && isfield(datatempE8,'alldistribsM')
            data{1} = datatempE8.alldistribsM;
            data{2} = datatempTA.alldistribsM;
        elseif strcmpi(varargin{2},'alldistribssub02') && isfield(datatempE8,'alldistribssub02')
            data{1} = datatempE8.alldistribssub02;
            data{2} = datatempTA.alldistribssub02;
        elseif strcmpi(varargin{2},'alldistribssub02B') && isfield(datatempE8,'alldistribssub02B')
            data{1} = datatempE8.alldistribssub02B;
            data{2} = datatempTA.alldistribssub02B;
        elseif strcmpi(varargin{2},'alldistribssub02M') && isfield(datatempE8,'alldistribssub02M')
            data{1} = datatempE8.alldistribssub02M;
            data{2} = datatempTA.alldistribssub02M;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc{1})
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub02') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'MeansB') || ...
                    strcmpi(varargin{2},'Meanssub02B') || strcmpi(varargin{2},'MeansM') || ...
                    strcmpi(varargin{2},'Meanssub02M')
                    data{1}(k,1) = datatemp4{1}{1}{k}(numsecs/1000);
                    data{1}(k,2) = datatemp4{1}{2}{k}(numsecs/1000);
                else
                    data{1}{k}=datatemp4{1}{k}{numsecs/1000}(:,2);
                end
            end
            for k = 1:length(conc{2})
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub02') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'MeansB') || ...
                    strcmpi(varargin{2},'Meanssub02B') || strcmpi(varargin{2},'MeansM') || ...
                    strcmpi(varargin{2},'Meanssub02M')
                    data{2}(k,1) = datatemp4{2}{1}{k}(numsecs/1000);
                    data{2}(k,2) = datatemp4{2}{2}{k}(numsecs/1000);
                else
                    data{2}{k}=datatemp4{2}{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'E898conccurve')
    conc = concE898;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E898conc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'ThreshFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E898conc{i},'Thresholding','ThreshFit.mat'));
        end   
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E898conc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'Threshstats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E898conc{i},'Thresholding','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E898conc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal') || strcmpi(kind,'ThreshAnal')
        if strcmpi(kind,'HistoAnal')
            disp('Not done yet!')
            return
            %datatemp = load(fullfile(path,'E898conccurveHistoAnal','E898conccurveThreshAnal.mat'));
        else
            datatemp = load(fullfile(path,'E898conccurveHistoAnal','E898conccurveThreshAnal.mat'));
        end
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'MeansB') && isfield(datatemp,'MeansB')
            datatemp2 = datatemp.MeansB;
            datatemp3 = datatemp.SEsB;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'MeansM') && isfield(datatemp,'MeansM')
            datatemp2 = datatemp.MeansM;
            datatemp3 = datatemp.SEsM;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0B') && isfield(datatemp,'Meanssub0B')
            datatemp2 = datatemp.Meanssub0B;
            datatemp3 = datatemp.SEssub0B;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0M') && isfield(datatemp,'Meanssub0M')
            datatemp2 = datatemp.Meanssub0M;
            datatemp3 = datatemp.SEssub0M;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02') && isfield(datatemp,'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02B') && isfield(datatemp,'Meanssub02B')
            datatemp2 = datatemp.Meanssub02B;
            datatemp3 = datatemp.SEssub02B;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02M') && isfield(datatemp,'Meanssub02M')
            datatemp2 = datatemp.Meanssub02M;
            datatemp3 = datatemp.SEssub02M;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribsB') && isfield(datatemp,'alldistribsB')
            data = datatemp.alldistribsB;
        elseif strcmpi(varargin{2},'alldistribsM') && isfield(datatemp,'alldistribsM')
            data = datatemp.alldistribsM;
        elseif strcmpi(varargin{2},'alldistribssub02') && isfield(datatemp,'alldistribssub02')
            data = datatemp.alldistribssub02;
        elseif strcmpi(varargin{2},'alldistribssub02B') && isfield(datatemp,'alldistribssub02B')
            data = datatemp.alldistribssub02B;
        elseif strcmpi(varargin{2},'alldistribssub02M') && isfield(datatemp,'alldistribssub02M')
            data = datatemp.alldistribssub02M;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'MeansB') || ...
                    strcmpi(varargin{2},'Meanssub0B') || strcmpi(varargin{2},'MeansM') || ...
                    strcmpi(varargin{2},'Meanssub0M')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'E8107conccurve')
    conc = concE8107;
    
    if strcmpi(kind,'names')
        data = E8107conc;
    elseif strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8107conc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'ThreshFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8107conc{i},'Thresholding','ThreshFit.mat'));
        end   
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8107conc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'Threshstats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8107conc{i},'Thresholding','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8107conc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal') || strcmpi(kind,'ThreshAnal')
        if strcmpi(kind,'HistoAnal')
            disp('Not done yet!')
            return
            %datatemp = load(fullfile(path,'E898conccurveHistoAnal','E898conccurveThreshAnal.mat'));
        else
            datatemp = load(fullfile(path,'E8107conccurveThreshAnal','E8107conccurveThreshAnal.mat'));
        end
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'MeansB') && isfield(datatemp,'MeansB')
            datatemp2 = datatemp.MeansB;
            datatemp3 = datatemp.SEsB;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'MeansM') && isfield(datatemp,'MeansM')
            datatemp2 = datatemp.MeansM;
            datatemp3 = datatemp.SEsM;
            data = {datatemp2 datatemp3};
%         elseif strcmpi(varargin{2},'Meanssub0')
%             datatemp2 = datatemp.Meanssub0;
%             datatemp3 = datatemp.SEssub0;
%             data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0B') && isfield(datatemp,'Meanssub0B')
            datatemp2 = datatemp.Meanssub0B;
            datatemp3 = datatemp.SEssub0B;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0M') && isfield(datatemp,'Meanssub0M')
            datatemp2 = datatemp.Meanssub0M;
            datatemp3 = datatemp.SEssub0M;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02') && isfield(datatemp,'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02B') && isfield(datatemp,'Meanssub02B')
            datatemp2 = datatemp.Meanssub02B;
            datatemp3 = datatemp.SEssub02B;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02M') && isfield(datatemp,'Meanssub02M')
            datatemp2 = datatemp.Meanssub02M;
            datatemp3 = datatemp.SEssub02M;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
%         elseif strcmpi(varargin{2},'Weightmeans')
%             data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
%         elseif strcmpi(varargin{2},'Weightmeanssub0')
%             data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
%         elseif strcmpi(varargin{2},'alldistribs')
%             data = datatemp.alldistribs;
%         elseif strcmpi(varargin{2},'alldistribsB') && isfield(datatemp,'alldistribsB')
%             data = datatemp.alldistribsB;
%         elseif strcmpi(varargin{2},'alldistribsM') && isfield(datatemp,'alldistribsM')
%             data = datatemp.alldistribsM;
        elseif strcmpi(varargin{2},'alldistribssub02') && isfield(datatemp,'alldistribssub02')
            data = datatemp.alldistribssub02;
        elseif strcmpi(varargin{2},'alldistribssub02B') && isfield(datatemp,'alldistribssub02B')
            data = datatemp.alldistribssub02B;
        elseif strcmpi(varargin{2},'alldistribssub02M') && isfield(datatemp,'alldistribssub02M')
            data = datatemp.alldistribssub02M;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub02') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'MeansB') || ...
                    strcmpi(varargin{2},'Meanssub02B') || strcmpi(varargin{2},'MeansM') || ...
                    strcmpi(varargin{2},'Meanssub02M')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'E8108conccurve')
    conc = concE8108;
    
    if strcmpi(kind,'names')
        data = E8108conc;
    elseif strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8108conc{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'ThreshFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8108conc{i},'Thresholding','ThreshFit.mat'));
        end   
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8108conc{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'Threshstats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8108conc{i},'Thresholding','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8108conc{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal') || strcmpi(kind,'ThreshAnal')
        if strcmpi(kind,'HistoAnal')
            disp('Not done yet!')
            return
            %datatemp = load(fullfile(path,'E898conccurveHistoAnal','E898conccurveThreshAnal.mat'));
        else
            datatemp = load(fullfile(path,'E8108conccurveThreshAnal','E8108conccurveThreshAnal.mat'));
        end
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'MeansB') && isfield(datatemp,'MeansB')
            datatemp2 = datatemp.MeansB;
            datatemp3 = datatemp.SEsB;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'MeansM') && isfield(datatemp,'MeansM')
            datatemp2 = datatemp.MeansM;
            datatemp3 = datatemp.SEsM;
            data = {datatemp2 datatemp3};
%         elseif strcmpi(varargin{2},'Meanssub0')
%             datatemp2 = datatemp.Meanssub0;
%             datatemp3 = datatemp.SEssub0;
%             data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0B') && isfield(datatemp,'Meanssub0B')
            datatemp2 = datatemp.Meanssub0B;
            datatemp3 = datatemp.SEssub0B;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0M') && isfield(datatemp,'Meanssub0M')
            datatemp2 = datatemp.Meanssub0M;
            datatemp3 = datatemp.SEssub0M;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02') && isfield(datatemp,'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02B') && isfield(datatemp,'Meanssub02B')
            datatemp2 = datatemp.Meanssub02B;
            datatemp3 = datatemp.SEssub02B;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02M') && isfield(datatemp,'Meanssub02M')
            datatemp2 = datatemp.Meanssub02M;
            datatemp3 = datatemp.SEssub02M;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
%         elseif strcmpi(varargin{2},'Weightmeans')
%             data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
%         elseif strcmpi(varargin{2},'Weightmeanssub0')
%             data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
%         elseif strcmpi(varargin{2},'alldistribs')
%             data = datatemp.alldistribs;
%         elseif strcmpi(varargin{2},'alldistribsB') && isfield(datatemp,'alldistribsB')
%             data = datatemp.alldistribsB;
%         elseif strcmpi(varargin{2},'alldistribsM') && isfield(datatemp,'alldistribsM')
%             data = datatemp.alldistribsM;
        elseif strcmpi(varargin{2},'alldistribssub02') && isfield(datatemp,'alldistribssub02')
            data = datatemp.alldistribssub02;
        elseif strcmpi(varargin{2},'alldistribssub02B') && isfield(datatemp,'alldistribssub02B')
            data = datatemp.alldistribssub02B;
        elseif strcmpi(varargin{2},'alldistribssub02M') && isfield(datatemp,'alldistribssub02M')
            data = datatemp.alldistribssub02M;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub02') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'MeansB') || ...
                    strcmpi(varargin{2},'Meanssub02B') || strcmpi(varargin{2},'MeansM') || ...
                    strcmpi(varargin{2},'Meanssub02M')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end    
    
elseif strcmpi(DNAset,'HGseqsWProm')
    conc{1} = lengthsEHG;
    conc{2} = lengthsTHG;
    
    if strcmpi(kind,'names')
        data{1} = ElengthHG;
        data{2} = TlengthHG;
    elseif strcmpi(kind,'GaussFit')
        data{1} = cell(length(conc{1}),1);
        data{2} = cell(length(conc{2}),1);
        for i=1:length(conc{1})
            data{1}{i} = load(fullfile(path,ElengthHG{i},'Single bead analysis','GaussFit.mat'));
        end
        for i=1:length(conc{2})
            if ~strcmpi(TlengthHG{i},'OidTA80O2wprom_100pMSJLacI') && ~strcmpi(TlengthHG{i},'OidTA79O2wprom_100pMSJLacI')
                data{2}{i} = load(fullfile(path,TlengthHG{i},'Single bead analysis','GaussFit.mat'));
            else
                data{2}{i} = load(fullfile(path,TlengthHG{i},'Thresholding','ThreshFit.mat'));
            end
        end 
    elseif strcmpi(kind,'ThreshFit')
        data{1} = cell(length(conc{1}),1);
        data{2} = cell(length(conc{2}),1);
        for i=1:length(conc{1})
            data{1}{i} = load(fullfile(path,ElengthHG{i},'Thresholding','ThreshFit.mat'));
        end
        for i=1:length(conc{2})
            data{2}{i} = load(fullfile(path,TlengthHG{i},'Thresholding','ThreshFit.mat'));
        end 
    elseif strcmpi(kind,'stats')
        data{1} = cell(length(conc{1}),1);
        data{2} = cell(length(conc{2}),1);
        for i=1:length(conc{1})
            data{1}{i} = load(fullfile(path,ElengthHG{i},'Single bead analysis','stats.mat'));
        end
        for i=1:length(conc{2})
            if ~strcmpi(TlengthHG{i},'OidTA80O2wprom_100pMSJLacI') && ~strcmpi(TlengthHG{i},'OidTA79O2wprom_100pMSJLacI')
                data{2}{i} = load(fullfile(path,TlengthHG{i},'Single bead analysis','stats.mat'));
            else
                data{2}{i} = load(fullfile(path,TlengthHG{i},'Thresholding','stats.mat'));
            end
        end
    elseif strcmpi(kind,'lacdataconcat')
        data{1} = cell(length(conc{1}),1);
        data{2} = cell(length(conc{2}),1);
        for i=1:length(conc{1})
            data{1}{i} = load(fullfile(path,ElengthHG{i},'lacdataconcat.mat'));
        end
        for i=1:length(conc{2})
            data{2}{i} = load(fullfile(path,TlengthHG{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal') || strcmpi(kind,'ThreshAnal')
        if strcmpi(kind,'HistoAnal')
            disp('Should use thresh anal.')
            datatempE8 = load(fullfile(path,'HGseqsE8HistoAnal','HGseqsE8HistoAnal.mat'));
            datatempTA = load(fullfile(path,'HGseqsTAHistoAnal','HGseqsTAHistoAnal.mat'));
        else
            datatempE8 = load(fullfile(path,'HGseqsE8ThreshAnal','HGseqsE8ThreshAnal.mat'));
            datatempTA = load(fullfile(path,'HGseqsTAThreshAnal','HGseqsTAThreshAnal.mat'));
        end

        if strcmpi(varargin{2},'Means')
            datatemp2E8 = datatempE8.Means;
            datatemp3E8 = datatempE8.SEs;
            datatemp2TA = datatempTA.Means;
            datatemp3TA = datatempTA.SEs;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'MeansB') && isfield(datatempE8,'MeansB')
            datatemp2E8 = datatempE8.MeansB;
            datatemp3E8 = datatempE8.SEsB;
            datatemp2TA = datatempTA.MeansB;
            datatemp3TA = datatempTA.SEsB;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'MeansM') && isfield(datatempE8,'MeansM')
            datatemp2E8 = datatempE8.MeansM;
            datatemp3E8 = datatempE8.SEsM;
            datatemp2TA = datatempTA.MeansM;
            datatemp3TA = datatempTA.SEsM;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Meanssub02') && isfield(datatempE8,'Meanssub02')
            datatemp2E8 = datatempE8.Meanssub02;
            datatemp3E8 = datatempE8.SEssub02;
            datatemp2TA = datatempTA.Meanssub02;
            datatemp3TA = datatempTA.SEssub02;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Meanssub02B') && isfield(datatempE8,'Meanssub02B')
            datatemp2E8 = datatempE8.Meanssub02B;
            datatemp3E8 = datatempE8.SEssub02B;
            datatemp2TA = datatempTA.Meanssub02B;
            datatemp3TA = datatempTA.SEssub02B;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Meanssub02M') && isfield(datatempE8,'Meanssub02M')
            datatemp2E8 = datatempE8.Meanssub02M;
            datatemp3E8 = datatempE8.SEssub02M;
            datatemp2TA = datatempTA.Meanssub02M;
            datatemp3TA = datatempTA.SEssub02M;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2E8 = datatempE8.Gaussmeans;
            datatemp3E8 = datatempE8.GaussSEs;
            datatemp2TA = datatempTA.Gaussmeans;
            datatemp3TA = datatempTA.GaussSEs;
            data{1} = {datatemp2E8 datatemp3E8};
            data{2} = {datatemp2TA datatemp3TA};
        elseif strcmpi(varargin{2},'alldistribs')
            data{1} = datatempE8.alldistribs;
            data{2} = datatempTA.alldistribs;
        elseif strcmpi(varargin{2},'alldistribsB') && isfield(datatempE8,'alldistribsB')
            data{1} = datatempE8.alldistribsB;
            data{2} = datatempTA.alldistribsB;
        elseif strcmpi(varargin{2},'alldistribsM') && isfield(datatempE8,'alldistribsM')
            data{1} = datatempE8.alldistribsM;
            data{2} = datatempTA.alldistribsM;
        elseif strcmpi(varargin{2},'alldistribssub02') && isfield(datatempE8,'alldistribssub02')
            data{1} = datatempE8.alldistribssub02;
            data{2} = datatempTA.alldistribssub02;
        elseif strcmpi(varargin{2},'alldistribssub02B') && isfield(datatempE8,'alldistribssub02B')
            data{1} = datatempE8.alldistribssub02B;
            data{2} = datatempTA.alldistribssub02B;
        elseif strcmpi(varargin{2},'alldistribssub02M') && isfield(datatempE8,'alldistribssub02M')
            data{1} = datatempE8.alldistribssub02M;
            data{2} = datatempTA.alldistribssub02M;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc{1})
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub02') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'MeansB') || ...
                    strcmpi(varargin{2},'Meanssub02B') || strcmpi(varargin{2},'MeansM') || ...
                    strcmpi(varargin{2},'Meanssub02M')
                    data{1}(k,1) = datatemp4{1}{1}{k}(numsecs/1000);
                    data{1}(k,2) = datatemp4{1}{2}{k}(numsecs/1000);
                else
                    data{1}{k}=datatemp4{1}{k}{numsecs/1000}(:,2);
                end
            end
            for k = 1:length(conc{2})
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub02') || ...
                    strcmpi(varargin{2},'Gaussmeans') || strcmpi(varargin{2},'MeansB') || ...
                    strcmpi(varargin{2},'Meanssub02B') || strcmpi(varargin{2},'MeansM') || ...
                    strcmpi(varargin{2},'Meanssub02M')
                    data{2}(k,1) = datatemp4{2}{1}{k}(numsecs/1000);
                    data{2}(k,2) = datatemp4{2}{2}{k}(numsecs/1000);
                else
                    data{2}{k}=datatemp4{2}{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'linked channels')
    conc = concLinkedChannels;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,linkedchannels{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,linkedchannels{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,linkedchannels{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'LCHistoAnal','LCHistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        elseif strcmpi(varargin{2},'alldistribssub02')
            data = datatemp.alldistribssub02;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Meanssub02') || strcmpi(varargin{2},'Gaussmeans')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
elseif strcmpi(DNAset,'dry ice')
    conc = concDI;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,DI{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,DI{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,DI{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'DIHistoAnal','DIHistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Gaussmeans')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'Indicia Beads')
    conc = concsmallbds;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,smallbds{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,smallbds{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,smallbds{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'SmallBdsHistoAnal','SmallBdsHistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        elseif strcmpi(varargin{2},'alldistribssub02')
            data = datatemp.alldistribssub02;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Meanssub02') || strcmpi(varargin{2},'Gaussmeans')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end    
    
elseif strcmpi(DNAset,'heparin')
    conc = conchep;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,hep{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,hep{i},'Single bead analysis','stats.mat'));
        end
    end    
    
elseif strcmpi(DNAset,'E894 conc curve') && ~strcmpi(repbatch,'SJLacI') && ...
        strcmpi(repbatch,'SJLacI2') && ~strcmpi(repbatch,'Old') && ~strcmpi(repbatch,'New')
    conc = concSJLac2E8;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8SJLac2{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8SJLac2{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,E8SJLac2{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'E8SJLac2HistoAnal','E8SJLac2HistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};    
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Meanssub02') || strcmpi(varargin{2},'Gaussmeans')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'TA94 conc curve') && ~strcmpi(repbatch,'SJLacI') && ...
        strcmpi(repbatch,'SJLacI2') && ~strcmpi(repbatch,'Old') && ~strcmpi(repbatch,'New')
    conc = concSJLac2TA;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,TASJLac2{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,TASJLac2{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,TASJLac2{i},'lacdataconcat.mat'));
        end  
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'E894 conc curve') && ~strcmpi(repbatch,'SJLacI') &&  ...
        ~strcmpi(repbatch,'SJLacI2') && ~strcmpi(repbatch,'Old') && strcmpi(repbatch,'New')
    conc = concNew;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,New{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,New{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,New{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'E8NewLacHistoAnal','E8NewLacHistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        elseif strcmpi(varargin{2},'alldistribssub02')
            data = datatemp.alldistribssub02;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') ||...
                        strcmpi(varargin{2},'Meanssub02') || ...
                    strcmpi(varargin{2},'Gaussmeans')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end
    
elseif strcmpi(DNAset,'E894 conc curve') && ~strcmpi(repbatch,'SJLacI') && ...
        ~strcmpi(repbatch,'SJLacI2') && strcmpi(repbatch,'Old') && ~strcmpi(repbatch,'New')
    conc = concOld;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,Old{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,Old{i},'Single bead analysis','stats.mat'));
        end
    elseif strcmpi(kind,'lacdataconcat')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,Old{i},'lacdataconcat.mat'));
        end
    elseif strcmpi(kind,'HistoAnal')
        datatemp = load(fullfile(path,'E8OldLacHistoAnal','E8OldLacHistoAnal.mat'));
        if strcmpi(varargin{2},'Means')
            datatemp2 = datatemp.Means;
            datatemp3 = datatemp.SEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub0')
            datatemp2 = datatemp.Meanssub0;
            datatemp3 = datatemp.SEssub0;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Meanssub02')
            datatemp2 = datatemp.Meanssub02;
            datatemp3 = datatemp.SEssub02;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Gaussmeans')
            datatemp2 = datatemp.Gaussmeans;
            datatemp3 = datatemp.GaussSEs;
            data = {datatemp2 datatemp3};
        elseif strcmpi(varargin{2},'Weightmeans')
            data = [datatemp.pLoopWeightMeans,datatemp.SEWeightMean];
        elseif strcmpi(varargin{2},'Weightmeanssub0')
            data = [datatemp.pLoopWeightMeanssub0,datatemp.SEWeightMeansub0];
        elseif strcmpi(varargin{2},'alldistribs')
            data = datatemp.alldistribs;
        elseif strcmpi(varargin{2},'alldistribssub0')
            data = datatemp.alldistribssub0;
        elseif strcmpi(varargin{2},'alldistribssub02')
            data = datatemp.alldistribssub02;
        else
            disp('Invalid analmethod.')
            return
        end
        
        if nargin>4 && ~strcmpi(varargin{2},'Weightmeans') && ~strcmpi(varargin{2},'Weightmeanssub0') %User gave numsecs
            datatemp4=data;
            clear data
            numsecs=varargin{3};
            for k = 1:length(conc)
                if strcmpi(varargin{2},'Means') || strcmpi(varargin{2},'Meanssub0') || ...
                    strcmpi(varargin{2},'Meanssub02') || strcmpi(varargin{2},'Gaussmeans')
                    data(k,1) = datatemp4{1}{k}(numsecs/1000);
                    data(k,2) = datatemp4{2}{k}(numsecs/1000);
                else
                    data{k}=datatemp4{k}{numsecs/1000}(:,2);
                end
            end
        end
    else
        disp('Invalid data kind entered.')
        return
    end    
    
elseif strcmpi(DNAset,'PUC check')
    conc = concPUCcheck;
    
    if strcmpi(kind,'GaussFit')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,PUCcheck{i},'Single bead analysis','GaussFit.mat'));
        end
    elseif strcmpi(kind,'stats')
        data = cell(length(conc),1);
        for i=1:length(conc)
            data{i} = load(fullfile(path,PUCcheck{i},'Single bead analysis','stats.mat'));
        end
    end
    
end
 
    
    