%(script) AnalyzeTetherLengths
%
%Runs CalculateTetherLength_Gauss and _Thresh on my current data sets, and
%makes plots of the results.
%
%This version makes handling different data sets easier
%
%Steph 1/2012

%RMS vs length (see below for concentration curves.  This uses
%thresholding)

HGseqs = 1; %Make this 1 to analyze with-promoter lengths, 0 for no-promoter
Scatter = 0; %Make this 1 for a scatter plot as well as means+SE's.  
    %Note that the scatter plots can throw errors if some lengths have no
    %middle or bottom state.  I haven't fixed that yet.

if ~HGseqs
    [lengths,paths] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','names');

    mainpath = '/Volumes/dumbo/stephj/TPM data analysis/SJLacI';

else
    [lengths,paths] = LoadDataAnalysis('HGseqsWProm','names');

    mainpath = '/Volumes/dumbo/stephj/TPM data analysis/SJLacI/HGseqs (Oid-E,T 80 to 90-O2 w prom)';
    
end

NLlensE = cell(length(lengths{1}),1);
UlensE =  cell(length(lengths{1}),1);
MlensE =  cell(length(lengths{1}),1);
BlensE =  cell(length(lengths{1}),1);
LlensE =  cell(length(lengths{1}),1);
UlensrelE =  cell(length(lengths{1}),1);
MlensrelE = cell(length(lengths{1}),1);
BlensrelE = cell(length(lengths{1}),1);
LlensrelE = cell(length(lengths{1}),1);

avgNLE = zeros(length(lengths{1}),1);
avgUE = zeros(length(lengths{1}),1);
avgME = zeros(length(lengths{1}),1);
avgBE = zeros(length(lengths{1}),1);
avgLE = zeros(length(lengths{1}),1);
avgUrelE = zeros(length(lengths{1}),1);
avgMrelE = zeros(length(lengths{1}),1);
avgBrelE = zeros(length(lengths{1}),1);
avgLrelE = zeros(length(lengths{1}),1);

SENLE = zeros(length(lengths{1}),1);
SEUE = zeros(length(lengths{1}),1);
SEME = zeros(length(lengths{1}),1);
SEBE = zeros(length(lengths{1}),1);
SELE = zeros(length(lengths{1}),1);
SEUrelE = zeros(length(lengths{1}),1);
SEMrelE = zeros(length(lengths{1}),1);
SEBrelE = zeros(length(lengths{1}),1);
SELrelE = zeros(length(lengths{1}),1);

%%
if Scatter
    figure
    hold on
end

for i = 1:length(lengths{1})
    [NLlensE{i}, UlensE{i},MlensE{i},BlensE{i},LlensE{i}] = CalculateTetherLength_Thresh(fullfile(mainpath,paths{1}{i}));
    UlensrelE{i} = abs(NLlensE{i}-UlensE{i}); 
    MlensrelE{i} = abs(NLlensE{i}-MlensE{i}); 
    BlensrelE{i} = abs(NLlensE{i}-BlensE{i}); 
    LlensrelE{i} = abs(NLlensE{i}-LlensE{i}); 
    
    avgNLE(i) = mean(NLlensE{i});
    SENLE(i) = std(NLlensE{i})/sqrt(length(NLlensE{i})-1);
    avgUE(i) = sum(UlensE{i})/length(find(UlensE{i})); %Don't want to include any zeros, because those indicate no U state for that bead
    SEUE(i) = std(UlensE{i}(find(UlensE{i})))/sqrt(length(find(UlensE{i}))-1);
    avgME(i) = sum(MlensE{i})/length(find(MlensE{i}));
    SEME(i) = std(MlensE{i}(find(MlensE{i})))/sqrt(length(find(MlensE{i}))-1);
    avgBE(i) = sum(BlensE{i})/length(find(BlensE{i}));
    SEBE(i) = std(BlensE{i}(find(BlensE{i})))/sqrt(length(find(BlensE{i}))-1);
    avgLE(i) = sum(LlensE{i})/length(find(LlensE{i}));
    SELE(i) = std(LlensE{i}(find(LlensE{i})))/sqrt(length(find(LlensE{i}))-1);
    
    %For the relative lengths, if the M or B entry for a bead is zero, then
    %the relative M or B will be the NL value for that bead.  So, can't
    %just take the average of UlensrelE, etc.  Those entries that I don't
    %want are well-separated from the real ones; so two ways to do this:
    %either use find to pick out all the ones less than, say, 70 nm; or do
    %it the right way and use find and logical indexing to pick out entries
    %that differ from that bead's NL length. Using the latter.
    avgUrelE(i) = mean(UlensrelE{i}(UlensrelE{i}~=NLlensE{i}));
    SEUrelE(i) = std(UlensrelE{i}(UlensrelE{i}~=NLlensE{i}))/sqrt(length(UlensrelE{i}(UlensrelE{i}~=NLlensE{i}))-1);
    avgMrelE(i) = mean(MlensrelE{i}(MlensrelE{i}~=NLlensE{i}));
    SEMrelE(i) = std(MlensrelE{i}(MlensrelE{i}~=NLlensE{i}))/sqrt(length(MlensrelE{i}(MlensrelE{i}~=NLlensE{i}))-1);
    avgBrelE(i) = mean(BlensrelE{i}(BlensrelE{i}~=NLlensE{i}));
    SEBrelE(i) = std(BlensrelE{i}(BlensrelE{i}~=NLlensE{i}))/sqrt(length(BlensrelE{i}(BlensrelE{i}~=NLlensE{i}))-1);
    avgLrelE(i) = mean(LlensrelE{i}(LlensrelE{i}~=NLlensE{i}));
    SELrelE(i) = std(LlensrelE{i}(LlensrelE{i}~=NLlensE{i}))/sqrt(length(LlensrelE{i}(LlensrelE{i}~=NLlensE{i}))-1);
    
%     if i==1
%         HandleERel = plot(lengths{1}(i).*ones(length(UlensrelE{i}),1),UlensrelE{i},'.k',...
%             lengths{1}(i).*ones(length(MlensrelE{i}),1),MlensrelE{i},'.r',...
%             lengths{1}(i).*ones(length(BlensrelE{i}),1),BlensrelE{i},'.b');
%     else
%         HandleERel(end+1:end+3) = plot(lengths{1}(i).*ones(length(UlensrelE{i}),1),UlensrelE{i},'.k',...
%             lengths{1}(i).*ones(length(MlensrelE{i}),1),MlensrelE{i},'.r',...
%             lengths{1}(i).*ones(length(BlensrelE{i}),1),BlensrelE{i},'.b');
%     end
    
    %better way:
    if Scatter
        if i==1
            HandleERel = plot(lengths{1}(i).*ones(length(UlensrelE{i}(UlensrelE{i}~=NLlensE{i})),1),...
                UlensrelE{i}(UlensrelE{i}~=NLlensE{i}),'.k',...
                lengths{1}(i).*ones(length(MlensrelE{i}(MlensrelE{i}~=NLlensE{i})),1),MlensrelE{i}(MlensrelE{i}~=NLlensE{i}),'.r',...
                lengths{1}(i).*ones(length(BlensrelE{i}(BlensrelE{i}~=NLlensE{i})),1),BlensrelE{i}(BlensrelE{i}~=NLlensE{i}),'.b');
        else
            HandleERel(end+1:end+3) = plot(lengths{1}(i).*ones(length(UlensrelE{i}(UlensrelE{i}~=NLlensE{i})),1),UlensrelE{i}(UlensrelE{i}~=NLlensE{i}),'.k',...
                lengths{1}(i).*ones(length(MlensrelE{i}(MlensrelE{i}~=NLlensE{i})),1),MlensrelE{i}(MlensrelE{i}~=NLlensE{i}),'.r',...
                lengths{1}(i).*ones(length(BlensrelE{i}(BlensrelE{i}~=NLlensE{i})),1),BlensrelE{i}(BlensrelE{i}~=NLlensE{i}),'.b');
        end

    end
end

if Scatter
    if ~HGseqs
        xlim([88 117])
    else
        xlim([55 89])
    end
    ylim([0 70])
    xlabel('Loop length (bp)')
    ylabel('RMS of state relative to No Lac (nm)')
    legend('Unlooped','Middle','Bottom')
    StandardFigure(HandleERel,gca)
end

NLlensT = cell(length(lengths{2}),1);
UlensT =  cell(length(lengths{2}),1);
MlensT =  cell(length(lengths{2}),1);
BlensT =  cell(length(lengths{2}),1);
LlensT =  cell(length(lengths{2}),1);
UlensrelT =  cell(length(lengths{2}),1);
MlensrelT = cell(length(lengths{2}),1);
BlensrelT = cell(length(lengths{2}),1);
LlensrelT = cell(length(lengths{2}),1);

avgNLT = zeros(length(lengths{2}),1);
avgUT = zeros(length(lengths{2}),1);
avgMT = zeros(length(lengths{2}),1);
avgBT = zeros(length(lengths{2}),1);
avgLT = zeros(length(lengths{2}),1);
avgUrelT = zeros(length(lengths{2}),1);
avgMrelT = zeros(length(lengths{2}),1);
avgBrelT = zeros(length(lengths{2}),1);
avgLrelT = zeros(length(lengths{2}),1);

SENLT = zeros(length(lengths{2}),1);
SEUT = zeros(length(lengths{2}),1);
SEMT = zeros(length(lengths{2}),1);
SEBT = zeros(length(lengths{2}),1);
SELT = zeros(length(lengths{2}),1);
SEUrelT = zeros(length(lengths{2}),1);
SEMrelT = zeros(length(lengths{2}),1);
SEBrelT = zeros(length(lengths{2}),1);
SELrelT = zeros(length(lengths{2}),1);

if Scatter
    figure
    hold on
end

for j = 1:length(lengths{2})
    [NLlensT{j}, UlensT{j},MlensT{j},BlensT{j},LlensT{j}] = CalculateTetherLength_Thresh(fullfile(mainpath,paths{2}{j}));
    UlensrelT{j} = abs(NLlensT{j}-UlensT{j}); 
    MlensrelT{j} = abs(NLlensT{j}-MlensT{j}); 
    BlensrelT{j} = abs(NLlensT{j}-BlensT{j}); 
    LlensrelT{j} = abs(NLlensT{j}-LlensT{j}); 
    
    avgNLT(j) = mean(NLlensT{j});
    SENLT(j) = std(NLlensT{j})/sqrt(length(NLlensT{j})-1);
    avgUT(j) = sum(UlensT{j})/length(find(UlensT{j})); %Don't want to include any zeros, because those indicate no U state for that bead
    SEUT(j) = std(UlensT{j}(find(UlensT{j})))/sqrt(length(find(UlensT{j}))-1);
    avgMT(j) = sum(MlensT{j})/length(find(MlensT{j}));
    SEMT(j) = std(MlensT{j}(find(MlensT{j})))/sqrt(length(find(MlensT{j}))-1);
    avgBT(j) = sum(BlensT{j})/length(find(BlensT{j}));
    SEBT(j) = std(BlensT{j}(find(BlensT{j})))/sqrt(length(find(BlensT{j}))-1);
    avgLT(j) = sum(LlensT{j})/length(find(LlensT{j}));
    SELT(j) = std(LlensT{j}(find(LlensT{j})))/sqrt(length(find(LlensT{j}))-1);
    
    avgUrelT(j) = mean(UlensrelT{j}(UlensrelT{j}~=NLlensT{j}));
    SEUrelT(j) = std(UlensrelT{j}(UlensrelT{j}~=NLlensT{j}))/sqrt(length(UlensrelT{j}(UlensrelT{j}~=NLlensT{j}))-1);
    avgMrelT(j) = mean(MlensrelT{j}(MlensrelT{j}~=NLlensT{j}));
    SEMrelT(j) = std(MlensrelT{j}(MlensrelT{j}~=NLlensT{j}))/sqrt(length(MlensrelT{j}(MlensrelT{j}~=NLlensT{j}))-1);
    avgBrelT(j) = mean(BlensrelT{j}(BlensrelT{j}~=NLlensT{j}));
    SEBrelT(j) = std(BlensrelT{j}(BlensrelT{j}~=NLlensT{j}))/sqrt(length(BlensrelT{j}(BlensrelT{j}~=NLlensT{j}))-1);
    avgLrelT(j) = mean(LlensrelT{j}(LlensrelT{j}~=NLlensT{j}));
    SELrelT(j) = std(LlensrelT{j}(LlensrelT{j}~=NLlensT{j}))/sqrt(length(LlensrelT{j}(LlensrelT{j}~=NLlensT{j}))-1);
    
    %warning, this throws an error because no bottom state for T95! I've
    %been doing the for-loop manually so can say HandleTRel(end+1:end+2)
    %for j=6
    if Scatter
        if j==1
            HandleTRel = plot(lengths{2}(j).*ones(length(UlensrelT{j}(UlensrelT{j}~=NLlensT{j})),1),...
                UlensrelT{j}(UlensrelT{j}~=NLlensT{j}),'.k',...
                lengths{2}(j).*ones(length(MlensrelT{j}(MlensrelT{j}~=NLlensT{j})),1),MlensrelT{j}(MlensrelT{j}~=NLlensT{j}),'.r',...
                lengths{2}(j).*ones(length(BlensrelT{j}(BlensrelT{j}~=NLlensT{j})),1),BlensrelT{j}(BlensrelT{j}~=NLlensT{j}),'.b');
        else
            HandleTRel(end+1:end+3) = plot(lengths{2}(j).*ones(length(UlensrelT{j}(UlensrelT{j}~=NLlensT{j})),1),...
                UlensrelT{j}(UlensrelT{j}~=NLlensT{j}),'.k',...
                lengths{2}(j).*ones(length(MlensrelT{j}(MlensrelT{j}~=NLlensT{j})),1),MlensrelT{j}(MlensrelT{j}~=NLlensT{j}),'.r',...
                lengths{2}(j).*ones(length(BlensrelT{j}(BlensrelT{j}~=NLlensT{j})),1),BlensrelT{j}(BlensrelT{j}~=NLlensT{j}),'.b');
        end
    end
    
end
% 
if Scatter
    if ~HGseqs
        xlim([88 117])
    else
        xlim([55 89])
    end
    ylim([0 70])
    xlabel('Loop length (bp)')
    ylabel('RMS of state relative to No Lac (nm)')
    legend('Unlooped','Middle','Bottom')
    StandardFigure(HandleTRel,gca)
end

%save('/Users/Steph/Desktop/110809lengthstemp.mat')

if HGseqs
    xlens = lengths{1}+36.*ones(1,length(lengths{1}));
else
    xlens = lengths{1};
end

figure
HandleAvgRelE = plot(xlens,0-avgUrelE,'.k',...
    xlens,0-avgMrelE,'.r',...
    xlens,0-avgBrelE,'.b');
hold on
HandleAvgRelE(end+1:end+length(lengths{1})+1) = errorbarxyHG(xlens,...
    0-avgUrelE,[],SEUrelE,[],[],'.k','k');
HandleAvgRelE(end+1:end+length(lengths{1})+1) = errorbarxyHG(xlens,...
    0-avgMrelE,[],SEMrelE,[],[],'.r','r');
HandleAvgRelE(end+1:end+length(lengths{1})+1) = errorbarxyHG(xlens,...
    0-avgBrelE,[],SEBrelE,[],[],'.b','b');
if ~HGseqs
    xlim([88 117])
else
    xlim([55+36 89+36])
end
ylim([-50 0])
xlabel('Loop length (bp)')
ylabel('RMS of state relative to No Lac (nm)')
legend('Unlooped','Middle','Bottom')
StandardFigure(HandleAvgRelE,gca)

if HGseqs
    xlens2 = lengths{2}+36.*ones(1,length(lengths{2}));
else
    xlens2 = lengths{2};
end

figure
HandleAvgRelT = plot(xlens2,0-avgUrelT,'.k',...
    xlens2,0-avgMrelT,'.r',...
    xlens2,0-avgBrelT,'.b');
hold on
HandleAvgRelT(end+1:end+length(lengths{2})+1) = errorbarxyHG(xlens2,...
    0-avgUrelT,[],SEUrelT,[],[],'.k','k');
HandleAvgRelT(end+1:end+length(lengths{2})+1) = errorbarxyHG(xlens2,...
    0-avgMrelT,[],SEMrelT,[],[],'.r','r');
HandleAvgRelT(end+1:end+length(lengths{2})+1) = errorbarxyHG(xlens2,...
    0-avgBrelT,[],SEBrelT,[],[],'.b','b');
if ~HGseqs
    xlim([88 117])
else
    xlim([55+36 89+36])
end
ylim([-50 0])
xlabel('Loop length (bp)')
ylabel('RMS of state relative to No Lac (nm)')
legend('Unlooped','Middle','Bottom')
StandardFigure(HandleAvgRelT,gca)

%%
%Concentration curves, Gauss
%I want 3 plots, one for unlooped, one for bottom, one for middle state,
%where the tether lengths of these states for each of the following data
%sets are plotted vs concentration on the same axes:
%Oid-E894-O1
%O1-E894-O1
%O2-E894-O1
%Oid-TA94-O1
%Oid-E8107-O1 -- This will need to be analyzed with the Thresh version in the above cells
%noOid-E872-O1
%noOid-E872-noO1

%For the operator KO controls, this analysis has already been done and is
%unaffected by the code errors that means I'll have to redo the rest:
[concsnOidE72nO1,lensnOidE72nO1] = LoadDataAnalysis('noOidE872noO1','Loop Lengths');
[concsnOidE72O1,lensnOidE72O1] = LoadDataAnalysis('noOidE872O1','Loop Lengths');

NLlensNoOps{1} = lensnOidE72nO1{1}.lengths.nolac;
NLlensNoOps{2} = lensnOidE72nO1{2}.lengths.nolac;
UlensNoOps{1} = lensnOidE72nO1{1}.lengths.unlooped;
UlensNoOps{2} = lensnOidE72nO1{2}.lengths.unlooped;
UlensrelNoOps{1} = lensnOidE72nO1{1}.lengths.NoLacMinusU;
UlensrelNoOps{2} = lensnOidE72nO1{2}.lengths.NoLacMinusU;

NLlensNoOid{1} = lensnOidE72O1{1}.lengths.nolac;
NLlensNoOid{2} = lensnOidE72O1{2}.lengths.nolac;
UlensNoOid{1} = lensnOidE72O1{1}.lengths.unlooped;
UlensNoOid{2} = lensnOidE72O1{2}.lengths.unlooped;
UlensrelNoOid{1} = lensnOidE72O1{1}.lengths.NoLacMinusU;
UlensrelNoOid{2} = lensnOidE72O1{2}.lengths.NoLacMinusU;

%Compute averages from the above cell arrays

for i = 1:length(concsnOidE72nO1)
    avgNLNoOps(i) = mean(NLlensNoOps{i});
    SENLNoOps(i) = std(NLlensNoOps{i})/sqrt(length(NLlensNoOps{i})-1);
    avgUNoOps(i) = sum(UlensNoOps{i})/length(find(UlensNoOps{i})); %Don't want to include any zeros, because those indicate no U state for that bead
    SEUNoOps(i) = std(UlensNoOps{i}(find(UlensNoOps{i})))/sqrt(length(find(UlensNoOps{i}))-1);
    
    %For the relative lengths, if the M or B entry for a bead is zero, then
    %the relative M or B will be the NL value for that bead.  So, can't
    %just take the average of UlensrelE, etc.  Those entries that I don't
    %want are well-separated from the real ones; so two ways to do this:
    %either use find to pick out all the ones less than, say, 70 nm; or do
    %it the right way and use find and logical indexing to pick out entries
    %that differ from that bead's NL length. Using the latter.
    avgUrelNoOps(i) = mean(UlensrelNoOps{i}(UlensrelNoOps{i}~=NLlensNoOps{i}));
    SEUrelNoOps(i) = std(UlensrelNoOps{i}(UlensrelNoOps{i}~=NLlensNoOps{i}))/sqrt(length(UlensrelNoOps{i}(UlensrelNoOps{i}~=NLlensNoOps{i}))-1);      
end
clear i
for i = 1:length(concsnOidE72O1)
    avgNLNoOid(i) = mean(NLlensNoOid{i});
    SENLNoOid(i) = std(NLlensNoOid{i})/sqrt(length(NLlensNoOid{i})-1);
    avgUNoOid(i) = sum(UlensNoOid{i})/length(find(UlensNoOid{i})); %Don't want to include any zeros, because those indicate no U state for that bead
    SEUNoOid(i) = std(UlensNoOid{i}(find(UlensNoOid{i})))/sqrt(length(find(UlensNoOid{i}))-1);
    
    %For the relative lengths, if the M or B entry for a bead is zero, then
    %the relative M or B will be the NL value for that bead.  So, can't
    %just take the average of UlensrelE, etc.  Those entries that I don't
    %want are well-separated from the real ones; so two ways to do this:
    %either use find to pick out all the ones less than, say, 70 nm; or do
    %it the right way and use find and logical indexing to pick out entries
    %that differ from that bead's NL length. Using the latter.
    avgUrelNoOid(i) = mean(UlensrelNoOid{i}(UlensrelNoOid{i}~=NLlensNoOid{i}));
    SEUrelNoOid(i) = std(UlensrelNoOid{i}(UlensrelNoOid{i}~=NLlensNoOid{i}))/sqrt(length(UlensrelNoOid{i}(UlensrelNoOid{i}~=NLlensNoOid{i}))-1);      
end
clear i

%Get lengths for E8107 conc curve--shoot all this time I've been using
%E108!  Won't change variable names but this is E10*7*
%[E108concs,E108paths] = LoadDataAnalysis('E8108conccurve','names');
[E108concs,E108paths] = LoadDataAnalysis('E8107conccurve','names');

mainpath = '/Volumes/dumbo/stephj/TPM data analysis/SJLacI/Oid-E,T 89 to 100-O1';

%Made a function to contain what's commented out below this call:
[NLlensE108, UlensE108,MlensE108,BlensE108,LlensE108,UlensrelE108,MlensrelE108,...
    BlensrelE108,LlensrelE108,avgNLE108,avgUE108,avgME108,avgBE108,avgLE108,...
    SENLE108,SEUE108,SEME108,SEBE108,SELE108,avgUrelE108,avgMrelE108,avgBrelE108,...
    avgLrelE108,SEUrelE108,SEMrelE108,SEBrelE108,SELrelE108] = GetAvgLengths(mainpath,E108paths,E108concs,0);

% NLlensE108 = cell(length(E108concs),1);
% UlensE108 =  cell(length(E108concs),1);
% MlensE108 =  cell(length(E108concs),1);
% BlensE108 =  cell(length(E108concs),1);
% LlensE108 =  cell(length(E108concs),1);
% UlensrelE108 =  cell(length(E108concs),1);
% MlensrelE108 = cell(length(E108concs),1);
% BlensrelE108 = cell(length(E108concs),1);
% LlensrelE108 = cell(length(E108concs),1);
% 
% avgNLE108 = zeros(length(E108concs),1);
% avgUE108 = zeros(length(E108concs),1);
% avgME108 = zeros(length(E108concs),1);
% avgBE108 = zeros(length(E108concs),1);
% avgLE108 = zeros(length(E108concs),1);
% avgUrelE108 = zeros(length(E108concs),1);
% avgMrelE108 = zeros(length(E108concs),1);
% avgBrelE108 = zeros(length(E108concs),1);
% avgLrelE108 = zeros(length(E108concs),1);
% 
% SENLE108 = zeros(length(E108concs),1);
% SEUE108 = zeros(length(E108concs),1);
% SEME108 = zeros(length(E108concs),1);
% SEBE108 = zeros(length(E108concs),1);
% SELE108 = zeros(length(E108concs),1);
% SEUrelE108 = zeros(length(E108concs),1);
% SEMrelE108 = zeros(length(E108concs),1);
% SEBrelE108 = zeros(length(E108concs),1);
% SELrelE108 = zeros(length(E108concs),1);
% 
% figure, hold on
% 
% for i = 1:length(E108concs)
%     [NLlensE108{i}, UlensE108{i},MlensE108{i},BlensE108{i},LlensE108{i}] = CalculateTetherLength_Thresh(fullfile(mainpath,E108paths{i}));
%     UlensrelE108{i} = abs(NLlensE108{i}-UlensE108{i}); 
%     MlensrelE108{i} = abs(NLlensE108{i}-MlensE108{i}); 
%     BlensrelE108{i} = abs(NLlensE108{i}-BlensE108{i}); 
%     LlensrelE108{i} = abs(NLlensE108{i}-LlensE108{i}); 
%     
%     avgNLE108(i) = mean(NLlensE108{i});
%     SENLE108(i) = std(NLlensE108{i})/sqrt(length(NLlensE108{i})-1);
%     avgUE108(i) = sum(UlensE108{i})/length(find(UlensE108{i})); %Don't want to include any zeros, because those indicate no U state for that bead
%     SEUE108(i) = std(UlensE108{i}(find(UlensE108{i})))/sqrt(length(find(UlensE108{i}))-1);
%     avgME108(i) = sum(MlensE108{i})/length(find(MlensE108{i}));
%     SEME108(i) = std(MlensE108{i}(find(MlensE108{i})))/sqrt(length(find(MlensE108{i}))-1);
%     avgBE108(i) = sum(BlensE108{i})/length(find(BlensE108{i}));
%     SEBE108(i) = std(BlensE108{i}(find(BlensE108{i})))/sqrt(length(find(BlensE108{i}))-1);
%     avgLE(i) = sum(LlensE108{i})/length(find(LlensE108{i}));
%     SELE(i) = std(LlensE108{i}(find(LlensE108{i})))/sqrt(length(find(LlensE108{i}))-1);
%     
%     %For the relative lengths, if the M or B entry for a bead is zero, then
%     %the relative M or B will be the NL value for that bead.  So, can't
%     %just take the average of UlensrelE, etc.  Those entries that I don't
%     %want are well-separated from the real ones; so two ways to do this:
%     %either use find to pick out all the ones less than, say, 70 nm; or do
%     %it the right way and use find and logical indexing to pick out entries
%     %that differ from that bead's NL length. Using the latter.
%     avgUrelE108(i) = mean(UlensrelE108{i}(UlensrelE108{i}~=NLlensE108{i}));
%     SEUrelE108(i) = std(UlensrelE108{i}(UlensrelE108{i}~=NLlensE108{i}))/sqrt(length(UlensrelE108{i}(UlensrelE108{i}~=NLlensE108{i}))-1);
%     avgMrelE108(i) = mean(MlensrelE108{i}(MlensrelE108{i}~=NLlensE108{i}));
%     SEMrelE108(i) = std(MlensrelE108{i}(MlensrelE108{i}~=NLlensE108{i}))/sqrt(length(MlensrelE108{i}(MlensrelE108{i}~=NLlensE108{i}))-1);
%     avgBrelE108(i) = mean(BlensrelE108{i}(BlensrelE108{i}~=NLlensE108{i}));
%     SEBrelE108(i) = std(BlensrelE108{i}(BlensrelE108{i}~=NLlensE108{i}))/sqrt(length(BlensrelE108{i}(BlensrelE108{i}~=NLlensE108{i}))-1);
%     avgLrelE108(i) = mean(LlensrelE108{i}(LlensrelE108{i}~=NLlensE108{i}));
%     SELrelE108(i) = std(LlensrelE108{i}(LlensrelE108{i}~=NLlensE108{i}))/sqrt(length(LlensrelE108{i}(LlensrelE108{i}~=NLlensE108{i}))-1);
%     
%     %To check state assignments:
% %     if i==1
% %          HandleERel = plot(E108concs(i).*10^12.*ones(length(UlensrelE108{i}(UlensrelE108{i}~=NLlensE108{i})),1),...
% %              UlensrelE108{i}(UlensrelE108{i}~=NLlensE108{i}),'.k',...
% %              E108concs(i).*10^12.*ones(length(MlensrelE108{i}(MlensrelE108{i}~=NLlensE108{i})),1),MlensrelE108{i}(MlensrelE108{i}~=NLlensE108{i}),'.r',...
% %              E108concs(i).*10^12.*ones(length(BlensrelE108{i}(BlensrelE108{i}~=NLlensE108{i})),1),BlensrelE108{i}(BlensrelE108{i}~=NLlensE108{i}),'.b');
% %      else
% %          HandleERel(end+1:end+3) = plot(E108concs(i).*10^12.*ones(length(UlensrelE108{i}(UlensrelE108{i}~=NLlensE108{i})),1),UlensrelE108{i}(UlensrelE108{i}~=NLlensE108{i}),'.k',...
% %              E108concs(i).*10^12.*ones(length(MlensrelE108{i}(MlensrelE108{i}~=NLlensE108{i})),1),MlensrelE108{i}(MlensrelE108{i}~=NLlensE108{i}),'.r',...
% %              E108concs(i).*10^12.*ones(length(BlensrelE108{i}(BlensrelE108{i}~=NLlensE108{i})),1),BlensrelE108{i}(BlensrelE108{i}~=NLlensE108{i}),'.b');
% %     end
% %     set(gca,'XScale','log')
% %     xlim([10^-15 10^-5])
% %     xlim([0.1 10000])
%     
% end

%The rest can be analyzed with the gauss version.
%Actually thresholding works much better, so thresholded all the
%concentration curve data:
[EOidO1concs,EOidO1paths] = LoadDataAnalysis('E894 conc curve','names');
[EO1O1concs,EO1O1paths] = LoadDataAnalysis('O1E894O1 conc curve','names');
[EO2O1concs,EO2O1paths] = LoadDataAnalysis('O2E894O1 conc curve','names');
[T94concs,T94paths] = LoadDataAnalysis('TA94 conc curve','names');

mainpathEOidO1 = '/Volumes/dumbo/stephj/TPM data analysis/SJLacI/E894 conc curve';
mainpathEO1O1 = '/Volumes/dumbo/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve';
mainpathEO2O1 = '/Volumes/dumbo/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve';
mainpathT94 = '/Volumes/dumbo/stephj/TPM data analysis/SJLacI/TA94 conc curve';

[NLlensEOidO1, UlensEOidO1,MlensEOidO1,BlensEOidO1,LlensEOidO1,UlensrelEOidO1,MlensrelEOidO1,...
    BlensrelEOidO1,LlensrelEOidO1,avgNLEOidO1,avgUEOidO1,avgMEOidO1,avgBEOidO1,avgLEOidO1,...
    SENLEOidO1,SEUEOidO1,SEMEOidO1,SEBEOidO1,SELEOidO1,avgUrelEOidO1,avgMrelEOidO1,avgBrelEOidO1,...
    avgLrelEOidO1,SEUrelEOidO1,SEMrelEOidO1,SEBrelEOidO1,SELrelEOidO1] = GetAvgLengths(mainpathEOidO1,EOidO1paths,EOidO1concs,0);

[NLlensEO1O1, UlensEO1O1,MlensEO1O1,BlensEO1O1,LlensEO1O1,UlensrelEO1O1,MlensrelEO1O1,...
    BlensrelEO1O1,LlensrelEO1O1,avgNLEO1O1,avgUEO1O1,avgMEO1O1,avgBEO1O1,avgLEO1O1,...
    SENLEO1O1,SEUEO1O1,SEMEO1O1,SEBEO1O1,SELEO1O1,avgUrelEO1O1,avgMrelEO1O1,avgBrelEO1O1,...
    avgLrelEO1O1,SEUrelEO1O1,SEMrelEO1O1,SEBrelEO1O1,SELrelEO1O1] = GetAvgLengths(mainpathEO1O1,EO1O1paths,EO1O1concs,0);

[NLlensEO2O1, UlensEO2O1,MlensEO2O1,BlensEO2O1,LlensEO2O1,UlensrelEO2O1,MlensrelEO2O1,...
    BlensrelEO2O1,LlensrelEO2O1,avgNLEO2O1,avgUEO2O1,avgMEO2O1,avgBEO2O1,avgLEO2O1,...
    SENLEO2O1,SEUEO2O1,SEMEO2O1,SEBEO2O1,SELEO2O1,avgUrelEO2O1,avgMrelEO2O1,avgBrelEO2O1,...
    avgLrelEO2O1,SEUrelEO2O1,SEMrelEO2O1,SEBrelEO2O1,SELrelEO2O1] = GetAvgLengths(mainpathEO2O1,EO2O1paths,EO2O1concs,0);

[NLlensT94, UlensT94,MlensT94,BlensT94,LlensT94,UlensrelT94,MlensrelT94,...
    BlensrelT94,LlensrelT94,avgNLT94,avgUT94,avgMT94,avgBT94,avgLT94,...
    SENLT94,SEUT94,SEMT94,SEBT94,SELT94,avgUrelT94,avgMrelT94,avgBrelT94,...
    avgLrelT94,SEUrelT94,SEMrelT94,SEBrelT94,SELrelT94] = GetAvgLengths(mainpathT94,T94paths,T94concs,0);


%Make the plots

%Relative Unlooped
plotconccurve({EOidO1concs,EO1O1concs,EO2O1concs,T94concs,E108concs,...
    concsnOidE72O1,concsnOidE72nO1},...
    {0-avgUrelEOidO1,0-avgUrelEO1O1,0-avgUrelEO2O1,0-avgUrelT94,...
    0-avgUrelE108,0-avgUrelNoOid,0-avgUrelNoOps},...
    {SEUrelEOidO1,SEUrelEO1O1,SEUrelEO2O1,SEUrelT94,SEUrelE108,...
    SEUrelNoOid,SEUrelNoOps},...
    {'Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1','Oid-E8107-O1','O1 only','No operators'},...
    [],[],{'.k','.b','.m','.r','.g','ok','xk'});
set(gca,'Ylim',[-6 1])
set(gca,'Xlim',[10*10^-15 10^-6])
ylabel('RMS relative to No Lac (nm)')

%Relative M and B
plotconccurve({EOidO1concs,EO1O1concs,EO2O1concs,T94concs,E108concs,E108concs},...
    {0-avgMrelEOidO1,0-avgMrelEO1O1,0-avgMrelEO2O1,0-avgMrelT94,0-avgMrelE108,0-avgBrelE108},...
    {SEMrelEOidO1,SEMrelEO1O1,SEMrelEO2O1,SEMrelT94,SEMrelE108,SEBrelE108},...
    {'Oid-E894-O1, M','O1-E894-O1, M','O2-E894-O1, M','Oid-TA94-O1, M','Oid-E8107-O1, M','Oid-E8107-O1, B'},...
    [],[],{'.k','.b','.m','.r','.g','og'});
set(gca,'Ylim',[-50 2]) %Rob and Martin thought it would make more sense to go negative
set(gca,'Xlim',[10*10^-15 10^-6])
ylabel('RMS relative to No Lac (nm)')

% mean(avgNLT94)
% temp=[];
% for x = 1:length(NLlensT94)
% temp(end+1:end+length(NLlensT94{x})) = NLlensT94{x};
% end
% std(temp)/sqrt(length(temp)-1)
