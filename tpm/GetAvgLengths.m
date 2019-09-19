%function GetAvgLengths()
%
%This came up a lot in AnalyzeTetherLengths.m.  Mainpath is a path to the
%directory where the data folders are stored; for ex, mainpath = 
%'/Volumes/dumbo/stephj/TPM data analysis/SJLacI/Oid-E,T 89 to 100-O1'.
%Loadpaths is a cell array of folder names to load; for ex, the second
%output of LoadDataAnalysis when 'names' is passed into it. Concs is what the x-axis
%should be for the summary figure, and is a vector of the same length as loadpaths.
%Gauss is 1 if CalculateTetherLength_Gauss should be used, 0 if CalculateTetherLength_Thresh
%should be used.
%
%Outputs: the ones that don't say avg or SE are cell arrays of the same size as
%loadpaths.  The ones that say avg or SE are vectors of the same length as
%loadpaths.
%
%Steph 11/2011

function [NLlens, Ulens,Mlens,Blens,Llens,Ulensrel,Mlensrel,Blensrel,Llensrel,...
    avgNL,avgU,avgM,avgB,avgL,SENL,SEU,SEM,SEB,SEL,avgUrel,avgMrel,avgBrel,avgLrel,...
    SEUrel,SEMrel,SEBrel,SELrel] = GetAvgLengths(mainpath,loadpaths,concs,Gauss)

NLlens = cell(length(loadpaths),1);
Ulens =  cell(length(loadpaths),1);
Mlens =  cell(length(loadpaths),1);
Blens =  cell(length(loadpaths),1);
Llens =  cell(length(loadpaths),1);
Ulensrel =  cell(length(loadpaths),1);
Mlensrel = cell(length(loadpaths),1);
Blensrel = cell(length(loadpaths),1);
Llensrel = cell(length(loadpaths),1);

avgNL = zeros(length(loadpaths),1);
avgU = zeros(length(loadpaths),1);
avgM = zeros(length(loadpaths),1);
avgB = zeros(length(loadpaths),1);
avgL = zeros(length(loadpaths),1);
avgUrel = zeros(length(loadpaths),1);
avgMrel = zeros(length(loadpaths),1);
avgBrel = zeros(length(loadpaths),1);
avgLrel = zeros(length(loadpaths),1);

SENL = zeros(length(loadpaths),1);
SEU = zeros(length(loadpaths),1);
SEM = zeros(length(loadpaths),1);
SEB = zeros(length(loadpaths),1);
SEL = zeros(length(loadpaths),1);
SEUrel = zeros(length(loadpaths),1);
SEMrel = zeros(length(loadpaths),1);
SEBrel = zeros(length(loadpaths),1);
SELrel = zeros(length(loadpaths),1);

figure, hold on

for i = 1:length(loadpaths)
    if Gauss == 1
        [NLlens{i}, Ulens{i},Mlens{i},Blens{i},Llens{i}] = CalculateTetherLength_Gauss(fullfile(mainpath,loadpaths{i}));
    else
        [NLlens{i}, Ulens{i},Mlens{i},Blens{i},Llens{i}] = CalculateTetherLength_Thresh(fullfile(mainpath,loadpaths{i}));
    end
    
    Ulensrel{i} = abs(NLlens{i}-Ulens{i}); 
    Mlensrel{i} = abs(NLlens{i}-Mlens{i}); 
    Blensrel{i} = abs(NLlens{i}-Blens{i}); 
    Llensrel{i} = abs(NLlens{i}-Llens{i}); 
    
    avgNL(i) = mean(NLlens{i});
    SENL(i) = std(NLlens{i})/sqrt(length(NLlens{i})-1);
    avgU(i) = sum(Ulens{i})/length(find(Ulens{i})); %Don't want to include any zeros, because those indicate no U state for that bead
    SEU(i) = std(Ulens{i}(find(Ulens{i})))/sqrt(length(find(Ulens{i}))-1);
    avgM(i) = sum(Mlens{i})/length(find(Mlens{i}));
    SEM(i) = std(Mlens{i}(find(Mlens{i})))/sqrt(length(find(Mlens{i}))-1);
    avgB(i) = sum(Blens{i})/length(find(Blens{i}));
    SEB(i) = std(Blens{i}(find(Blens{i})))/sqrt(length(find(Blens{i}))-1);
    avgL(i) = sum(Llens{i})/length(find(Llens{i}));
    SEL(i) = std(Llens{i}(find(Llens{i})))/sqrt(length(find(Llens{i}))-1);
    
    %For the relative lengths, if the M or B entry for a bead is zero, then
    %the relative M or B will be the NL value for that bead.  So, can't
    %just take the average of UlensrelE, etc.  Those entries that I don't
    %want are well-separated from the real ones; so two ways to do this:
    %either use find to pick out all the ones less than, say, 70 nm; or do
    %it the right way and use find and logical indexing to pick out entries
    %that differ from that bead's NL length. Using the latter.
    avgUrel(i) = mean(Ulensrel{i}(Ulensrel{i}~=NLlens{i}));
    SEUrel(i) = std(Ulensrel{i}(Ulensrel{i}~=NLlens{i}))/sqrt(length(Ulensrel{i}(Ulensrel{i}~=NLlens{i}))-1);
    avgMrel(i) = mean(Mlensrel{i}(Mlensrel{i}~=NLlens{i}));
    SEMrel(i) = std(Mlensrel{i}(Mlensrel{i}~=NLlens{i}))/sqrt(length(Mlensrel{i}(Mlensrel{i}~=NLlens{i}))-1);
    avgBrel(i) = mean(Blensrel{i}(Blensrel{i}~=NLlens{i}));
    SEBrel(i) = std(Blensrel{i}(Blensrel{i}~=NLlens{i}))/sqrt(length(Blensrel{i}(Blensrel{i}~=NLlens{i}))-1);
    avgLrel(i) = mean(Llensrel{i}(Llensrel{i}~=NLlens{i}));
    SELrel(i) = std(Llensrel{i}(Llensrel{i}~=NLlens{i}))/sqrt(length(Llensrel{i}(Llensrel{i}~=NLlens{i}))-1);
    
    %To check state assignments:
%     if i==1
%          plot(concs(i).*ones(length(Ulensrel{i}(Ulensrel{i}~=NLlens{i})),1),...
%              Ulensrel{i}(Ulensrel{i}~=NLlens{i}),'.k',...
%              concs(i).*ones(length(Mlensrel{i}(Mlensrel{i}~=NLlens{i})),1),Mlensrel{i}(Mlensrel{i}~=NLlens{i}),'.r',...
%              concs(i).*ones(length(Blensrel{i}(Blensrel{i}~=NLlens{i})),1),Blensrel{i}(Blensrel{i}~=NLlens{i}),'.b');
%     else
% %         HandleERel(end+1:end+3) = plot(concs(i).*ones(length(Ulensrel{i}(Ulensrel{i}~=NLlens{i})),1),Ulensrel{i}(Ulensrel{i}~=NLlens{i}),'.k',...
% %              concs(i).*ones(length(Mlensrel{i}(Mlensrel{i}~=NLlens{i})),1),Mlensrel{i}(Mlensrel{i}~=NLlens{i}),'.r',...
% %              concs(i).*ones(length(Blensrel{i}(Blensrel{i}~=NLlens{i})),1),Blensrel{i}(Blensrel{i}~=NLlens{i}),'.b');
%         plot(concs(i).*ones(length(Ulensrel{i}(Ulensrel{i}~=NLlens{i})),1),Ulensrel{i}(Ulensrel{i}~=NLlens{i}),'.k',...
%              concs(i).*ones(length(Mlensrel{i}(Mlensrel{i}~=NLlens{i})),1),Mlensrel{i}(Mlensrel{i}~=NLlens{i}),'.r',...
%              concs(i).*ones(length(Blensrel{i}(Blensrel{i}~=NLlens{i})),1),Blensrel{i}(Blensrel{i}~=NLlens{i}),'.b');
%     end
    
    disp(strcat('Completed conc ',int2str(i),' of ',int2str(length(loadpaths))))
    
end

%This is clunky but this part needs to be changed manually depending on
%what's being plotted ...
%set(gca,'XScale','log')
%xlim([10^-15 10^-5])

% 
% xlim([88 117])
% ylim([0 70])
% xlabel('Loop length (bp)')
% ylabel('RMS of state relative to No Lac (nm)')
% legend('Unlooped','Middle','Bottom')
% StandardFigure(HandleTRel,gca)


