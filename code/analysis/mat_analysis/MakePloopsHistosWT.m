%(script) MakePloopsHistosWT
%
%This is an example of a script that can be used to generate plots of
%looping probability distributions, and that can save these distributions
%(and means and SEs, etc).  This example is for Steph's wild-type (WT) lac
%operon data.  
%
%It's best to run a script like this on a collection of data sets which 
%have some commonality (eg., a concentration curve with a single DNA), for 
%which you'll want to use the same rules to subtract nonloopers.
%
%Steph 3/2012

%Main path where all the subpath folders are
path = '/Volumes/dumbo/stephj/TPM data analysis/SJLacI/WTLac conc curves';

%Folders for individual data sets
subpath{1} = 'WTLac_noO3/WTLac_noO3_10pMSJLacI';
subpath{2} = 'WTLac_noO3/WTLac_noO3_100pMSJLacI';
subpath{3} = 'WTLac_noO3/WTLac_noO3_1nMSJLacI';
subpath{4} = 'WTLac_noO1_O3toO1/WTLac_noO1_O3toO1_10pMSJLacI';
subpath{5} = 'WTLac_noO1_O3toO1/WTLac_noO1_O3toO1_100pMSJLacI';
subpath{6} = 'WTLac_noO1_O3toO1/WTLac_noO1_O3toO1_1nMSJLacI';
subpath{7} = 'WTLac/WTLac_10pMSJLacI';
subpath{8} = 'WTLac/WTLac_100pMSJLacI';
subpath{9} = 'WTLac/WTLac_1nMSJLacI';

%Names: these are what each figure file will be called
names{1} = 'WTLacNoO310pM';
names{2} = 'WTLacNoO3100pM';
names{3} = 'WTLacNoO31nM';
names{4} = 'WTLacNoO1O3toO110pM';
names{5} = 'WTLacNoO1O3toO1100pM';
names{6} = 'WTLacNoO1O3toO11nM';
names{7} = 'WTLac10pM';
names{8} = 'WTLac100pM';
names{9} = 'WTLac1nM';

%%
%Do this the first time you calculate looping probabilities.

%Iterates through all the data sets, calculates distributions of looping
%probabilities keeping only those beads longer than 1000 sec, 2000 sec,
%etc, and extracts means, a Gaussian-fitted parameter (mostly useless),
%standard errors (SEs), etc. Saves figures and a mat file with the means
%etc in a new folder called "ThreshAnal" in each data folder.
for i=1:length(names)
    mkdir(fullfile(path,subpath{i},'ThreshAnal'))
    Gfit = load(fullfile(path,subpath{i},'Thresholding','ThreshFit.mat'));
    %If you want distributions of total pLoops only:
    [fighandle, Gaussmeans, GaussSEs, Means, SEs, ...
        alldistribs] = pLoopsvslengthDistrib(Gfit,names{i},fullfile(path,subpath{i},'ThreshAnal'));
    %If you want distributions of the two looped states separately as well:
    %note this does not make figures so you might want to run both the
    %previous line and this line every time.  The Gaussmeans, GaussSEs,
    %Means, SEs, and alldistribs from this are the same as the above.
    %Change the save command in the next line depending on which of these
    %you run!
    [notneeded, Gaussmeans, GaussSEs, Means, SEs, alldistribs,...
           MeansB,SEsB,alldistribsB,MeansM,SEsM,...
           alldistribsM] = pLoopsvslengthDistribBvsM(GfitsT{t},name,savenameT);
    %save(fullfile(path,subpath{i},'ThreshAnal',strcat(names{i},'ThreshAnal')),...
    %   'Gaussmeans','GaussSEs','Means','SEs','alldistribs')
    %If you run the DistribBvsM version,
    save(fullfile(path,subpath{i},'ThreshAnal',strcat(names{i},'ThreshAnal')),...
        'Gaussmeans','GaussSEs','Means','SEs','alldistribs','MeansB','SEsB',...
        'alldistribsB','MeansM','SEsM','alldistribsM')
    clear name Gaussmeans GaussSEs Means SEs alldistribs
       clear MeansB SEsB alldistribsB MeansM SEsM alldistribsM
       clear Gfit 
end

%%
%Do this any time you want to re-calculate means using different
%combinations of nonloopers.

%Use this for-loop to load the data saved by the above loop
for i = 1:length(names)
    temp = load(fullfile(path,subpath{i},'ThreshAnal',strcat(names{i},'ThreshAnal')));
    allalldistribs{i} = temp.alldistribs;
    allMeans{i} = temp.Means;
    %If you ran the DistribBvsM version in the previous cell, and want to
    %run the next cell that makes nice figures, you'll also need to load
    %MeansB, etc ... eg, allMeansB{i}=temp.MeansB;
end

%okrange are those data sets for which nonloopers are a sufficiently
%well-separated peak from the looping peak.  The numbers correspond to the
%numbers of the names array above: here we're subtracting all nonloopers
%(all beads with pLoop=0) from noO3_100pM and WT_100pM, and using the
%average value of the nonloopers from those two data sets to use as the
%fraction of nonloopers to subtract from everything else.
okrange = [2 8];
[alldistribssub02,Meanssub02,SEssub02,alldistribssub02B,Meanssub02B,SEssub02B,...
    alldistribssub02M,Meanssub02M,SEssub02M] = ComputeMeansMinusNonLoopers3BvsM(allalldistribs,...
    okrange,allalldistribsB,allalldistribsM);

save(fullfile(path,'WTLacThreshAnal'),'alldistribssub02','Meanssub02','SEssub02')

%%
%This last section makes prettier figures, using only beads that lasted
%>3000 sec and with nonloopers subtracted.  
for i = 1:length(names)
    pLoopsDistribBvsM(allalldistribs{i}{3}(:,2),allMeans{i}(3), Meanssub02{i}(3), ...
        alldistribsB{i}{3}(:,2),allMeansB{i}(3), Meanssub02B{i}(3), alldistribsM{i}{3}(:,2),...
        allMeansM{i}(3), Meanssub02M{i}(3), strcat(names{i},'ThreshAnal'),...
        fullfile(path,subpath{i},'ThreshAnal',strcat(names{i},'ThreshAnal')))
end
