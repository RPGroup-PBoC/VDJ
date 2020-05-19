%function masterscriptV3(name,unloopedlength,varargin)
%
%Function that handles the analysis from the data saved by realtime, to
%looping probabilities for each bead.
%
%Inputs: 
%Name: name of the folder to be created where all of the analysis will be
%   saved.  If the folder already exists, you can save over old analysis or
%   if the analysis process was interrupted, can add to the saved analysis.
%Unloopedlength: length, in bp, of the tether (for plotting guides to the
%eye)
%Optional input: an already-saved GaussFit to edit
%
%Stephanie Johnson, 1/11

function masterscriptV3(name,unloopedlength,varargin)

%%% PARAMETERS %%%
fcB = 0.05; %The cutoff frequency for the Butterworth filter; we use 0.05 Hz
    %for 0.49 um diameter beads, or 0.07 Hz for the smaller 0.27 um
    %diameter beads
%fcB = 0.07
fcG = 0.0166; %8sec
%fcG = 0.0326;4sec
%The cutoff freq for Gaussian filtering; we use 0.0326 Hz, 
    %which corresponds to a standard deviation of sigma_c = 4 sec, according
    %to the formula sigma_c = sqrt(ln2)/(2*pi*fcG)); or for small beads,
    %fcG = corresponding to sigma_c = 2.87 sec.
%fcG = 0.0461

%%% ANALYSIS TYPE %%%
%Allows user to load pre-screened data and just perform a single-bead
%analysis
analtype=input('Start with raw data (1) or perform single-bead analysis (2)?');

%%% SINGLE-BEAD ANALYSIS %%%
if analtype==2
    storedir=uigetdir('/Volumes/dumbo/stephj/TPM data analysis','Choose folder with stored lacdataconcat.mat:');

    if exist(fullfile(storedir,'lacdataconcat.mat'))
        lacdata = load(fullfile(storedir,'lacdataconcat.mat'));
        nolacdata = load(fullfile(storedir,'nolacdata.mat'));
    else
        disp('Please select folder containing lacdataconcat.mat file.')
        return
    end

    if isfield(lacdata,'fpsL')
        fps = lacdata.fpsL{1}; %Have to assume all data was taken at the same frame rate, because all beads
            %get put into one cell array that is indifferent to set number
    else
        fps = input('Camera frame rate? (fps): ');
    end
    
    if nargin==2
        [GaussFit,type,pLoop,SEpLoop,rest] = singlebdanalysisV3(nolacdata.nolacANAL_RMS,nolacdata.nolacnames,...
            lacdata.newlacANAL_RMS,lacdata.corres3,lacdata.newlacnames,fps,storedir);
        
    else
        [GaussFit,type,pLoop,SEpLoop,rest] = singlebdanalysisV3(nolacdata.nolacANAL_RMS,nolacdata.nolacnames,...
            lacdata.newlacANAL_RMS,lacdata.corres3,lacdata.newlacnames,fps,storedir,varargin{1});
    end
                        
    close all
    
    mkdir(fullfile(storedir,'Single bead analysis'))

    %Save GaussFIt
    save(fullfile(storedir,'Single bead analysis','GaussFit'),'GaussFit')
    %Save stats
    if ~strcmpi(type,'WT')
        pB = rest{1};
        SEpB = rest{2};
        pM = rest{3};
        SEpM = rest{4};
        save(fullfile(storedir,'Single bead analysis','stats'),'pLoop','SEpLoop','pB','SEpB','pM','SEpM')
    else %means wt data i.e. 5 states
        pO23loop = rest{1};
        SEO23 = rest{2};
        pO12loop = rest{3};
        SEO12 = rest{4};
        pOtherloop = rest{5};
        SEotherloop = rest{6};
        save(fullfile(storedir,'Single bead analysis','stats'),'pLoop','SEpLoop','pO23loop',...
            'SEO23','pO12loop','SEO12','pOtherloop','SEotherloop')
    end

    makeallreptracesV3('Full Analysis Figures',storedir,lacdata.newlacANAL_RMS,lacdata.newlacnames,nolacdata.nolacANAL_RMS,...
        nolacdata.nolacnames,lacdata.corres3,GaussFit,fps,type)
        

else
%%% STARTING WITH RAW DATA %%%

    %%% LOAD FILES %%%
    %Load parameters that user has entered into masterscriptV2params (it didn't
    %change when we changed to masterscriptV3, and it contains a list of all
    %data analyzed, so keep the same file)
    [path, nolacnames, lacnames, lacnames2, corres] = masterscriptV2params();
    %Check that parameters make sense
    if length(lacnames2)>1 && length(lacnames)~=length(lacnames2)
        disp('Lacnames or Lacnames2 entered incorrectly.')
        return 
    end
    if length(corres)~=length(lacnames)
        disp('Corres or Lacnames entered incorrectly.')
        return
    end
    if length(nolacnames)~=length(lacnames) 
        disp('Number of nolac and lac sets not equal.')
        return
    end
    %Ask user where to create a folder for analyzed data
    savedir = uigetdir('/Volumes/dumbo/stephj/TPM data analysis','Save analyzed data where?');
    %Check if user will be saving over other data
    if exist(fullfile(savedir,name))
        saveover=input('Save directory already exists: save over existing directory? (y/n)','s');
        
        if strcmpi(saveover,'n') && exist(fullfile(savedir,name,'lacdata.mat'))
            saveover2=input('Continue with a previous analysis? (y/n)','s');
            if strcmpi(saveover2,'y')
                savedir = fullfile(savedir,name);
                tracesavedir = fullfile(savedir,strcat(name,' traces'));
                alllacnames = lacnames;
                %Resume an interrupted analysis, or add additional sets. Work
                %forwards to find out where the previous analysis ended.
                oldnolac = load(fullfile(savedir,'nolacdata.mat'));
                concat = 0;

                %If analysis was interrupted during analysis of first lac sets,
                %or user wants to add more sets: the best indicator of this is
                %the length of nolacANAL_RMS, since the length of lacANAL_RMS
                %is complicated by the possibility of second sets.
                if length(oldnolac.nolacANAL_RMS)<length(nolacnames)
                    concat = 1; %If go through either this or the next section, need to redo lacdataconcat and corres2 and corres3 files
                    startset = length(oldnolac.nolacANAL_RMS)+1;
                    for i=startset:length(nolacnames)
                        [lacXtemp,lacYtemp,fpsLtemp,success,wind,clA]=loadrawsV3(path, char(lacnames(i)), strcat('LacSet',int2str(i)));
                        if ~success
                            return
                        end
                        %Check that the corres matrix was entered correctly in
                        %masterscriptV2params:
                        if length(corres{i})~=size(lacXtemp,1)
                            disp(strcat('Corres matrix for set ',int2str(i),' entered incorrectly.'))
                            return
                        end

                        %Drift correct and Gaussian filter the first-set lac data:
                        [laccorrXtemp,laccorrYtemp, lacr1]= meandriftcorrfunctionV3(lacXtemp,lacYtemp,fpsLtemp,fcB);
                        %**
                        %To drift correct with Martin's code instead:
                        %laccorrXtemp=BWdriftcorrect_ML(lacXtemp',fcB,fpsLtemp);
                        %laccorrYtemp=BWdriftcorrect_ML(lacYtemp',fcB,fpsLtemp);
                        %**
                        lacRMS = sqrt(gaussfiltV3(lacr1.^2, fcG, fpsLtemp));
                        %**
                        %To Gaussian filter with Martin's code instead: Note his code takes
                        %the cutoff frequency in sec.
                        %[lacRMS,notneeded,notneeded,notneeded]=gaussFilter_ML([laccorrXtemp' laccorrYtemp'],4,fpsLtemp);
                        %**

                        %Load corresponding no-lac sets, discarding any beads in the no-lac sets
                        %for which no lac data was taken.  The resulting nolacX and nolacY 
                        %arrays have the same length as lacnames
                        [nolacXtemp2,nolacYtemp2,fpsNLtemp,success,notneeded,notneeded]=loadrawV3(path, char(nolacnames(i)), strcat('NoLacSet',int2str(i)));
                        if ~success
                            return
                        end
                        %nolacXtemp and nolacYtemp are bds-by-frames matrices.  That is, each
                        %row is a bead. Corres{i} is a vector where the value of each element
                        %is the number of a bead to keep.  So, keep all the rows whose number
                        %is in corres{i}:
                        nolacXtemp = nolacXtemp2(corres{i},:);
                        nolacYtemp = nolacYtemp2(corres{i},:);
                        clear nolacXtemp2 nolacYtemp2

                        %Drift correct and Gaussian filter the no lac data:
                        [nolaccorrXtemp,nolaccorrYtemp, nolacr1]= meandriftcorrfunctionV3(nolacXtemp,nolacYtemp,fpsNLtemp,fcB);
                        %To drift correct with Martin's code instead:
                        %nolaccorrXtemp=BWdriftcorrect_ML(nolacXtemp',fcB,fpsNLtemp);
                        %nolaccorrYtemp=BWdriftcorrect_ML(nolacYtemp',fcB,fpsNLtemp);
                        %**
                        nolacANAL_RMStemp = sqrt(gaussfiltV3(nolacr1.^2, fcG, fpsNLtemp));
                        %**
                        %To Gaussian filter with Martin's code instead:
                        %[nolacANAL_RMStemp,notneeded,notneeded,notneeded]=gaussFilter_ML([nolaccorrXtemp' nolaccorrYtemp'],4,fpsNLtemp);
                        %**

                        %Allow user to screen the data by watching movies:
                        pxlpath = fullfile(path,char(lacnames(i)));
                        [lacANAL_RMStemp,recordtemp,cutptstemp] = screenbeadsV3(nolacANAL_RMStemp,...
                                        lacRMS,fpsLtemp, fpsNLtemp,char(nolacnames(i)),char(lacnames(i)),...
                                        pxlpath,tracesavedir, unloopedlength,wind,clA);  %Note that this allows the lac and 
                                            %nolac sets to have different fps's.  That's
                                            %because old versions of the acquisition code have
                                            %a bug that does this in some cases.

                        %Update the lacdata and nolacdata files--done this way to save space
                        %for very large sets to analyze
                        oldlac = load(fullfile(savedir,'lacdata.mat'));
                        if i > startset %For the first set the previous nolac analysis is already loaded
                            oldnolac = load(fullfile(savedir,'nolacdata.mat'));
                            %oldlac = load(fullfile(savedir,'lacdata.mat'));
                        end
                        nolacX = [oldnolac.nolacX, nolacXtemp];
                        nolacY = [oldnolac.nolacY, nolacYtemp];
                        nolacANAL_RMS = [oldnolac.nolacANAL_RMS, nolacANAL_RMStemp];
                        nolaccorrX = [oldnolac.nolaccorrX, nolaccorrXtemp];
                        nolaccorrY = [oldnolac.nolaccorrY, nolaccorrYtemp];
                        fpsNL = [oldnolac.fpsNL fpsNLtemp];
                        %If the user is adding data, and there were originally
                        %some second sets, can't just concatenate these results
                        %on the ends of the original arrays!
                        if length(oldlac.lacANAL_RMS)>length(oldnolac.nolacANAL_RMS)
                            lacX = [oldlac.lacX(1:i-1), lacXtemp,oldlac.lacX(i:end)];
                            lacY = [oldlac.lacY(1:i-1), lacYtemp,oldlac.lacY(i:end)];
                            lacANAL_RMS = [oldlac.lacANAL_RMS(1:i-1), lacANAL_RMStemp,oldlac.lacANAL_RMS(i:end)];
                            record = [oldlac.record(1:i-1), recordtemp,oldlac.record(i:end)];
                            cutpts = [oldlac.cutpts(1:i-1), cutptstemp,oldlac.cutpts(i:end)];
                            laccorrX = [oldlac.laccorrX(1:i-1), laccorrXtemp,oldlac.laccorrX(i:end)];
                            laccorrY = [oldlac.laccorrY(1:i-1), laccorrYtemp,oldlac.laccorrY(i:end)];
                            fpsL = [oldlac.fpsL(1:i-1), fpsLtemp,oldlac.fpsL(i:end)];
                        else
                            lacX = [oldlac.lacX, lacXtemp];
                            lacY = [oldlac.lacY, lacYtemp];
                            lacANAL_RMS = [oldlac.lacANAL_RMS, lacANAL_RMStemp];
                            record = [oldlac.record, recordtemp];
                            cutpts = [oldlac.cutpts, cutptstemp];
                            laccorrX = [oldlac.laccorrX, laccorrXtemp];
                            laccorrY = [oldlac.laccorrY, laccorrYtemp];
                            fpsL = [oldlac.fpsL, fpsLtemp];
                        end

                        save(fullfile(savedir,'nolacdata'),'nolacANAL_RMS','nolacX', 'nolacY',...
                            'nolaccorrX', 'nolaccorrY', 'nolacnames','fpsNL');
                        save(fullfile(savedir,'lacdata'),'lacANAL_RMS', 'lacX', 'lacY', 'laccorrX',...
                            'laccorrY', 'alllacnames', 'corres', 'record', 'cutpts','fpsL')

                        disp(strcat('Lac Set ', int2str(i),'/',int2str(length(lacnames)),' Finished'))
                        close all
                        if i<length(nolacnames) %If this is the last time through this loop, don't delete these things so can call them below
                            clear lacX lacY laccorrX laccorrY nolacX  nolacY nolaccorrX nolaccorrY
                            clear nolacANAL_RMS lacANAL_RMS record cutpts fpsL fpsNL
                        end
                        clear fpsNLtemp fpsLtemp success lacr1 lacRMS nolacr1 pxlpath lacXtemp ...
                            lacYtemp laccorrXtemp laccorrYtemp nolacXtemp nolacYtemp nolaccorrXtemp nolaccorrYtemp...
                            nolacANAL_RMStemp lacANAL_RMStemp recordtemp cutptstemp oldlac oldnolac

                    end
                end
                
                if ~concat %if didn't go through the first if statement, need the information in lacdata.mat
                    oldlac = load(fullfile(savedir,'lacdata.mat'));
                    oldnolac = load(fullfile(savedir,'nolacdata.mat'));
                    lacANAL_RMS = oldlac.lacANAL_RMS; %Necessary for the second if statement, if the first one is skipped
                    record = oldlac.record;
                    alllacnames = oldlac.alllacnames;
                    nolacANAL_RMS = oldnolac.nolacANAL_RMS;
                    lacX = oldlac.lacX;
                    lacY = oldlac.lacY;
                    laccorrX = oldlac.laccorrX;
                    laccorrY = oldlac.laccorrY;
                    corres = oldlac.corres;
                    cutpts = oldlac.cutpts;
                    fpsL = oldlac.fpsL;
                    fpsNL = oldnolac.fpsNL;
                end

                %Second sets: Having the analysis interupted and adding second
                %sets look the same at this point: lacANAL_RMS will have fewer
                %cells than lacnames+lacnames2
                if (~(length(lacnames2)==1 && strcmpi(lacnames2{1},''))) && length(lacANAL_RMS)<(length(lacnames)+length(lacnames2))%Second set on at least one area, and analysis not finished
                    concat=1;
                    alllacnames = [lacnames lacnames2];
                    startset = length(lacANAL_RMS)-length(nolacANAL_RMS)+1;
                    for i=startset:length(lacnames2)
                        if strcmpi(lacnames2{i},'')
                            lacXtemp = [];
                            lacYtemp = [];
                            laccorrXtemp=[];
                            laccorrYtemp=[];
                            lacANAL_RMStemp=[];
                            recordtemp = [];
                            cutptstemp = [];
                            fpsLtemp = [];
                        else
                            [lacXtemp,lacYtemp,fpsLtemp,success,wind,clA]=loadrawV3(path, char(lacnames2(i)), strcat('Lac2Set',int2str(i)));
                            if ~success
                                return
                            end

                            %Drift correct and Gaussian filter the first-set lac data:
                            [laccorrXtemp,laccorrYtemp, lacr1]= meandriftcorrfunctionV3(lacXtemp,lacYtemp,fpsLtemp,fcB);
                            %**
                            %To drift correct with Martin's code instead:
                            %laccorrXtemp=BWdriftcorrect_ML(lacXtemp',fcB,fpsLtemp);
                            %laccorrYtemp=BWdriftcorrect_ML(lacYtemp',fcB,fpsLtemp);
                            %**
                            lacRMS = sqrt(gaussfiltV3(lacr1.^2, fcG, fpsLtemp));
                            %**
                            %To Gaussian filter with Martin's code instead:
                            %[lacRMS,notneeded,notneeded,notneeded]=gaussFilter_ML([laccorrXtemp' laccorrYtemp'],4,fpsLtemp);
                            %**

                            %Allow user to screen the data by watching movies:
                            pxlpath = fullfile(path,char(lacnames2(i)));
                            [lacANAL_RMStemp,recordtemp,cutptstemp] = screenbeadsV3(nolacANAL_RMS{i},...
                                            lacRMS,fpsLtemp, fpsNL{i},char(nolacnames(i)),char(lacnames2(i)),...
                                            pxlpath,tracesavedir, unloopedlength,wind,clA);  %Here we just assume the same fps for lac and nolac out of laziness ...

                        end

                        %Update the lacdata and nolacdata files--note that there won't ever
                        %not be a file already created for these second sets, but don't
                        %need to reload for the first lac2 set
                        if i>startset
                            oldlac = load(fullfile(savedir,'lacdata.mat'));
                            lacX = oldlac.lacX;
                            lacY = oldlac.lacY;
                            lacANAL_RMS = oldlac.lacANAL_RMS;
                            record = oldlac.record;
                            cutpts = oldlac.cutpts;
                            laccorrX = oldlac.laccorrX;
                            laccorrY = oldlac.laccorrY;
                            fpsL = oldlac.fpsL;
                        end
                        
                        lacX{end+1} = lacXtemp; %Can't concatenate empty matrices onto cells, so need to do this in case any of the lacX, etc are empty
                        lacY{end+1} = lacYtemp;
                        lacANAL_RMS{end+1} = lacANAL_RMStemp;
                        record{end+1} = recordtemp;
                        cutpts{end+1} = cutptstemp;
                        laccorrX{end+1} = laccorrXtemp;
                        laccorrY{end+1} = laccorrYtemp;
                        fpsL{end+1} = fpsLtemp;

                        save(fullfile(savedir,'lacdata'),'lacANAL_RMS', 'lacX', 'lacY', 'laccorrX',...
                                'laccorrY', 'alllacnames', 'corres', 'record', 'cutpts','fpsL')

                        disp(strcat('Lac 2 Set ', int2str(i),'/',int2str(length(lacnames2)),' Finished'))
                        close all   

                        if i<length(lacnames2) %Unless this is the last set, clear all variables
                            clear lacX lacY laccorrX laccorrY lacANAL_RMS record cutpts fpsL
                        end
                        clear fpsLtemp success lacr1 lacRMS pxlpath lacXtemp ...
                            lacYtemp laccorrXtemp laccorrYtemp lacANAL_RMStemp recordtemp cutptstemp oldlac

                    end
                end

                if ~exist(fullfile(savedir,'lacdataconcat.mat')) || concat
                    delete(fullfile(savedir,'lacdataconcat.mat'))

                    if ~concat %means no lacdataconcat file but had previously finished all first and second lac set analysis
                        record = oldlac.record;
                        alllacnames = oldlac.alllacnames;
                        nolacANAL_RMS = oldnolac.nolacANAL_RMS;
                        lacANAL_RMS = oldlac.lacANAL_RMS;
                        lacX = oldlac.lacX;
                        lacY = oldlac.lacY;
                        laccorrX = oldlac.laccorrX;
                        laccorrY = oldlac.laccorrY;
                        corres = oldlac.corres;
                        cutpts = oldlac.cutpts;
                        fpsL = oldlac.fpsL;
                    end

                    %Assume the final version of the lacdata.mat file, with
                    %corres2 and corres3, didn't get saved or needs to be
                    %updated.  (Rare to have the analysis interrupted between
                    %finishing all analysis and this final lacdata save).

                    %Create the corres2 and corres3 matrices.  Corres2
                    %relates the beads numbers for multiple sets on the same area and is used
                    %primarily for bdconcat; corres3 relates the bead numbers in lacdataconcat
                    %to those in the saved nolac data file.

                    [corres2, corres3] = createCorres2and3(lacnames,alllacnames,record,nolacANAL_RMS);

                    %Save the lacdata.mat file a final time, this time with corres2 and corres3
                    save(fullfile(savedir,'lacdata'),'lacANAL_RMS', 'lacX', 'lacY', 'laccorrX',...
                                    'laccorrY', 'alllacnames', 'corres', 'record', 'cutpts','fpsL',...
                                    'corres2','corres3')
                    disp('Lac Data Saved')

                    %We used to make a histogram of all the beads summed, but Steph never uses it.
                    %Use probhistfunctionV2 if necessary.

                    clear lacX lacY laccorrX lacorrY nolacX nolacY nolaccorrX nolaccorrY record cutpts
                    clear tracesavedir

                    %If there were second sets, concatenate data from different sets that are on the same beads:
                    if ~(isequal(length(lacnames),length(alllacnames)))
                        [newlacANAL_RMS,newlacnames]=bdconcatV3(lacANAL_RMS,corres2,alllacnames);
                    else
                        newlacANAL_RMS = lacANAL_RMS;
                        newlacnames = alllacnames;
                    end

                    save(fullfile(savedir,'lacdataconcat'),'newlacANAL_RMS', 'newlacnames', 'corres2','corres3','fpsL')
                    disp('Lacdataconcat saved')

                    concat = 1;
                end

                if ~exist(fullfile(savedir,'Single bead analysis')) || concat
                    if ~concat
                            lacdataconcat = load(fullfile(savedir,'lacdataconcat.mat'));
                            nolacANAL_RMS = oldnolac.nolacANAL_RMS;
                            newlacANAL_RMS = lacdataconcat.newlacANAL_RMS;
                            corres3=oldlac.corres3;
                            newlacnames = lacdataconcat.newlacnames;
                            fpsL = oldlac.fpsL;
                    end
                    
                    if exist(fullfile(savedir,'Single bead analysis'))
                        oldGfit = load(fullfile(savedir,'Single bead analysis','GaussFit.mat'));
                        oldGfit = oldGfit.GaussFit;
                        %Figure out which beads have already been fit so
                        %the user doesn't have to redo those
                        startset=1;
                        totbds = 0;
                        for i=1:length(newlacANAL_RMS)
                            if length(oldGfit)>=totbds+size(newlacANAL_RMS{i},1)
                                startset = startset+1;
                                totbds = totbds+size(newlacANAL_RMS{i},1);
                            end
                        end
                        [newGaussFit,type,pLoop,SEpLoop,rest] = singlebdanalysisV3(nolacANAL_RMS(startset:end),nolacnames(startset:end),newlacANAL_RMS(startset:end),...
                            corres3(startset:end),newlacnames(startset:end),fpsL{1},savedir);   
                        GaussFit = oldGfit;
                        for i=1:length(newGaussFit)
                            GaussFit(end+1) = newGaussFit(i);
                        end
                        
                    else
                        [GaussFit,type,pLoop,SEpLoop,rest] = singlebdanalysisV3(nolacANAL_RMS,nolacnames,newlacANAL_RMS,...
                            corres3,newlacnames,fpsL{1},savedir);
                    end
                    
                    close all
                    
                    mkdir(fullfile(savedir,'Single bead analysis'))
                    
                    %Save GaussFIt
                    save(fullfile(savedir,'Single bead analysis','GaussFit'),'GaussFit')
                    %Save stats
                    if ~strcmpi(type,'WT')
                        pB = rest{1};
                        SEpB = rest{2};
                        pM = rest{3};
                        SEpM = rest{4};
                        save(fullfile(savedir,'Single bead analysis','stats'),'pLoop','SEpLoop','pB','SEpB','pM','SEpM')
                    else %means wt data i.e. 5 states
                        pO23loop = rest{1};
                        SEO23 = rest{2};
                        pO12loop = rest{3};
                        SEO12 = rest{4};
                        pOtherloop = rest{5};
                        SEotherloop = rest{6};
                        save(fullfile(savedir,'Single bead analysis','stats'),'pLoop','SEpLoop','pO23loop',...
                            'SEO23','pO12loop','SEO12','pOtherloop','SEotherloop')
                    end
                    
                    %Give the option of creating nice, more final figures.
                    %In _V3, not going to bother to ask ...
                    %cont2=input('Make pretty figures? (y/n)','s');
                    %if strcmpi(cont2,'n')
                    %    return
                    %else
                    makeallreptracesV3('Full Analysis Figures',savedir,newlacANAL_RMS,newlacnames,nolacANAL_RMS,...
                        nolacnames,corres3,GaussFit,fpsL{1},type)
                    %end
                
                end
            else
                return
            end
            return
        end
        
    end
        
    %elseif ~exist(fullfile(savedir,name)) || strcmpi(saveover,'y')
        %Delete these files that get re-loaded when redoing the analysis.
        %This just produces a warning if they don't exist.
        delete(fullfile(savedir,name,'nolacdata.mat'))
        delete(fullfile(savedir,name,'lacdata.mat'))
        delete(fullfile(savedir,name,'lacdataconcat.mat')) 

        %Set up the folders for saving things:
        %Make the parent directory:
        mkdir(savedir,name);
        savedir = fullfile(savedir,name);
        mkdir(savedir,strcat(name,' traces'));
        tracesavedir = fullfile(savedir,strcat(name,' traces'));

        %Load each data set and analyze.  In contrast to previous versions, this
        %version executes a for-loop so that each set is loaded, analyzed, saved,
        %and removed from memory.  This means even if masterscriptV3 crashes or
        %analysis is executed, analysis of previous sets is saved; and this way
        %very large data sets can be analyzed.

        %Load and analyze all first-lac sets
        for i=1:length(lacnames)
            [lacXtemp,lacYtemp,fpsLtemp,success,wind,clA]=loadrawV3(path, char(lacnames(i)), strcat('LacSet',int2str(i)));
            if ~success
                return
            end
            %Check that the corres matrix was entered correctly in
            %masterscriptV2params:
            if length(corres{i})~=size(lacXtemp,1)
                disp(strcat('Corres matrix for set ',int2str(i),' entered incorrectly.'))
                return
            end

            %Drift correct and Gaussian filter the first-set lac data:
            [laccorrXtemp,laccorrYtemp, lacr1]= meandriftcorrfunctionV3(lacXtemp,lacYtemp,fpsLtemp,fcB);
            %**
            %To drift correct with Martin's code instead:
            %laccorrXtemp=BWdriftcorrect_ML(lacXtemp',fcB,fpsLtemp);
            %laccorrYtemp=BWdriftcorrect_ML(lacYtemp',fcB,fpsLtemp);
            %**
            lacRMS = sqrt(gaussfiltV3(lacr1.^2, fcG, fpsLtemp));
            %**
            %To Gaussian filter with Martin's code instead: Note his code takes
            %the cutoff frequency in sec.
            %[lacRMS,notneeded,notneeded,notneeded]=gaussFilter_ML([laccorrXtemp' laccorrYtemp'],4,fpsLtemp);
            %**

            %Load corresponding no-lac sets, discarding any beads in the no-lac sets
            %for which no lac data was taken.  The resulting nolacX and nolacY 
            %arrays have the same length as lacnames
            [nolacXtemp2,nolacYtemp2,fpsNLtemp,success,notneeded,notneeded]=loadrawV3(path, char(nolacnames(i)), strcat('NoLacSet',int2str(i)));
            if ~success
                return
            end
            %nolacXtemp and nolacYtemp are bds-by-frames matrices.  That is, each
            %row is a bead. Corres{i} is a vector where the value of each element
            %is the number of a bead to keep.  So, keep all the rows whose number
            %is in corres{i}:
            nolacXtemp = nolacXtemp2(corres{i},:);
            nolacYtemp = nolacYtemp2(corres{i},:);
            clear nolacXtemp2 nolacYtemp2

            %Drift correct and Gaussian filter the no lac data:
            [nolaccorrXtemp,nolaccorrYtemp, nolacr1]= meandriftcorrfunctionV3(nolacXtemp,nolacYtemp,fpsNLtemp,fcB);
            %To drift correct with Martin's code instead:
            %nolaccorrXtemp=BWdriftcorrect_ML(nolacXtemp',fcB,fpsNLtemp);
            %nolaccorrYtemp=BWdriftcorrect_ML(nolacYtemp',fcB,fpsNLtemp);
            %**
            nolacANAL_RMStemp = sqrt(gaussfiltV3(nolacr1.^2, fcG, fpsNLtemp));
            %**
            %To Gaussian filter with Martin's code instead:
            %[nolacANAL_RMStemp,notneeded,notneeded,notneeded]=gaussFilter_ML([nolaccorrXtemp' nolaccorrYtemp'],4,fpsNLtemp);
            %**

            %Allow user to screen the data by watching movies:
            pxlpath = fullfile(path,char(lacnames(i)));
            [lacANAL_RMStemp,recordtemp,cutptstemp] = screenbeadsV3(nolacANAL_RMStemp,...
                            lacRMS,fpsLtemp, fpsNLtemp,char(nolacnames(i)),char(lacnames(i)),...
                            pxlpath,tracesavedir, unloopedlength,wind,clA);  %Note that this allows the lac and 
                                %nolac sets to have different fps's.  That's
                                %because old versions of the acquisition code have
                                %a bug that does this in some cases.

            %Update the lacdata and nolacdata files--done this way to save space
            %for very large sets to analyze
            if ~exist(fullfile(savedir,'nolacdata.mat')) %Haven't already saved 1 or more sets
                nolacANAL_RMS{1} = nolacANAL_RMStemp;
                nolacX{1} = nolacXtemp;
                nolacY{1} = nolacYtemp;
                nolaccorrX{1} = nolaccorrXtemp;
                nolaccorrY{1} = nolaccorrYtemp;
                fpsNL{1} = fpsNLtemp;
                save(fullfile(savedir,'nolacdata'),'nolacANAL_RMS','nolacX', 'nolacY',...
                    'nolaccorrX', 'nolaccorrY', 'nolacnames','fpsNL');
                alllacnames = lacnames; %for consistency with previous versions of this code
                lacANAL_RMS{1}=lacANAL_RMStemp;
                lacX{1}=lacXtemp;
                lacY{1}=lacYtemp;
                laccorrX{1}=laccorrXtemp;
                laccorrY{1}=laccorrYtemp;
                record{1}=recordtemp;
                cutpts{1}=cutptstemp;
                fpsL{1}=fpsLtemp;
                save(fullfile(savedir,'lacdata'),'lacANAL_RMS', 'lacX', 'lacY', 'laccorrX',...
                    'laccorrY', 'alllacnames', 'corres', 'record', 'cutpts','fpsL')
            else
                oldnolac = load(fullfile(savedir,'nolacdata.mat'));
                oldlac = load(fullfile(savedir,'lacdata.mat'));

                nolacX = [oldnolac.nolacX, nolacXtemp];
                nolacY = [oldnolac.nolacY, nolacYtemp];
                nolacANAL_RMS = [oldnolac.nolacANAL_RMS, nolacANAL_RMStemp];
                nolaccorrX = [oldnolac.nolaccorrX, nolaccorrXtemp];
                nolaccorrY = [oldnolac.nolaccorrY, nolaccorrYtemp];
                fpsNL = [oldnolac.fpsNL fpsNLtemp];

                lacX = [oldlac.lacX, lacXtemp];
                lacY = [oldlac.lacY, lacYtemp];
                lacANAL_RMS = [oldlac.lacANAL_RMS, lacANAL_RMStemp];
                record = [oldlac.record, recordtemp];
                cutpts = [oldlac.cutpts, cutptstemp];
                laccorrX = [oldlac.laccorrX, laccorrXtemp];
                laccorrY = [oldlac.laccorrY, laccorrYtemp];
                fpsL = [oldlac.fpsL, fpsLtemp];

                save(fullfile(savedir,'nolacdata'),'nolacANAL_RMS','nolacX', 'nolacY',...
                    'nolaccorrX', 'nolaccorrY', 'nolacnames','fpsNL');
                save(fullfile(savedir,'lacdata'),'lacANAL_RMS', 'lacX', 'lacY', 'laccorrX',...
                    'laccorrY', 'alllacnames', 'corres', 'record', 'cutpts','fpsL')

                clear oldnolac oldlac
            end

            disp(strcat('Lac Set ', int2str(i),'/',int2str(length(lacnames)),' Finished'))
            close all
            if i<length(lacnames) %If this is the last time through this loop, don't delete these things so can call them below
                clear lacX lacY laccorrX laccorrY nolacX  nolacY nolaccorrX nolaccorrY
                clear nolacANAL_RMS lacANAL_RMS record cutpts fpsL fpsNL
            end
            clear fpsNLtemp fpsLtemp success lacr1 lacRMS nolacr1 pxlpath lacXtemp ...
                lacYtemp laccorrXtemp laccorrYtemp nolacXtemp nolacYtemp nolaccorrXtemp nolaccorrYtemp...
                nolacANAL_RMStemp lacANAL_RMStemp recordtemp cutptstemp

        end

        %Load any additional sets on same areas
        if ~(length(lacnames2)==1 && strcmpi(lacnames2{1},'')) %Second set on at least one area
            alllacnames = [lacnames lacnames2];
            %oldnolac = load(fullfile(savedir,'nolacdata.mat')); %The lac2 sets have the same nolac data as the first sets,
                %so just use previous analysis which isn't cleared after the
                %last set above
            %nolacANAL_RMS = oldnolac.nolacANAL_RMS;
            for i=1:length(lacnames2)
                if strcmpi(lacnames2{i},'')
                    lacXtemp = [];
                    lacYtemp = [];
                    laccorrXtemp=[];
                    laccorrYtemp=[];
                    lacANAL_RMStemp=[];
                    recordtemp = [];
                    cutptstemp = [];
                    fpsLtemp = [];
                else
                    [lacXtemp,lacYtemp,fpsLtemp,success,wind,clA]=loadrawV3(path, char(lacnames2(i)), strcat('Lac2Set',int2str(i)));
                    if ~success
                        return
                    end

                    %Drift correct and Gaussian filter the first-set lac data:
                    [laccorrXtemp,laccorrYtemp, lacr1]= meandriftcorrfunctionV3(lacXtemp,lacYtemp,fpsLtemp,fcB);
                    %**
                    %To drift correct with Martin's code instead:
                    %laccorrXtemp=BWdriftcorrect_ML(lacXtemp',fcB,fpsLtemp);
                    %laccorrYtemp=BWdriftcorrect_ML(lacYtemp',fcB,fpsLtemp);
                    %**
                    lacRMS = sqrt(gaussfiltV3(lacr1.^2, fcG, fpsLtemp));
                    %**
                    %To Gaussian filter with Martin's code instead:
                    %[lacRMS,notneeded,notneeded,notneeded]=gaussFilter_ML([laccorrXtemp' laccorrYtemp'],4,fpsLtemp);
                    %**

                    %Allow user to screen the data by watching movies:
                    pxlpath = fullfile(path,char(lacnames2(i)));
                    [lacANAL_RMStemp,recordtemp,cutptstemp] = screenbeadsV3(nolacANAL_RMS{i},...
                                    lacRMS,fpsLtemp, fpsNL{i},char(nolacnames(i)),char(lacnames2(i)),...
                                    pxlpath,tracesavedir, unloopedlength,wind,clA);  %Here we just assume the same fps for lac and nolac out of laziness ...

                end

                %Update the lacdata and nolacdata files--note that there won't ever
                %not be a file already created for these second sets, but don't
                %need to reload for the first lac2 set
                if i>1
                    oldlac = load(fullfile(savedir,'lacdata.mat'));
                    lacX = oldlac.lacX;
                    lacY = oldlac.lacY;
                    lacANAL_RMS = oldlac.lacANAL_RMS;
                    record = oldlac.record;
                    cutpts = oldlac.cutpts;
                    laccorrX = oldlac.laccorrX;
                    laccorrY = oldlac.laccorrY;
                    fpsL = oldlac.fpsL;
                end

                lacX{end+1} = lacXtemp; %Can't concatenate empty matrices onto cells, so need to do this in case any of the lacX, etc are empty
                lacY{end+1} = lacYtemp;
                lacANAL_RMS{end+1} = lacANAL_RMStemp;
                record{end+1} = recordtemp;
                cutpts{end+1} = cutptstemp;
                laccorrX{end+1} = laccorrXtemp;
                laccorrY{end+1} = laccorrYtemp;
                fpsL{end+1} = fpsLtemp;

                save(fullfile(savedir,'lacdata'),'lacANAL_RMS', 'lacX', 'lacY', 'laccorrX',...
                        'laccorrY', 'alllacnames', 'corres', 'record', 'cutpts','fpsL')

                disp(strcat('Lac 2 Set ', int2str(i),'/',int2str(length(lacnames2)),' Finished'))
                close all   

                if i<length(lacnames2) %Unless this is the last set, clear all variables
                    clear lacX lacY laccorrX laccorrY lacANAL_RMS record cutpts fpsL
                end
                clear fpsLtemp success lacr1 lacRMS pxlpath lacXtemp ...
                    lacYtemp laccorrXtemp laccorrYtemp lacANAL_RMStemp recordtemp cutptstemp

            end
        end

        %Create the corres2 and corres3 matrices.  Corres2
        %relates the beads numbers for multiple sets on the same area and is used
        %primarily for bdconcat; corres3 relates the bead numbers in lacdataconcat
        %to those in the saved nolac data file.

        [corres2, corres3] = createCorres2and3(lacnames,alllacnames,record,nolacANAL_RMS);

        %Save the lacdata.mat file a final time, this time with corres2 and corres3
        save(fullfile(savedir,'lacdata'),'lacANAL_RMS', 'lacX', 'lacY', 'laccorrX',...
                        'laccorrY', 'alllacnames', 'corres', 'record', 'cutpts','fpsL',...
                        'corres2','corres3')
        disp('Lac Data Saved')

        %We used to make a histogram of all the beads summed, but Steph never uses it.
        %Use probhistfunctionV2 if necessary.

        clear lacX lacY laccorrX lacorrY nolacX nolacY nolaccorrX nolaccorrY record cutpts
        clear tracesavedir

        %If there were second sets, concatenate data from different sets that are on the same beads:
        if ~(isequal(length(lacnames),length(alllacnames)))
            [newlacANAL_RMS,newlacnames]=bdconcatV3(lacANAL_RMS,corres2,alllacnames);
        else
            newlacANAL_RMS = lacANAL_RMS;
            newlacnames = alllacnames;
        end

        save(fullfile(savedir,'lacdataconcat'),'newlacANAL_RMS', 'newlacnames', 'corres2','corres3','fpsL')
        disp('Lacdataconcat saved')

        %Give the option of continuing with single-bead analysis:
        cont=input('Perform single-bead anaylsis now? (y/n)','s');
        if strcmpi(cont,'n')
            return
        else
            [GaussFit,type,pLoop,SEpLoop,rest] = singlebdanalysisV3(nolacANAL_RMS,nolacnames,...
                newlacANAL_RMS,corres3,newlacnames,fpsL{1},savedir);

            close all

             mkdir(fullfile(savedir,'Single bead analysis'))

            %Save GaussFIt
            save(fullfile(savedir,'Single bead analysis','GaussFit'),'GaussFit')
            %Save stats
            if ~strcmpi(type,'WT')
                pB = rest{1};
                SEpB = rest{2};
                pM = rest{3};
                SEpM = rest{4};
                save(fullfile(savedir,'Single bead analysis','stats'),'pLoop','SEpLoop','pB','SEpB','pM','SEpM')
            else %means wt data i.e. 5 states
                pO23loop = rest{1};
                SEO23 = rest{2};
                pO12loop = rest{3};
                SEO12 = rest{4};
                pOtherloop = rest{5};
                SEotherloop = rest{6};
                save(fullfile(savedir,'Single bead analysis','stats'),'pLoop','SEpLoop','pO23loop',...
                    'SEO23','pO12loop','SEO12','pOtherloop','SEotherloop')
            end

            makeallreptracesV3('Full Analysis Figures',savedir,newlacANAL_RMS,newlacnames,nolacANAL_RMS,...
                nolacnames,corres3,GaussFit,fpsL{1},type)
        end
    %end

end