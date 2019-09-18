%[GaussFit,type,pLoop,SEpLoop,varargout] = singlebdanalysisV3(nolacANAL_RMS,nolacnames,newlacANAL_RMS,corres3,newlacnames,fps,storedir,varargin)
%
%Allows user to fit Gaussians to each bead's trace in a lacdataconcat.mat
%file and calculates the looping probability.  Called by masterscriptV3.
%
%Stephanie Johnson 2/11

function [GaussFit,type,pLoop,SEpLoop,varargout] = singlebdanalysisV3(nolacANAL_RMS,nolacnames,newlacANAL_RMS,corres3,newlacnames,fps,storedir,varargin)

%For better guides to the eye for singlebdanalysis: calculate each bead's
        %individual unlooped length:
unloopedlengths = calc_unlooped_lengths(nolacANAL_RMS,corres3);

%Perform single-bead analysis.  This used to be done through
%singlebdanalysis which converted to an intermediate format, a holdover
%from Lin's code which we no longer use.  V3 eliminates that
%intermediate format and calls FitIndivBds from here directly

%First put all the beads into one cell array:
k =1; %bead counter
allbeads=[]; %This will contain all the bds RMS values, so that FitIndivBds can create a summed histogram and get starting fit parameters from that
for i = 1:length(newlacANAL_RMS)
    [beads frames] = size(newlacANAL_RMS{i});
    for j = 1:beads
        newr{k}=transpose(newlacANAL_RMS{i}(j,:));
        allbeads=[allbeads;newr{k}];
        allnames{k}=strcat(newlacnames(i),'Bead', int2str(k));
        k = k+1;
    end
end

if isempty(varargin)
    [GaussFit,type]=FitIndivBdsSJ(newr,allbeads,allnames,fps,unloopedlengths);
    %mkdir(fullfile(storedir,'Single bead analysis')) %Put this in a separate folder for consistency with old code
    %save(fullfile(storedir,'Single bead analysis','GaussFit'),'GaussFit')%this is now done in masterscriptV3
    %   main file
    %Call calc_ploop, which calculates pLoop and SEpLoop and saves them in a
    %separate file called stats--the saving is now done in masterscriptV3
    %   main code
    if ~strcmpi(type,'WT')
        [pLoop,SEpLoop,pB,SEpB,pM,SEpM]=calc_ploop(GaussFit,0);
        %save(fullfile(storedir,'Single bead analysis','stats'),'pLoop','SEpLoop','pB','SEpB','pM','SEpM')
        rest{1} = pB;
        rest{2} = SEpB;
        rest{3} = pM;
        rest{4} = SEpM;
        varargout{1} = rest;
    else %means wt data i.e. 5 states
        [pLoop,SEpLoop,pO23loop,SEO23,pO12loop,SEO12,pOtherloop,SEotherloop]=calc_ploop(GaussFit,1);
        %save(fullfile(storedir,'Single bead analysis','stats'),'pLoop','SEpLoop','pO23loop',...
        %    'SEO23','pO12loop','SEO12','pOtherloop','SEotherloop')
        rest{1} = pO23loop;
        rest{2} = SEO23;
        rest{3} = pO12loop;
        rest{4} = SEO12;
        rest{5} = pOtherloop;
        rest{6} = SEotherloop;
        varargout{1} = rest;
    end
else 
    [GaussFit,type]=FitIndivBds(newr,allbeads,allnames,fps,unloopedlengths,varargin{1});
    %save(fullfile(storedir,'Single bead analysis','newGaussFit'),'GaussFit')%So don't save over the old one; the user 
    %will have to change the file's name after it's been saved in order to
    %not save over it again
    if ~strcmpi(type,'WT')
        [pLoop,SEpLoop,pB,SEpB,pM,SEpM]=calc_ploop(GaussFit,0);
        %save(fullfile(storedir,'Single bead analysis','newstats'),'pLoop','SEpLoop','pB','SEpB','pM','SEpM')
        rest{1} = pB;
        rest{2} = SEpB;
        rest{3} = pM;
        rest{4} = SEpM;
        varargout{1} = rest;
    else %means wt data i.e. 5 states
        [pLoop,SEpLoop,pO23loop,SEO23,pO12loop,SEO12,pOtherloop,SEotherloop]=calc_ploop(GaussFit,1);
        %save(fullfile(storedir,'Single bead analysis','newstats'),'pLoop','SEpLoop','pO23loop',...
         %   'SEO23','pO12loop','SEO12','pOtherloop','SEotherloop')
         rest{1} = pO23loop;
        rest{2} = SEO23;
        rest{3} = pO12loop;
        rest{4} = SEO12;
        rest{5} = pOtherloop;
        rest{6} = SEotherloop;
        varargout{1} = rest;
    end
end