%function calcpLoop_threshold()
%
%Some traces are not easily fit by Gaussians, so this code allows the user
%to threshold the histogrammed trace in order to calculate the looping
%probability.  NOTE: this code depends on the output of masterscriptV3 and 
%is not backwards compatible with previous versions of that code!
%
%No input necessary; it will ask the user to choose the folder to analyze
%(this is the folder that contains lacdataconcat, etc.)
%
%Updated 8/2011 to load a previous ThreshFit object so that edits can be
%made.  Pass anything (doesn't matter what) into varargin to load a previous ThreshFit.
%
%Steph 2/11

function calcpLoop_threshold(varargin)

loaddir = uigetdir('/Volumes/Phillips Lab Data/TPM data analysis/SJLacI','Choose folder of data to analyze');

%Set up the paramters that will control how the histogram is binned and
%displayed and the start parameters for the fits:
ok=0;
while ok==0
    datatype = input('Shortloops, PUC, smallbds, WTKO, other?','s');
    if strcmpi(datatype,'shortloops') || strcmpi(datatype,'short loops')
        %xx=[50:2:200];
        TraceRange=[75,175];
        autotopline = 150;
        ok=1;
    elseif strcmpi(datatype,'PUC')
        %xx=[100:2:260];
        TraceRange=[100,260];
        autotopline = 200;
        ok=1;
    elseif strcmpi(datatype,'smallbds')
        %xx=[100:2:260];
        TraceRange=[50,125];
        autotopline = 120;
        ok=1;
    elseif strcmpi(datatype,'WTKO')
        %xx=[80:2:240];
        TraceRange=[75,240];
        autotopline = 200;
        ok=1;
    elseif strcmpi(datatype,'other') %Should do some error-handling here ...
        %xx = input('Enter x range (= start:increment:end) :');
        TraceRange = input('Enter TraceRange (= [start end]) :');
        autotopline = TraceRange(end)-1;
        ok=1;
    else
        disp('Invalid input.')
    end
end
xxRange=linspace(TraceRange(1),TraceRange(2));
binvect = [TraceRange(1):2:TraceRange(2)];

if nargin == 0
    mkdir(fullfile(loaddir,'Thresholding'))
else
    if exist(fullfile(loaddir,'Thresholding','ThreshFit.mat'))
        oldThresh = load(fullfile(loaddir,'Thresholding','ThreshFit.mat'));
        oldThresh = oldThresh.ThreshFit;
    else
        disp('No ThreshFit to load?')
        return
    end
end
savedir = fullfile(loaddir,'Thresholding');

nolacdata = load(fullfile(loaddir,'nolacdata.mat')); %This is only needed for one line, below
load(fullfile(loaddir,'lacdataconcat.mat'));
load(fullfile(loaddir,'Single bead analysis','GaussFit.mat'));
%Backwards compatibility with old code
if ~exist('fpsL')
    needautofps = 1;
    disp('Using 30 fps.')
end

%I will plot a line corresponding to the unlooped length of each
%trace, so need to calculate the average unlooped RMS for each bead.  Will
%use the variance to estimate the threshold values:
[unloopedlengths,unloopedwidths] = calc_unlooped_lengths(nolacdata.nolacANAL_RMS,corres3);
clear nolacdata

numbds = length(GaussFit);

%First use the fit results in GaussFit to put initial thresholding
%lines on the plot:
b=1; %Indexes GaussFit and unloopedlengths
for s = 1:length(newlacnames)
    for j = 1:size(newlacANAL_RMS{s},1)
        ThreshFit(b).name = strcat(newlacnames{s},'_Bead',int2str(j));
        if nargin==0
            %First, an upper bound above the unlooped state: choose this to be
            %three times the standard dev. of the Gaussian
            topline = unloopedlengths(b)+3*unloopedwidths(b)/sqrt(2);
            %Second, a lower bound to exclude sticking: set this at 80 nm,
            %though there can be stuck states closer to 100 nm ...
            %Steph added 5/2011: for dealing with small beads, where the
            %sticking states are closer to 50 nm:
            if ~strcmpi(datatype,'smallbds')
                bottomline = 80;
            else
                bottomline = 55;
            end
            %If for some reason the unloopedlength is faulty, topline could be
            %below bottomline which leads to Inf's and NaN's in the calculated
            %p's for this bead:
            if topline<=bottomline
                topline = autotopline;
            end
            %Also a problem if the tether is abnormally long ...
            if topline>=TraceRange(2)
                topline = TraceRange(2)-10;
            end
            %Third, use the number of Gaussians from GaussFit to put additional
            %lines
            if GaussFit(b).Approved
                ThreshFit(b).Approved = 1;
                if isequal(GaussFit(b).p, [0 0 1]) || isequal(GaussFit(b).p,[1 0]) %User set this to all looped or all unlooped--nothing to threshold
                    midline = [];
                    ThreshFit(b).Identity = GaussFit(b).Identity;
                elseif length(GaussFit(b).p)==3
                    %Put lines at the minima between the Gaussians
                    temprange = sort([GaussFit(b).fp.b1,GaussFit(b).fp.b2,GaussFit(b).fp.b3]);
                    minRange1 = linspace(temprange(1),temprange(2));
                    minRange2 = linspace(temprange(2),temprange(3));
                    summedGauss1 = GaussFit(b).fp.a1*exp(-((minRange1-GaussFit(b).fp.b1)./GaussFit(b).fp.c1).^2)+...
                        GaussFit(b).fp.a2*exp(-((minRange1-GaussFit(b).fp.b2)./GaussFit(b).fp.c2).^2)+...
                        GaussFit(b).fp.a3*exp(-((minRange1-GaussFit(b).fp.b3)./GaussFit(b).fp.c3).^2);
                    summedGauss2 = GaussFit(b).fp.a1*exp(-((minRange2-GaussFit(b).fp.b1)./GaussFit(b).fp.c1).^2)+...
                        GaussFit(b).fp.a2*exp(-((minRange2-GaussFit(b).fp.b2)./GaussFit(b).fp.c2).^2)+...
                        GaussFit(b).fp.a3*exp(-((minRange2-GaussFit(b).fp.b3)./GaussFit(b).fp.c3).^2);
                    [notneeded,ind1] = min(summedGauss1); 
                    midline(1) = minRange1(ind1);
                    [notneeded,ind2] = min(summedGauss2); 
                    midline(2) = minRange2(ind2);
                    clear summedGauss1 summedGauss2 temprange minRange1 minRange2
                    ThreshFit(b).Identity = {'B','M','U'};
                else
                    %Put a line at the minimum between the two Gaussians
                    temprange = sort([GaussFit(b).fp.b1,GaussFit(b).fp.b2]);
                    minRange = linspace(temprange(1),temprange(2));
                    summedGauss = GaussFit(b).fp.a1*exp(-((minRange-GaussFit(b).fp.b1)./GaussFit(b).fp.c1).^2)+...
                        GaussFit(b).fp.a2*exp(-((minRange-GaussFit(b).fp.b2)./GaussFit(b).fp.c2).^2);
                    [notneeded,ind] = min(summedGauss); 
                    midline = minRange(ind);
                    clear summedGauss temprange minRange
                    %Guess that if this is the middle looped state the midline will
                    %be closer to topline than bottom line
                    if topline-midline < (topline-bottomline)/2
                        ThreshFit(b).Identity = {'M','U'};
                    else
                        ThreshFit(b).Identity = {'B','U'};
                    end
                end
            else
                ThreshFit(b).Approved = 0;
                midline=[];
            end
        else %use values in oldThresh
            topline = oldThresh(b).topline;
            bottomline = oldThresh(b).bottomline;
            midline = oldThresh(b).midline;
            ThreshFit(b).Approved = oldThresh(b).Approved;
            ThreshFit(b).Identity = oldThresh(b).Identity;
        end
        ThreshFit(b).trace = newlacANAL_RMS{s}(j,:);
        ThreshFit(b).topline = topline;
        ThreshFit(b).bottomline = bottomline;
        ThreshFit(b).midline = midline;
        clear topline bottomline midline
        b=b+1;
    end
end

b=1; %index the current bead--this is for GaussFit/Threshfit and unloopedlengths, which don't know anything about sets
s = 1; %indexes the set--this is for newlacnames and newlacANAL_RMS
j = 1; %indexes the bead in the current set
LineSelected = 0;

%Interactive section
while b<=numbds
    
    %Backwards compatibility with old code
    if exist('needautofps')
        fpsL{s} = 30;
    end
    
    name=strrep(ThreshFit(b).name,'_','\_');    
    %Plot the trace and its histogram
        figure('Position',[150 300 950 500])
            subplot('Position', [.08,.15,.62,.75])
            title(strcat(name,',',int2str(b),'/',int2str(numbds)),'Fontsize',16);
            hold on
            ylim(TraceRange)
            xlabel('Time (sec)','Fontsize',16);
            ylabel('<R> (nm)','Fontsize',16)
            set(gca,'FontSize',14)

            plot([1/fpsL{s}:1/fpsL{s}:size(newlacANAL_RMS{s},2)/fpsL{s}],newlacANAL_RMS{s}(j,:),'-b')
            plot([0 size(newlacANAL_RMS{s},2)/fpsL{s}],[unloopedlengths(b) unloopedlengths(b)],'k--')%Plot the unlooped length
            xlim([0 size(newlacANAL_RMS{s},2)/fpsL{s}])
            %Plot the thresholding lines, only if this bead is kept
            if ThreshFit(b).Approved
                if LineSelected == 1
                    plot([0 size(newlacANAL_RMS{s},2)/fpsL{s}],[ThreshFit(b).topline ThreshFit(b).topline], 'r--','Linewidth',2)
                else
                    plot([0 size(newlacANAL_RMS{s},2)/fpsL{s}],[ThreshFit(b).topline ThreshFit(b).topline], 'r--')
                end
                if LineSelected == 2 && isempty(ThreshFit(b).midline) || ...
                        (LineSelected == 3 && ~isempty(ThreshFit(b).midline) && length(ThreshFit(b).midline) == 1) || ...
                        (LineSelected == 4 && ~isempty(ThreshFit(b).midline) && length(ThreshFit(b).midline) == 2)
                    plot([0 size(newlacANAL_RMS{s},2)/fpsL{s}],[ThreshFit(b).bottomline ThreshFit(b).bottomline], 'r--','Linewidth',2)
                else
                    plot([0 size(newlacANAL_RMS{s},2)/fpsL{s}],[ThreshFit(b).bottomline ThreshFit(b).bottomline], 'r--')
                end
                if ~isempty(ThreshFit(b).midline)
                    for l = 1:length(ThreshFit(b).midline) 
                        if LineSelected == 2 && length(ThreshFit(b).midline)==1
                            plot([0 size(newlacANAL_RMS{s},2)/fpsL{s}],[ThreshFit(b).midline(l) ThreshFit(b).midline(l)], 'r--','Linewidth',2)
                        elseif (LineSelected == 2 && l==2) || (LineSelected == 3 && l==1 && length(ThreshFit(b).midline)==2)
                            plot([0 size(newlacANAL_RMS{s},2)/fpsL{s}],[ThreshFit(b).midline(l) ThreshFit(b).midline(l)], 'r--','Linewidth',2)
                        else
                            plot([0 size(newlacANAL_RMS{s},2)/fpsL{s}],[ThreshFit(b).midline(l) ThreshFit(b).midline(l)], 'r--')
                        end
                    end
                end
            end

            hold off

        subplot('Position', [.71,.15,.25,.75])
            [n, xout] = hist(newlacANAL_RMS{s}(j,:), binvect);
            prob = n./(sum(n(2:end))); %normalize by the total number of counts
            prob(1) = 0; %Because screenbeads sets 'bad' data to zero, the first bin will contain all these points
            hold on
            if ThreshFit(b).Approved
                barh(xout, prob)
                %plot([0 max(prob)+0.01],[unloopedlengths(b) unloopedlengths(b)],'k--')
                if LineSelected == 1
                    plot([0 max(prob)+0.01],[ThreshFit(b).topline ThreshFit(b).topline], 'r--','Linewidth',2)
                else
                    plot([0 max(prob)+0.01],[ThreshFit(b).topline ThreshFit(b).topline], 'r--')
                end
                if LineSelected == 2 && isempty(ThreshFit(b).midline)
                    plot([0 max(prob)+0.01],[ThreshFit(b).bottomline ThreshFit(b).bottomline], 'r--','Linewidth',2)
                else
                    plot([0 max(prob)+0.01],[ThreshFit(b).bottomline ThreshFit(b).bottomline], 'r--')
                end
                if ~isempty(ThreshFit(b).midline)
                    for l = 1:length(ThreshFit(b).midline)
                        if LineSelected ~=0 && LineSelected == l-1
                            plot([0 max(prob)+0.01],[ThreshFit(b).midline(l) ThreshFit(b).midline(l)], 'r--','Linewidth',2)
                        else
                            plot([0 max(prob)+0.01],[ThreshFit(b).midline(l) ThreshFit(b).midline(l)], 'r--')
                        end
                    end
                end
                %This is a tacky way to tell the user what the state labels
                %are, but it's the best I've come up with so far ...
                if length(ThreshFit(b).Identity) == 3
                    title(strcat('(',ThreshFit(b).Identity{1},',',ThreshFit(b).Identity{2},',',ThreshFit(b).Identity{3},')'))
                else
                    title(strcat('(',ThreshFit(b).Identity{1},',',ThreshFit(b).Identity{2},')'))
                end
            else
                barh(xout, prob,'r')
            end
            ylim(TraceRange)
            set(gca,'YTick',[]);
            set(gca,'XTick',linspace(0,max(prob)+0.01,3));
            xlim([0 max(prob)+0.01])
            xlabel('Probability','Fontsize',16)
            set(gca,'FontSize',14)


    %User-interactive section--this uses similar commands to the
    %user-interaction section in FitIndivBds.  Cycle through while loop
    %till user gives an input.
    cc=1;
    while cc~=13
        ct=waitforbuttonpress;
        cc=get(gcf,'currentcharacter');

        if ct==1
            %Go forward to the next bead
            if cc=='.' 
                %First calculate this bead's looping probability
                if isempty(ThreshFit(b).midline) || ~ThreshFit(b).Approved %All looped or unlooped, or bead discarded
                    ThreshFit(b).p = GaussFit(b).p;
                else
                    %Compute the total counts between topline and bottom line--this is
                    %the area under the whole histogram.  First find the bins
                    %corresponding to the top and bottom lines:
                    temp = find(xout>ThreshFit(b).topline);
                    if isempty(temp)
                        topbin = length(xout);
                    else
                        topbin = temp(1);
                    end
                    clear temp
                    temp = find(xout>ThreshFit(b).bottomline);
                    bottombin = temp(1)-1;
                    clear temp
                    totarea = sum(prob(bottombin:topbin));
                    if length(ThreshFit(b).midline)==1 %only two states
                        %Find the bin in the bead's histogram that corresponds to the
                        %threshold:
                        temp = find(xout>ThreshFit(b).midline(1));
                        midbin = temp(1);
                        clear temp
                        ThreshFit(b).p = [sum(prob(bottombin:midbin-1))/totarea, sum(prob(midbin:topbin))/totarea];

                    else %two looped states
                        temp1 = find(xout>ThreshFit(b).midline(1));
                        midbin1 = temp1(1);
                        temp2 = find(xout>ThreshFit(b).midline(2));
                        midbin2 = temp2(1);
                        clear temp1 temp2
                        ThreshFit(b).p = [sum(prob(bottombin:midbin1-1))/totarea, ...
                            sum(prob(midbin1:midbin2-1))/totarea, sum(prob(midbin2:topbin))/totarea];
                     end
                end

                disp(ThreshFit(b).p)
                disp(ThreshFit(b).Identity)
                pause
                
                b = b+1;
                if j==size(newlacANAL_RMS{s},1) %means going to the next set
                    j=1;
                    s=s+1;
                else
                    j = j+1;
                end
                LineSelected = 0;
                cc=13;
            %Go back one bead
            elseif cc==',' 
                if b>1
                    b=b-1;
                    if s>1 && j==1 %means going back one set
                        s = s-1;
                        j=size(newlacANAL_RMS{s},1);
                    else
                        j=j-1;
                    end
                end
                LineSelected = 0;
                cc=13;
            %Select a line.  There can be between 2 and 4 lines
            elseif cc == '1'
                LineSelected = 1;
                cc=13;
            elseif cc == '2'
                LineSelected = 2;
                cc=13;
            elseif cc == '3' && ~isempty(ThreshFit(b).midline)
                LineSelected = 3;
                cc=13;
            elseif cc == '4' && ~isempty(ThreshFit(b).midline) && length(ThreshFit(b).midline) == 2
                LineSelected = 4;
                cc=13;                
            %Move a line.  I've stupidly numbered the lines top to bottom,
            %   but midline indexes them shortest to longest RMS ...
            elseif cc == 'm' && LineSelected ~= 0
                [x,y] = ginput(1);
                if LineSelected == 1 %Move the topline
                    ThreshFit(b).topline = y;
                elseif (LineSelected == 2 && isempty(ThreshFit(b).midline)) ||...
                        (LineSelected == 3 && ~isempty(ThreshFit(b).midline) && length(ThreshFit(b).midline) == 1) || ...
                        (LineSelected == 4 && ~isempty(ThreshFit(b).midline) && length(ThreshFit(b).midline) == 2)%move the bottomline
                    ThreshFit(b).bottomline = y;
                elseif LineSelected == 2 %%move the top of one or two midlines
                    if length(ThreshFit(b).midline) == 1
                        ThreshFit(b).midline(1) = y;
                    else
                        ThreshFit(b).midline(2) = y;
                    end
                elseif LineSelected == 3 %move the bottom of two midlines
                    ThreshFit(b).midline(1) = y;
                end
                cc=13;
            %Add a line--can only choose a new line between topline and
            %   bottomline, and can only add until there are 4 lines total.
            elseif cc == 'a' && (isempty(ThreshFit(b).midline) || length(ThreshFit(b).midline) == 1) %Can only add a line if there aren't already 4
                [x,y] = ginput(1);
                if y<ThreshFit(b).topline && y>ThreshFit(b).bottomline
                    if isempty(ThreshFit(b).midline)
                        ThreshFit(b).midline = y;
                        %Need to adjust the identity matrix.  As before,
                        %Guess that if this is the middle looped state the midline will
                        %be closer to topline than bottom line
                        if ThreshFit(b).topline-ThreshFit(b).midline < (ThreshFit(b).topline-ThreshFit(b).bottomline)/2
                            ThreshFit(b).Identity = {'M','U'};
                        else
                            ThreshFit(b).Identity = {'B','U'};
                        end
                    else
                        ThreshFit(b).Identity = {'B','M','U'};
                        if y<ThreshFit(b).midline(1)
                            ThreshFit(b).midline = [y ThreshFit(b).midline];
                        else
                            ThreshFit(b).midline = [ThreshFit(b).midline y];
                        end
                    end
                end
                LineSelected = 0;
                cc=13;
            %Remove a line--can't remove topline or bottomline (can move
            %   them though)
            elseif cc == 'r' && LineSelected ~= 0 && LineSelected ~= 1
                if LineSelected == 2 && ~isempty(ThreshFit(b).midline)
                    if length(ThreshFit(b).midline) == 2
                        ThreshFit(b).midline = ThreshFit(b).midline(1);
                        %Again adjust identity matrix
                        if ThreshFit(b).topline-ThreshFit(b).midline < (ThreshFit(b).topline-ThreshFit(b).bottomline)/2
                            ThreshFit(b).Identity = {'M','U'};
                        else
                            ThreshFit(b).Identity = {'B','U'};
                        end
                    else 
                        ThreshFit(b).midline = [];
                        %What to do about adjusting the identity matrix.
                        %This line executes if the user has removed all
                        %middle lines, leaving the bead totally looped or
                        %totally unlooped.  I don't think it'll be the case
                        %that the user will make everything looped or
                        %unlooped if they didn't already do that in the
                        %GaussFit ...
                    end
                elseif LineSelected == 3 && ~isempty(ThreshFit(b).midline) && length(ThreshFit(b).midline) == 2
                    ThreshFit(b).midline = ThreshFit(b).midline(2);
                    %Again adjust identity matrix
                    if ThreshFit(b).topline-ThreshFit(b).midline < (ThreshFit(b).topline-ThreshFit(b).bottomline)/2
                        ThreshFit(b).Identity = {'M','U'};
                    else
                        ThreshFit(b).Identity = {'B','U'};
                    end
                end
                LineSelected = 0;
                cc=13;
            %Discard or un-discard a bead
            elseif cc == 'd'
                ThreshFit(b).Approved = ~ThreshFit(b).Approved;
                cc = 13;
            %Rename the states
            elseif cc == 'n'
                ok = 0;
                while ok == 0
                    temp = input('Enter state labels from bottom to top: ','s');
                    if length(temp)==2 && length(ThreshFit(b).midline) == 1
                        ThreshFit(b).Identity = {temp(1),temp(2)};
                        ok = 1;
                    elseif length(temp) == 3 && length(ThreshFit(b).midline) == 2
                        ThreshFit(b).Identity = {temp(1),temp(2),temp(3)};
                        ok =1;
                    end
                end
            end

        end


    end
    
    close
end

if ~exist(fullfile(savedir,'ThreshFit'))
  save(fullfile(savedir,'ThreshFit'),'ThreshFit')
else
  save(fullfile(savedir,'ThreshFit'),'ThreshFit2')
end
%save('/Users/Steph/Desktop/ThreshFit.mat','ThreshFit')

%Calculate the total looping probability
allpB=[];
allpM=[];
allpL=[];

ThreshFit2=ThreshFit(logical([ThreshFit.Approved]));

%3/11: this has been wrong! (and is wrong in all Threshfit files for the
%length series).  allpB, etc need to have the same length as numbeads--if
%there's no bottom state, need to have allpB(b) = 0.  The means and SEs are
%all wrong!

for k=1:length(ThreshFit2)
    allpL(k) = 0;
    if sum(strcmp(ThreshFit2(k).Identity,'B'))
        allpB(end+1)=ThreshFit2(k).p(strcmp(ThreshFit2(k).Identity,'B'));
        allpL(k) = allpL(k)+allpB(end);
    end
    if sum(strcmp(ThreshFit2(k).Identity,'M'))
        allpM(end+1)=ThreshFit2(k).p(strcmp(ThreshFit2(k).Identity,'M'));
        allpL(k) = allpL(k)+allpM(end);
    end
    if sum(strcmp(ThreshFit2(k).Identity,'L'))
        allpL(k) = allpL(k)+ThreshFit2(k).p(strcmp(ThreshFit2(k).Identity,'L'));
    end
end

pLoop=mean(allpL)
SEpLoop=std(allpL)/sqrt(length(allpL)-1)
pB=mean(allpB);
SEpB=std(allpB)/sqrt(length(allpB)-1);
pM=mean(allpM);
SEpM=std(allpM)/sqrt(length(allpM)-1);

if ~exist(fullfile(savedir,'stats'))
  save(fullfile(savedir,'stats'),'pLoop','SEpLoop','pB','SEpB','pM','SEpM')
else
  save(fullfile(savedir,'stats2'),'pLoop','SEpLoop','pB','SEpB','pM','SEpM')
end



