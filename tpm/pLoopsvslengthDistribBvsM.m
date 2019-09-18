%[fighandle, Gaussmeans, GaussSEs, totmeans, totSEs, alldistribs,...
%    totmeansB,totSEsB,alldistribsB,totmeansM,totSEsM,alldistribsM] = pLoopsvslengthDistribBvsM(GaussFit,figtitle,savepath,varargin)
%
%Variant of pLoopsvslengthDistrib but deals with the two looped states.
%
%Right now doesn't do the plots or the Gaussian fit (commented out towards
%the end) to save time.
%
%Steph 3/11

function [fighandle, Gaussmeans, GaussSEs, totmeans, totSEs, alldistribs,...
    totmeansB,totSEsB,alldistribsB,totmeansM,totSEsM,alldistribsM] = pLoopsvslengthDistribBvsM(GaussFit,figtitle,savepath,varargin)

if length(GaussFit)==1
    if isfield(GaussFit,'GaussFit')
        GaussFit = GaussFit.GaussFit;
    else
        GaussFit = GaussFit.ThreshFit;
    end
end

histrange = [0:0.05:1];
binvect = [75:2:175];%FOR SHORTLOOPS OR HGSEQS ONLY!
disp('Binvect for Shortloops or HGseqs only!')

%%%%%This is the same as pLoopvslength

fps = 30 %Frame rate (could make this an input ...)

GaussFit=GaussFit(logical([GaussFit.Approved])); %Keep only beads that were approved by the user

data = zeros(length(GaussFit),2); %data will be a length by pLoop matrix, with each row corresponding to one bead
dataB = zeros(length(GaussFit),2);
dataM = zeros(length(GaussFit),2);

[data(:,2), totmean, totSE,dataB(:,2),dataM(:,2)]=extractIndivBdpLoops(GaussFit,'all'); %AllpLoops will have each bead's looping probability

for b = 1:length(GaussFit)
    if isfield(GaussFit,'n')
        %Each GaussFit element has two fields, n and all, which are related to the
        %histogram of the trace.  But all is normalized by the total number of
        %counts.  So to know how many non-zero data points were in the trace, we
        %want
        tempn = GaussFit(b).n;
        data(b,1) = sum(tempn(2:end));
        clear tempn
    else
        [tempn, xout] = hist(GaussFit(b).trace, binvect); %Want to find the number of nonzero points in this trace to estimate its length
        data(b,1) = sum(tempn(2:end));
        dataB(b,1) = sum(tempn(2:end));
        dataM(b,1) = sum(tempn(2:end));
        clear tempn 
    end
end

%%%%%%%%%%

%Now make histograms of all beads that last >= 1000 sec, >= 2000 sec, etc

%First sort the data matrix in order of increasing trace length:
data2 = sortrows(data,1); %Sort sorts the columns separately!
data2B = sortrows(dataB,1);
data2M = sortrows(dataM,1);

%Decide the largest bin according to the nearest maximal thousand seconds
steps = floor(data2(end,1)./fps/1000);
xticklabels=cell(length(steps*1000:-1000:1000),1);

hists = zeros(length(histrange),steps);
histsB = zeros(length(histrange),steps);
histsM = zeros(length(histrange),steps);

Gaussmeans = zeros(steps,1); %Only fit Gaussians to the whole distribution, since I don't really use this anyway
GaussSEs = zeros(steps,1);

totmeans = zeros(steps,1);
totSEs = zeros(steps,1);
alldistribs = cell(steps,1);
totmeansB = zeros(steps,1);
totSEsB = zeros(steps,1);
alldistribsB = cell(steps,1);
totmeansM = zeros(steps,1);
totSEsM = zeros(steps,1);
alldistribsM = cell(steps,1);

%Right now only making a figure with beads >= 5000 sec and shorter
% figure
% hindex=1;
% hands=cell(5,1);
% for h=1:5
%     subhands{1}=subplot(5,3,hindex);
%     subhands{2}=subplot(5,3,hindex+1);
%     subhands{3}=subplot(5,3,hindex+2);
%     hindex=hindex+3;
%     hands{h}=subhands;
%     clear subhands
% end
% 
% hindex = 5;

%Those figures are pretty ugly.  Now making separate figures for each 1000
%seconds: histogram, means, SEs

for i = steps*1000:-1000:1000
    mask = data2(:,1)./fps>=i; %This is a vector of logicals the same length as data2 which is 1 iff the length of that trace is greater than i
    
    thisdata = data2(mask,:); %This picks out the rows of data2 (i.e., traces) which have length > i seconds
    thisdataB = data2B(mask,:);
    thisdataM = data2M(mask,:);
    
    hists(:,i/1000) = hist(thisdata(:,2),histrange); %Each *column* of hists is the hist output for a time step
    counts = sum(hists(:,i/1000));
    hists(:,i/1000) = hists(:,i/1000)./(sum(hists(:,i/1000)));
    
    histsB(:,i/1000) = hist(thisdataB(:,2),histrange);
    %countsB = sum(histsB(:,i/1000));
    histsB(:,i/1000) = histsB(:,i/1000)./(sum(histsB(:,i/1000)));
    
    histsM(:,i/1000) = hist(thisdataM(:,2),histrange);
    %countsM = sum(histsM(:,i/1000));
    histsM(:,i/1000) = histsM(:,i/1000)./(sum(histsM(:,i/1000)));
    
    xticklabels{i/1000} = int2str(i);
    
    totmeans(i/1000) = mean(thisdata(:,2));
    totSEs(i/1000) = std(thisdata(:,2))/sqrt(size(thisdata,1)-1);
    totmeansB(i/1000) = mean(thisdataB(:,2));
    totSEsB(i/1000) = std(thisdataB(:,2))/sqrt(size(thisdataB,1)-1);
    totmeansM(i/1000) = mean(thisdataM(:,2));
    totSEsM(i/1000) = std(thisdataM(:,2))/sqrt(size(thisdataM,1)-1);
%     if i/1000<= 5
%         thishands = hands{hindex};
%         hindex=hindex-1;
%         [Gaussmeans(i/1000), GaussSEs(i/1000)] = fitpLoopdistrib(hists(:,i/1000)', histrange, counts,strcat(figtitle,'>=',int2str(i),'sec'),mean(thisdata(:,2)),thishands{1});
%         HowManyBeadsV2(thisdata(:,2),thishands{2},thishands{3});
%         drawnow
%         disp(strcat('Num bds = ',int2str(size(thisdata,1)),'for >=',int2str(i),' sec'))
%         pause(1)
%     else
%         [Gaussmeans(i/1000), GaussSEs(i/1000)] = fitpLoopdistrib(hists(:,i/1000)', histrange, counts,figtitle,mean(thisdata(:,2)));
%         HowManyBeadsV2(thisdata(:,2));
%         disp(strcat('Num bds = ',int2str(size(thisdata,1)),'for >=',int2str(i),' sec'))
%         pause(2)
%         close
%         close
%         close
%     end

    %Figures and Gaussian fitting commented out for now:
    %figure
    %h1=subplot(1,3,1);
    %h2=subplot(1,3,2);
    %h3=subplot(1,3,3);
    
    %[Gaussmeans(i/1000), GaussSEs(i/1000)] = fitpLoopdistrib(hists(:,i/1000)', histrange, counts,strcat(figtitle,'>=',int2str(i),'sec'),mean(thisdata(:,2)),h1);
    Gaussmeans(i/1000) = 0;
    GaussSEs(i/1000) = 0;
    %HowManyBeadsV2(thisdata(:,2),h2,h3);
    %set(gcf,'Position',[131, 227, 1056, 249])%This sets the figure to a nice size
    %set(gcf,'PaperPositionMode','auto') %This makes sure that when you save the figure (by using "print -depsc filename") it's the same nice size
    %drawnow
    %print('-depsc',strcat('/Users/Steph/Desktop/Figs/',figtitle,int2str(i),'sec'))
    %print('-djpeg',strcat('/Users/Steph/Desktop/Figs/',figtitle,int2str(i),'sec'))
    %print('-djpeg',strcat(savepath,'/',figtitle,int2str(i),'sec'))
    %disp(strcat('Num bds = ',int2str(size(thisdata,1)),'for >=',int2str(i),' sec'))
    %pause(1)
    %close
    %clear h1 h2 h3
    
    alldistribs{i/1000} = thisdata;
    alldistribsB{i/1000} = thisdataB;
    alldistribsM{i/1000} = thisdataM;
end

%As above, commented out:
fighandle = 1;
% if nargin==3
%     fighandle = figure;
% else
%     subplot(varargin{1})
%     fighandle = varargin{1};
% end

% PlotHandleb = bar3(histrange,hists,'b');
% ylim([0 1]);
% %ylabel('Looping Probability')
% %xlabel('Minimum Trace Length (seconds)')
% set(gca,'XTickLabel',xticklabels)
% %zlabel('Frequency')
% zlim([0 0.7])
% title(figtitle)
% 
% StandardFigureBar(PlotHandleb,gca)
% pause(4)
% close