%Steph 5/3/09

%Instead of just looking at the mean and standard deviations of the pLoops
%for each data set, can we learn anything from looking at the
%distributions?  E.g. number of beads that never loop.  Based loosely on
%PlotIndividualBeads_SJedit and related.

function plot_pLoop_distrib()

num=input('How many data sets?');

for a=1:num
    dir=uigetdir('/Volumes/dumbo-3/stephj/TPM data analysis','Choose directory with GaussFit.mat');

    file=load(fullfile(dir,'GaussFit.mat'));
    GaussFit=file.GaussFit;
    stats=load(fullfile(dir,'stats.mat'));
    
    pLoop=stats.pLoop;
    lower=pLoop-stats.SEpLoop;
    upper=pLoop+stats.SEpLoop;

    GaussFit=GaussFit(logical([GaussFit.Approved]));%removes any bead(s) whose Approved field was 0

    pLoops=zeros(length(GaussFit),1);

    for j=1:length(GaussFit)%for all the beads in GaussFit
        strcat('Bead ',int2str(j),'/',int2str(length(GaussFit))) %prints to tell you how many you've looked at so far
        clear thisbead
        if sum(strcmp(GaussFit(j).Identity,'B'))
            thisbead.pB=GaussFit(j).p(strcmp(GaussFit(j).Identity,'B'));%strcmp=string compare; 
        else
            thisbead.pB=0;
        end

        if sum(strcmp(GaussFit(j).Identity,'M'))
            thisbead.pM=GaussFit(j).p(strcmp(GaussFit(j).Identity,'M'));
        else
            thisbead.pM=0;
        end

        if sum(strcmp(GaussFit(j).Identity,'U'))
            thisbead.pU=GaussFit(j).p(strcmp(GaussFit(j).Identity,'U'));
        else
            thisbead.pU=0;
        end

        if sum(strcmp(GaussFit(j).Identity,'L'))
            thisbead.pL=GaussFit(j).p(strcmp(GaussFit(j).Identity,'L'));
        else
            thisbead.pL=0;
        end

        if (thisbead.pL==0)&&(thisbead.pM==0)&&(thisbead.pB==0)
            pLoops(j)=0;
        elseif thisbead.pL==0
            pLoops(j)=thisbead.pB+thisbead.pM;
        else
            pLoops(j)=thisbead.pL;
        end

    end

    %[n, xout]=hist(pLoops,length(GaussFit)); %The second argument here is pretty random ... what's a good bin number or size?
    [n, xout]=hist(pLoops,50);
    figure
    hold on
    bar(xout,n./sum(n))
    plot([pLoop pLoop],[min(n./sum(n)) max(n./sum(n))],'k',[lower lower],[min(n./sum(n)) max(n./sum(n))],'k--',[upper upper],[min(n./sum(n)) max(n./sum(n))],'k--','Linewidth',2)
    hold off
    xlim([0 1])
    ylim([0 max(n./sum(n))])
    pause
end
