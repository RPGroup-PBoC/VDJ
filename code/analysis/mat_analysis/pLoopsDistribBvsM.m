%function pLoopsDistribBvsM()
%
%Makes histograms like pLoopsvslengthDistrib (without the bootstrapped
%panels) but for the two looped states separately.  Inputs are the
%distributions for a particular cutoff time (eg, >=3000 sec) for the total
%pLoop, pB, and pM; the mean and mean minus nonloopers for all three
%distributions; the name that is the first part of the figure titles; and
%the path to a directory of where to save the figures.
%
%Steph 5/2011

function pLoopsDistribBvsM(totdistrib,totmean, totmeanMinusNL, ...
    Bdistrib,Bmean, BmeanMinusNL, Mdistrib, Mmean, MmeanMinusNL, name,savepath)

histrange = [0:0.05:1];

totdistribhisto = hist(totdistrib,histrange);
%totdistribhisto = totdistribhisto./sum(totdistribhisto); %Normalize by the
%total number of counts -- right now not going to normalize anything
totbeads = sum(totdistribhisto); %Total number of beads
NLs = totdistribhisto(1); %This is the number of beads that never show any looping
NLshisto = zeros(1,size(histrange,2));
NLshisto(1) = NLs;

Bdistribhisto = hist(Bdistrib,histrange);
%Bdistribhisto = Bdistribhisto./sum(Bdistribhisto); %Normalize by the total number of counts
%Bbeads = sum(Bdistribhisto); %Total number of beads
Bbeads = sum(Bdistribhisto)-Bdistribhisto(1); %Total number of beads that have nonzero B probabilities

Mdistribhisto = hist(Mdistrib,histrange);
%Mdistribhisto = Mdistribhisto./sum(Mdistribhisto); %Normalize by the total number of counts
%Mbeads = sum(Mdistribhisto); %Total number of beads
Mbeads = sum(Mdistribhisto)-Mdistribhisto(1); %Total number of beads that have nonzero M probabilities

figure
subplot(1,3,1);
    
    PlotHandle = bar(histrange,totdistribhisto,'k');
    hold on
    PlotHandleb=plot([totmean totmean],[0 50],'--m');
    PlotHandleb(end+1)=plot([totmeanMinusNL totmeanMinusNL],[0 50],'-m');
    title(name)
    xlim([0 1])
    set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1])
    %set(gca,'YTick',[0 0.2 0.4 0.6])
    %ylim([0 0.7])
    ylim([0 30])
    xlabel('Looping Probability')
    ylabel('Counts')
    
    legend(strcat('N_{tot}=',int2str(totbeads)))

    StandardFigureBar(PlotHandle,gca)
    StandardFigure(PlotHandleb,gca)
    
subplot(1,3,2);
    
    PlotHandleB = bar(histrange,Bdistribhisto,'b');
    hold on
    %Plot the total nonloopers in a different color
    PlotHandleB(end+1) = bar(histrange,NLshisto,'k');
    PlotHandlebB=plot([Bmean Bmean],[0 50],'--k');
    PlotHandlebB(end+1)=plot([BmeanMinusNL BmeanMinusNL],[0 50],'-k');
    title(strcat(name,', B'))
    xlim([0 1])
    set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1])
    %set(gca,'YTick',[0 0.2 0.4 0.6])
    %ylim([0 0.7])
    ylim([0 30])
    xlabel('Looping Probability')
    ylabel('Counts')
    
    legend(strcat('N_{loop}=',int2str(Bbeads)))

    StandardFigureBar(PlotHandleB,gca)
    StandardFigure(PlotHandlebB,gca)
    
subplot(1,3,3);

    PlotHandleM = bar(histrange,Mdistribhisto,'r');
    hold on
    %Plot the total nonloopers in a different color
    PlotHandleM(end+1) = bar(histrange,NLshisto,'k');
    PlotHandlebM=plot([Mmean Mmean],[0 50],'--k');
    PlotHandlebM(end+1)=plot([MmeanMinusNL MmeanMinusNL],[0 50],'-k');
    title(strcat(name,', M'))
    xlim([0 1])
    set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1])
    %set(gca,'YTick',[0 0.2 0.4 0.6])
    %ylim([0 0.7])
    ylim([0 30])
    xlabel('Looping Probability')
    ylabel('Counts')
    
    legend(strcat('N_{loop}=',int2str(Mbeads)))

    StandardFigureBar(PlotHandleM,gca)
    StandardFigure(PlotHandlebM,gca)
    
    
set(gcf,'Position',[131, 227, 1056, 249])%This sets the figure to a nice size
set(gcf,'PaperPositionMode','auto') %This makes sure that when you save the figure (by using "print -depsc filename") it's the same nice size
drawnow
print('-depsc',fullfile(savepath,strcat(name,'NiceHistos')))
close

    
    
    
    
    
    


