%Steph 9/10

%Called by AvgPLoopByGaussFit--fits Gaussians and plots ploop distributions

%Totmean is the result of just taking the mean of all the pLoops; mean is
%from the Gauss fit in this function.

function [mean, SE, varargout] = fitpLoopdistrib(distribhist, xout, counts,name,totmean,varargin)

    xRange = linspace(0,1);
    
    %instead of pLoopshist, range as inputs could have inputs be pLoops, range
    %[n,xout]=hist(pLoops,range);
    %distribhist = n./sum(n);
    opts = fitoptions('gauss1',...
        'Lower',[0 0 0],'Upper',[Inf 1 Inf],'Algorithm','Trust-Region');
    try 
        tempfit = fit(xout',distribhist','gauss1',opts);
    catch
        disp('Cant fit a gaussian.')
        tempfit.a1 = 0;
        tempfit.b1 = 0;
        tempfit.c1 = 0;
        keyboard
    end
    mean = tempfit.b1;
    SE = sqrt(tempfit.c1/2)/sqrt(counts-1); %Is this right???

    Gauss=tempfit.a1*exp(-((xRange-tempfit.b1)/tempfit.c1).^2);

    if nargin==6
    	subplot(varargin{1});
        varargout{1}=varargin{1};
    else
        varargout{1}=figure;
    end
    PlotHandle = bar(xout,distribhist,'b');
    hold on
    PlotHandleb=plot([mean mean],[0 0.7],'--r');
    PlotHandleb(end+1)=plot(xRange,Gauss,'-r');
    PlotHandleb(end+1)=plot([totmean totmean],[0 0.7],'--k');
    %title(name,'FontSize',12,'FontName','Lucida Sans')
    title(name)
    xlim([0 1])
    ylim([0 0.7])
    set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1])
    set(gca,'YTick',[0 0.2 0.4 0.6])
    xlabel('Looping Probability')
    ylabel('Frequency')

    StandardFigureBar(PlotHandle,gca)
    StandardFigure(PlotHandleb,gca)

end