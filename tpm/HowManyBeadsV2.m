%Steph 9/10

%Based off of HowManyBeads, but takes as its input a vector of looping
%probabilities, so that I can input different subsets of the pLoops from a
%given GaussFit.  So far only uses the 1st method described in
%HowManyBeads.  Optional inputs are handles to the figures to plot things in.
%Outputs are the handles to the figure (same as second and third inputs if they exist).

%10/10 Changed to sampling with replacement after talking with Martin (need
%each resampling to be independent ... )

function [fighandle1,fighandle2] = HowManyBeadsV2(pLoops,varargin)
    totmean=mean(pLoops);
    totSE = std(pLoops)/sqrt(length(pLoops)-1);
    
    %%%%% WITHOUT REPLACEMENT

	%Use randpermute to choose j random
    %pLoops from all pLoops, and do this some large number of times. 
    %Randpermute returns a random permutation of numbers from 1:n.  
    %Then just select the first 5, the first 10, etc.
    
%     incr = 2;
%     start = 4;
%     numperms = 10000; %Tried this with 100, 1000, 10000: as the number of beads included
%         %gets larger, the total number of possible permutations gets large
%         %very fast, so even though 100 length-5 permutations may be enough
%         %to see the distribution of means, 100 length-40 permutations may
%         %not sample the space enough.  100 vs. 1000 seemed to make a small
%         %difference for the spread in the means at large bead number (and
%         %of course at small), 1000 vs 10000 not so much ...
%     
%     index = 1;
%     
%     means = zeros(floor((length(pLoops)-start)/incr)+1,numperms);
%     SEs = zeros(floor((length(pLoops)-start)/incr)+1,numperms);
%     
%     for j=start:incr:length(pLoops)
%         newdistribs = zeros(numperms,j);
%         for k=1:numperms %Want numperms permutations of each length
%             newperm = randperm(length(pLoops));
%             newdistribs(k,:) = pLoops(newperm(1:j));
%         end
%         means(index,:) = mean(newdistribs,2);%Means of all the rows
%         SEs(index,:)=std(newdistribs,0,2)/sqrt(size(newdistribs,2)-1);%The function
%             %std computes according to one of two formulas, depending on
%             %whether the second argument is 0 or 1.  The default is 0, and
%             %since that's what Hernan used in calc_ploop, that's what I'm
%             %using here.  (I need to put that argument in so that I can put
%             %the dim argument in the third place).  So this takes std
%             %across the rows, according to the same formula in calc_ploop.
%             %Then to get the standard error I want to divide by the number
%             %of samples going into each SE, which means divide by the
%             %number of columns minus 1.
%         disp(strcat('Subsample:',int2str(j),'/',int2str(length(pLoops))))
%         index=index+1;
%     end
    
    %%%%%% WITH REPLACEMENT %%%%
    
    %Use randsample(n,k,1), which returns a 1-by-k vector of values sampled
    %uniformly at random from integers 1 to n; the third input is a logical
    %which says yes replacement if true.
    
    incr = 2;
    start = 4;
    numperms = 10000; %Tried this with 100, 1000, 10000: as the number of beads included
        %gets larger, the total number of possible permutations gets large
        %very fast, so even though 100 length-5 permutations may be enough
        %to see the distribution of means, 100 length-40 permutations may
        %not sample the space enough.  100 vs. 1000 seemed to make a small
        %difference for the spread in the means at large bead number (and
        %of course at small), 1000 vs 10000 not so much ...
    
    index = 1;
    
    means = zeros(floor((length(pLoops)-start)/incr)+1,numperms);
    SEs = zeros(floor((length(pLoops)-start)/incr)+1,numperms);
    bootstrpSD = zeros(floor((length(pLoops)-start)/incr)+1,1);
    
    for j=start:incr:length(pLoops)
        newdistribs = zeros(numperms,j);
        for k=1:numperms %Want numperms permutations of each length
            newperm = randsample(length(pLoops),j,1);
            newdistribs(k,:) = pLoops(newperm);
        end
        means(index,:) = mean(newdistribs,2);%Means of all the rows
        SEs(index,:)=std(newdistribs,0,2)/sqrt(size(newdistribs,2)-1);%The function
            %std computes according to one of two formulas, depending on
            %whether the second argument is 0 or 1.  The default is 0, and
            %since that's what Hernan used in calc_ploop, that's what I'm
            %using here.  (I need to put that argument in so that I can put
            %the dim argument in the third place).  So this takes std
            %across the rows, according to the same formula in calc_ploop.
            %Then to get the standard error I want to divide by the number
            %of samples going into each SE, which means divide by the
            %number of columns minus 1.
        bootstrpSD(index) = std(means(index,:));
        disp(strcat('Subsample:',int2str(j),'/',int2str(length(pLoops))))
        index=index+1;
    end
    
    %Plot the results
    if nargin==3
        subplot(varargin{1});
        fighandle1=varargin{1};
    else
        fighandle1=figure;
    end
    PlotHandle1=plot([start:incr:incr*(size(means,1)+1)],means,'ob');
    hold on
    PlotHandle1(end+1:end+2)=errorbarxyHG(length(pLoops),totmean,[],totSE,[],[],'ob','b');
    PlotHandle1(end+1:end+3)=plot([0 incr*ceil(length(pLoops)/incr)],[totmean totmean],'--k', ...
        [0 incr*ceil(length(pLoops)/incr)],[totmean+totSE totmean+totSE],'--k',...
        [0 incr*ceil(length(pLoops)/incr)],[totmean-totSE totmean-totSE],'--k');
    %plot(length(GaussFit),allpLoops,'ob')
    xlim([0 incr*ceil(length(pLoops)/incr)])
    ylim([0 1])
    xlabel('Number of Beads')
    ylabel('Mean Looping Probability')
    
    StandardFigure(PlotHandle1,gca)
    
    if nargin==3
        subplot(varargin{2});
        fighandle2=varargin{2};
    else
        fighandle2=figure;
    end
    PlotHandle2 = plot([start:incr:incr*(size(means,1)+1)],SEs,'xr');
    hold on
    PlotHandle2(end+1)=plot([0 incr*ceil(length(pLoops)/incr)],[totSE totSE],'--k');
    if length(pLoops)>3 %This is the only line that throws an error if there's <start data points ... should do a better fix than this at some point
        PlotHandle2(end+1) = plot([start:incr:incr*(size(means,1)+1)],bootstrpSD,'-g');
    end
    xlim([0 incr*ceil(length(pLoops)/incr)])
    %ylim([0 1])
    xlabel('Number of Beads')
    ylabel('Standard Error')
    
    StandardFigure(PlotHandle2,gca)
    