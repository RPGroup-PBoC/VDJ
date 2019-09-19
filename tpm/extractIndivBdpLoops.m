%[allpLoops, totmean, totSE]=extractIndivBdpLoops(GaussFit)
%
%Extracts each bead's looping probability from the GaussFit structure saved
%by masterscript.  If the optional input is 'all', varargout will contain
%pLoops vectors for the bottom and middle looped states (allpLoops is
%always the total looping probability).  Specifically varargout{1} is
%bottom, varargout{2} is middle.
%
%Steph 1/11

function [allpLoops, totmean, totSE, varargout]=extractIndivBdpLoops(GaussFit,varargin)

    GaussFit=GaussFit(logical([GaussFit.Approved]));

    %Put all the beads' pLoops into one vector
    allpLoops = zeros(1,length(GaussFit));
    allBpLoops = zeros(1,length(GaussFit));
    allMpLoops = zeros(1,length(GaussFit));

    for i=1:length(GaussFit)
        tempp = 0;
        tempB = 0;
        tempM = 0;
        if sum(strcmp(GaussFit(i).Identity,'B'))
            tempp=tempp+GaussFit(i).p(strcmp(GaussFit(i).Identity,'B')); 
            tempB = tempB+GaussFit(i).p(strcmp(GaussFit(i).Identity,'B')); 
        end
        
        if sum(strcmp(GaussFit(i).Identity,'M'))
            tempp=tempp+GaussFit(i).p(strcmp(GaussFit(i).Identity,'M'));
            tempM = tempM+GaussFit(i).p(strcmp(GaussFit(i).Identity,'M'));
        end
        
        if sum(strcmp(GaussFit(i).Identity,'L'))
            tempp=tempp+GaussFit(i).p(strcmp(GaussFit(i).Identity,'L'));
        end
        
        allpLoops(i) = tempp;
        allBpLoops(i) = tempB;
        allMpLoops(i) = tempM;
        
    end
    
    totmean = mean(allpLoops);
    totSE = std(allpLoops)/sqrt(length(allpLoops)-1);
    
    if nargin>1 && strcmpi(varargin{1},'all')
        varargout{1}=allBpLoops;
        varargout{2}=allMpLoops;
    elseif nargin>1
        disp('Invalid second input.')
        varargout{1}=0;
    end
end