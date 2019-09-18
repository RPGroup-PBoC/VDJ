%Steph 3/4/09

%Function based on plot_ploop.m that's called from singlebdanalysisV2 to
%calculate ploop and SEpLoop for each looped state in a three-state or
%five-state system (e.g. PUC/shortloops/WTKO vs. WT).  Five-state systems
%take WT=1 as their input, three-state systems have WT=0.

function [pLoop,SEpLoop,pB,SEpB,pM,SEpM,varargout]=calc_ploop(GaussFit,WT)

if WT==0
    i=1; %Artifact from when this code was part of a script that would calculate ploop for
            %for multiple GaussFits--should clean this up at some point
    
    pB=[];
    pM=[];
    pL=[];
    pU=[];

    GaussFit=GaussFit(logical([GaussFit.Approved]));%converts the array of approved 1's and zeros in GaussFit to logicals (trues and falses);
            %somehow this removes any beads whose approved fields were
            %0, meaning that the size of the structure GaussFit is
            %reduced
    %Steph 8/10: Note that this works even for cases where pLoop = 0 or 1,
    %because FitIndivBds assigns Identity as [L,U] if pLoop = 1, where
    %GaussFit.p = [1,0], and Identity as [B, M, U] if ploop = 0, where
    %GaussFit.p = [0,0,1].
    for j=1:length(GaussFit)%for all the beads in GaussFit
        strcat('Bead ',int2str(j),'/',int2str(length(GaussFit))) %prints to tell you how many you've looked at so far
        if sum(strcmp(GaussFit(j).Identity,'B'))
            SeqFit(i).pB(j)=GaussFit(j).p(strcmp(GaussFit(j).Identity,'B'));%strcmp=string compare; 
        else
            SeqFit(i).pB(j)=0;
        end
        
        if sum(strcmp(GaussFit(j).Identity,'M'))
            SeqFit(i).pM(j)=GaussFit(j).p(strcmp(GaussFit(j).Identity,'M'));
        else
            SeqFit(i).pM(j)=0;
        end
                
        if sum(strcmp(GaussFit(j).Identity,'U'))
            SeqFit(i).pU(j)=GaussFit(j).p(strcmp(GaussFit(j).Identity,'U'));
        else
            SeqFit(i).pU(j)=0;
        end
        
        if sum(strcmp(GaussFit(j).Identity,'L'))
            SeqFit(i).pL(j)=GaussFit(j).p(strcmp(GaussFit(j).Identity,'L'));
        else
            SeqFit(i).pL(j)=0;
        end
        
    end
    
    %In the code for calculating pLoop for multiple GuassFits at once
    %there's a loop over length(SeqFit) here, but since in this case
    %length(SeqFit) is always 1 ...
    for j=1:length(SeqFit(1).pB)
        if (SeqFit(1).pL(j)==0)&(SeqFit(1).pM(j)==0)&(SeqFit(1).pB(j)==0)
            SeqFit(1).pLoop(j)=0;
        elseif SeqFit(1).pL(j)==0
            SeqFit(1).pLoop(j)=SeqFit(1).pB(j)+SeqFit(1).pM(j);
        else
            SeqFit(1).pLoop(j)=SeqFit(1).pL(j);
        end
    end
    pLoop=mean(SeqFit(1).pLoop);
    SEpLoop=std(SeqFit(1).pLoop)/sqrt(length(SeqFit(1).pLoop)-1);
    pB=mean(SeqFit(1).pB);
    SEpB=std(SeqFit(1).pB)/sqrt(length(SeqFit(1).pB)-1);
    pM=mean(SeqFit(1).pM);
    SEpM=std(SeqFit(1).pM)/sqrt(length(SeqFit(1).pM)-1);
        
%For WT:    
else
    %For wtlac, there are 9 possible peak labels: O12, O12u, O12b, O23, O23u,
    %O23b, O13, l, u; only going to separate out O12+O12u+O12b, O23+O23u+O23b,
    %O13,l, and u, but need to extract them all separately:
    pO12=[];
    pO12u=[];
    pO12b=[];
    pO23=[];
    pO23u=[];
    pO23b=[];
    pO13=[];
    pL=[];
    pU=[];

    GaussFit=GaussFit(logical([GaussFit.Approved]));%converts the array of approved 1's and zeros in GaussFit to logicals (trues and falses);
            %somehow this removes any beads whose approved fields were
            %0, meaning that the size of the structure GaussFit is
            %reduced
    for j=1:length(GaussFit)%for all the beads in GaussFit
        strcat('Bead ',int2str(j),'/',int2str(length(GaussFit))) %prints to tell you how many you've looked at so far
        if sum(strcmpi(GaussFit(j).Identity,'U')) %Means there's a 'U' somewhere in the Identity matrix
            SeqFit(i).pU(j)=GaussFit(j).p(strcmpi(GaussFit(j).Identity,'U'));%strcmp=string compare; 
        else
            SeqFit(i).pU(j)=0;
        end
        
        if sum(strcmpi(GaussFit(j).Identity,'O13')) 
            SeqFit(i).pO13(j)=GaussFit(j).p(strcmpi(GaussFit(j).Identity,'O13'));%strcmp=string compare; 
        else
            SeqFit(i).pO13(j)=0;
        end
        
        if sum(strcmpi(GaussFit(j).Identity,'L')) 
            SeqFit(i).pL(j)=GaussFit(j).p(strcmpi(GaussFit(j).Identity,'L'));%strcmp=string compare; 
        else
            SeqFit(i).pL(j)=0;
        end
        
        if sum(strcmpi(GaussFit(j).Identity,'O12')) 
            SeqFit(i).pO12(j)=GaussFit(j).p(strcmpi(GaussFit(j).Identity,'O12'));%strcmp=string compare; 
        else
            SeqFit(i).pO12(j)=0;
        end
        
        if sum(strcmpi(GaussFit(j).Identity,'O12u')) 
            SeqFit(i).pO12u(j)=GaussFit(j).p(strcmpi(GaussFit(j).Identity,'O12u'));%strcmp=string compare; 
        else
            SeqFit(i).pO12u(j)=0;
        end
        
        if sum(strcmpi(GaussFit(j).Identity,'O12b')) 
            SeqFit(i).pO12b(j)=GaussFit(j).p(strcmpi(GaussFit(j).Identity,'O12b'));%strcmp=string compare; 
        else
            SeqFit(i).pO12b(j)=0;
        end
        
        if sum(strcmpi(GaussFit(j).Identity,'O23')) 
            SeqFit(i).pO23(j)=GaussFit(j).p(strcmpi(GaussFit(j).Identity,'O23'));%strcmp=string compare; 
        else
            SeqFit(i).pO23(j)=0;
        end
        
        if sum(strcmpi(GaussFit(j).Identity,'O23u')) 
            SeqFit(i).pO23u(j)=GaussFit(j).p(strcmpi(GaussFit(j).Identity,'O23u'));%strcmp=string compare; 
        else
            SeqFit(i).pO23u(j)=0;
        end
        
        if sum(strcmpi(GaussFit(j).Identity,'O23b')) 
            SeqFit(i).pO23b(j)=GaussFit(j).p(strcmpi(GaussFit(j).Identity,'O23b'));%strcmp=string compare; 
        else
            SeqFit(i).pO23b(j)=0;
        end
        
    end

    for j=1:length(SeqFit(i).pU) %Length of SeqFit(i).pWhatever = number of beads 
        %Don't completely understand Hernan's code here--need to check
        %this:
        SeqFit(i).pLoop(j)=SeqFit(i).pL(j)+SeqFit(i).pO13(j)+SeqFit(i).pO12(j)+...
            SeqFit(i).pO12b(j)+SeqFit(i).pO12u(j)+SeqFit(i).pO23(j)+SeqFit(i).pO23u(j)+...
            SeqFit(i).pO23b(j);

        SeqFit(i).pO12loop(j)=SeqFit(i).pO12(j)+SeqFit(i).pO12u(j)+SeqFit(i).pO12b(j);

        SeqFit(i).pO23loop(j)=SeqFit(i).pO23(j)+SeqFit(i).pO23u(j)+SeqFit(i).pO23b(j);

        SeqFit(i).pOtherloop(j)=SeqFit(i).pO13(j)+SeqFit(i).pL(j);

    end
    pLoop=mean(SeqFit(1).pLoop);
    SEpLoop=std(SeqFit(1).pLoop)/sqrt(length(SeqFit(1).pLoop)-1);
    pO12loop=mean(SeqFit(1).pO12loop);
    SEO12=std(SeqFit(1).pO12loop)/sqrt(length(SeqFit(1).pO12loop)-1);
    pO23loop=mean(SeqFit(1).pO23loop);
    SEO23=std(SeqFit(1).pO23loop)/sqrt(length(SeqFit(1).pO23loop)-1);
    pOtherloop=mean(SeqFit(1).pOtherloop);
    SEotherloop=std(SeqFit(1).pOtherloop)/sqrt(length(SeqFit(1).pOtherloop)-1);
    
    %To get the outputs right:
    pM=pO12loop;
    SEpM=SEO12;
    pB=pO23loop;
    SEpB=SEO23;
    varargout{1}=pOtherloop;
    varargout{2}=SEotherloop;
end



