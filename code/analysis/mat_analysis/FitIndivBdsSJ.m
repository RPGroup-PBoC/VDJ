%Steph 11/08

%Modification of HGIndividualBeads to interface with masterscriptV2: plot
%the histograms for each bead, fit to some number of Gaussians, store the
%fit parameters.  Optional input can be output from a previous run of
%FitIndivBds.

%Update 1/2012 to run as fast with R2011a as it used to with R2007a.

function [GaussFit,type]=FitIndivBdsSJ(newr,allbeads,names,fps,unloopedlengths,varargin)

%Set up the paramters that will control how the histogram is binned and
%displayed and the start parameters for the fits:
ok=0;
while ok==0
    datatype = input('Shortloops, HGseqs, PUC, WT, WTKO, smallbds,other?','s');
    if strcmpi(datatype,'shortloops') || strcmpi(datatype,'short loops')
        xx=[50:2:200];
        startp=[0.2 100 6,0.002 120 0.5,0.1 136 5];
        TraceRange=[90,160];
        type='SL';
        ok=1;
    elseif strcmpi(datatype,'smallbds')
        xx=[50:2:200];
        startp=[0.2 60 6,0.002 75 0.5,0.1 100 5];
        TraceRange=[50,160];
        type='SB';
        ok=1;
    elseif strcmpi(datatype,'HGseqs')
        if isempty(varargin)
            GaussFit=FitIndivBds_HGseqs(newr,allbeads,names,fps,unloopedlengths);
            type='HGseqs';
        else
            GaussFit=FitIndivBds_HGseqs(newr,allbeads,names,fps,unloopedlengths,varargin{1});
            type='HGseqs';
        end
        return
    elseif strcmpi(datatype,'PUC')
        xx=[100:2:360];
        startp=[ 0.2069 146.7278 8.0346, 0.2025 168.2301 8.8809, 0.0621 198.5643 9.1675];
        TraceRange=[100,360];
        type='PUC';
        ok=1;
    elseif strcmpi(datatype,'WT')
        if isempty(varargin)
            GaussFit=FitIndivBds_wt(newr,allbeads,names,fps,unloopedlengths);
            type='WT';
        else
            GaussFit=FitIndivBds_wt(newr,allbeads,names,fps,unloopedlengths,varargin{1});
            type='WT';
        end
        return
    elseif strcmpi(datatype,'WTKO')
        xx=[80:2:240];
        startp=[0.0042 115 14.89, 0.0012 129 63.3, 0.02449 192 9.085];
        TraceRange=[80,240];
        type='WTKO';
        ok=1;
    elseif strcmpi(datatype,'other') %Should do some error-handling here ...
        xx = input('Enter x range (= start:increment:end) :');
        startp= input('Enter start parameters (six-element vector): ');
        TraceRange = input('Enter TraceRange (= [start end]) :');
        type='O';
        ok=1;
    else
        disp('Invalid input.')
    end
end
xxRange=linspace(xx(1),xx(end));
StandardRange=[xx(1),xx(end),0,0.16];

%First add up all the Gaussians--used in order to guess where each peak is.
%Fit the sum of the distributions and uses the results as starting
%parameters.  Why are they paramenters in HGIndividualBeads then ... ?
GaussFitAll=Fit23Gaussians(allbeads','All',startp,[],[],xx);

startp=GaussFitAll.startp;

if isempty(varargin)
    for i=1:length(newr)
        display([num2str(i),'/',num2str(length(newr))]);
        GaussFit(i)=Fit23Gaussians(newr{i},names{i},startp,[],GaussFitAll.fp,xx);    
    end
else
    GaussFit=varargin{1};
end



%Interactive section:
%1, 2, 3: Select the respective peaks.

PeakSelected=0; %Which peak is selected for modification
PlotRange=StandardRange; %Default range for the plot.
VerticalLine=linspace(0,0.2);
i=1;

while (i<(length(newr)+1))
    close
    %load([Folder,filesep,D(i).name])
    newxyr=newr{i};
    Time=[1:length(newxyr)]/fps;
    %figure(2)
    BeadHist = figure('Position',[150 300 950 500]);
    subplot('Position', [.51,.15,.45,.75])
    %set(gcf,'Position',[177,37,692,151])
    
    %plot(Time,newxyr,'-k')
    
    if unloopedlengths==0
        plot(Time,newxyr,'-k')
    else
        plot(Time,newxyr,'-k',Time,ones(1,length(Time)).*unloopedlengths(i),'--b');
    end
    ylim(TraceRange)

    %BeadHist=figure(1);
    
    subplot('Position', [.08,.15,.35,.75]);
    %set(gcf,'Position',[229,270,560,420]);
    
    if GaussFit(i).Approved==1
        Color='b';
    elseif GaussFit(i).Approved==0
        Color='r';
    end
    
    
    if length(GaussFit(i).startp)==9
        
        IX=[1,2,3];
        
%         IX(1)=find(strcmp(GaussFit(i).Identity,'B'));
%         IX(2)=find(strcmp(GaussFit(i).Identity,'M'));
%         IX(3)=find(strcmp(GaussFit(i).Identity,'U'));
    
        GaussParameters=[GaussFit(i).fp.a1,GaussFit(i).fp.a2,GaussFit(i).fp.a3;...
            GaussFit(i).fp.b1,GaussFit(i).fp.b2,GaussFit(i).fp.b3;...
            GaussFit(i).fp.c1,GaussFit(i).fp.c2,GaussFit(i).fp.c3];

    
        %Plots everything
        Gauss1=GaussParameters(1,IX(1))*...
            exp(-((xxRange-GaussParameters(2,IX(1)))/GaussParameters(3,IX(1))).^2);
        Gauss2=GaussParameters(1,IX(2))*...
            exp(-((xxRange-GaussParameters(2,IX(2)))/GaussParameters(3,IX(2))).^2);
        Gauss3=GaussParameters(1,IX(3))*...
            exp(-((xxRange-GaussParameters(2,IX(3)))/GaussParameters(3,IX(3))).^2);
    
   

        bar(xx,GaussFit(i).all,Color)
        hold on
        plot(xxRange,GaussFit(i).fp(xxRange),'-r',...
            xxRange,Gauss1,'--g',...
            xxRange,Gauss2,'--g',...
            xxRange,Gauss3,'--g',...
            'LineWidth',1);
        
        % Plotting vertical lines to denote RMS of looped and unlooped
        % states.
        plot([255 255],[0 1],'r')
        plot([318 318],[0 1],'r')
          
        
    elseif length(GaussFit(i).startp)==6

        IX=[1,2];
%         
%         IX(1)=find(strcmp(GaussFit(i).Identity,'B'));
%         IX(2)=find(strcmp(GaussFit(i).Identity,'M'));
%         IX(3)=find(strcmp(GaussFit(i).Identity,'U'));
        
        GaussParameters=[GaussFit(i).fp.a1,GaussFit(i).fp.a2;...
            GaussFit(i).fp.b1,GaussFit(i).fp.b2;...
            GaussFit(i).fp.c1,GaussFit(i).fp.c2];

    
        %Plots everything
        Gauss1=GaussParameters(1,IX(1))*...
            exp(-((xxRange-GaussParameters(2,IX(1)))/GaussParameters(3,IX(1))).^2);
        Gauss2=GaussParameters(1,IX(2))*...
            exp(-((xxRange-GaussParameters(2,IX(2)))/GaussParameters(3,IX(2))).^2);
       

        bar(xx,GaussFit(i).all,Color)
        hold on
        plot(xxRange,GaussFit(i).fp(xxRange),'-r',...
            xxRange,Gauss1,'--g',...
            xxRange,Gauss2,'--g',...
            'LineWidth',1);
    
      
        
    end
    
    if PeakSelected~=0
        eval(['plot(xxRange,Gauss',num2str(PeakSelected),',''--k'',''LineWidth'',2.5)'])
    end
    
    plot(GaussFitAll.fp.b1*ones(length(VerticalLine),1),VerticalLine,'--r',...
        GaussFitAll.fp.b2*ones(length(VerticalLine),1),VerticalLine,'--r',...
        GaussFitAll.fp.b3*ones(length(VerticalLine),1),VerticalLine,'--r')
    
    axis(PlotRange)
    hold off

    
   
        
    for j=1:length(IX)
        eval(['PositionPeak=GaussFit(i).fp.b',num2str(IX(j)),';']);
        text(PositionPeak+5,GaussFit(i).fp(PositionPeak)+.007,[GaussFit(i).Identity{j},...
            '/',num2str(IX(j))],'BackgroundColor',[.7 .9 .7])
    end
        
    %Steph 2/3/11 changed: somehow with the new version of masterscript
    %(v3) this no longer works.  Replaced with the next line.
%     title(['Bead ',num2str(i),'/',num2str(length(GaussFit)),...
%         ', length: ',num2str(sum(GaussFit(i).n(2:end))),', ',...
%         GaussFit(i).name],'Interpreter','none')
    
    title(strcat('Bead ',num2str(i),'/',num2str(length(GaussFit)),...
        ', length: ',num2str(sum(GaussFit(i).n(2:end))),', ',...
        GaussFit(i).name),'Interpreter','none')
    
    cc=1;
    while (cc~=13)
        ct=waitforbuttonpress;
        cc=get(BeadHist,'currentcharacter');
        cm=get(gca,'CurrentPoint');

        if ct==1
            if cc=='.'
                %Steph 10/09 added the following so that you don't get an
                %error at the very end for mislabeled peaks
                if GaussFit(i).Approved==1
                    if length(GaussFit(i).Identity)==2
                        if strcmpi(GaussFit(i).Identity{1},GaussFit(i).Identity{2}) || sum(strcmpi(GaussFit(i).Identity{1},{'B','M','U','L'}))~=1 ...
                                || sum(strcmpi(GaussFit(i).Identity{2},{'B','M','U','L'}))~=1
                            'Peaks mislabeled.'
                        else
                            i=i+1;
                        end
                    else
                        if strcmpi(GaussFit(i).Identity{1},GaussFit(i).Identity{2}) || strcmpi(GaussFit(i).Identity{1},GaussFit(i).Identity{3}) || ...
                                strcmpi(GaussFit(i).Identity{2},GaussFit(i).Identity{3}) || sum(strcmpi(GaussFit(i).Identity{1},{'B','M','U','L'}))~=1 || ...
                                sum(strcmpi(GaussFit(i).Identity{2},{'B','M','U','L'}))~=1 || sum(strcmpi(GaussFit(i).Identity{3},{'B','M','U','L'}))~=1
                            'Peaks mislabeled.'
                        else
                            i=i+1;
                        end
                    end
                else
                    i=i+1;
                end
                cc=13;
                PeakSelected=0;
            elseif cc==','
                if i>1
                    i=i-1;
                end
                cc=13;
                PeakSelected=0;
            elseif cc=='g'
                if (length(GaussFit(i).startp)==9)&(PeakSelected~=0)
                    Filter=logical(ones(1,9));
                    Filter(((PeakSelected-1)*3+1):((PeakSelected-1)*3+3))=0;
                                    
                    GaussFit(i).startp=GaussFit(i).startp(Filter);
                    GaussFit(i)=Fit23Gaussians(newr{i},GaussFit(i).name,GaussFit(i).startp,[],GaussFitAll.fp,xx);
                elseif length(GaussFit(i).startp)==6
                    GaussFit(i).startp=[GaussFit(i).startp,startp(7:9)];
                    GaussFit(i)=Fit23Gaussians(newr{i},GaussFit(i).name,GaussFit(i).startp,[],GaussFitAll.fp,xx);
                end
                PeakSelected=0;
                cc=13;
            elseif (cc=='r')&(PeakSelected~=0)
                NewIdentity=upper(input('Enter new peak identity (B,M,U): ','s'));
                GaussFit(i).Identity{PeakSelected}=NewIdentity;
                cc=13;
            elseif cc=='d'  %Discard a bead
                GaussFit(i).Approved=~(GaussFit(i).Approved);
                cc=13;
            elseif cc=='1'; %Select peak 1
                PeakSelected=1;
                cc=13;
            elseif cc=='2'; %Select peak 2
                PeakSelected=2;
                cc=13;
            elseif (cc=='3')&(length(GaussFit(i).startp)==9);  %Select peak 3
                PeakSelected=3;
                cc=13;
            elseif cc=='z'; %Zoom in
                Limits=ginput(2);
                PlotRange=[min(Limits(:,1)),max(Limits(:,1)),...
                    min(Limits(:,2)),max(Limits(:,2))];
                cc=13;
            elseif cc=='u'; %Zoom out
                PlotRange=StandardRange;
                cc=13;
            elseif cc=='f'&(PeakSelected~=0) %Refit the current peak
                display('Select peak position and width');
                NewParam=ginput(2);
                
                GaussFit(i).startp((PeakSelected-1)*3+1)=NewParam(1,2);
                GaussFit(i).startp((PeakSelected-1)*3+2)=NewParam(1,1);
                GaussFit(i).startp((PeakSelected-1)*3+3)=...
                    abs(NewParam(1,1)-NewParam(2,1));
                
                GaussFit(i)=Fit23Gaussians(newr{i},GaussFit(i).name,GaussFit(i).startp,[],GaussFitAll.fp,xx);
                
                cc=13;
            %added by Steph 6/9/08: for dealing with histograms with only
            %one peak (right now this is irreversible, kinda):
            elseif cc=='a' %Only a looped state ("_a_ll looped")
                if (length(GaussFit(i).startp)==9)%Want to just have looped and unlooped
                    Filter=logical(ones(1,9));
                    Filter(1:3)=0;%Arbitrarily getting rid of B
                                    
                    GaussFit(i).startp=GaussFit(i).startp(Filter);
                    GaussFit(i)=Fit23Gaussians(newr{i},GaussFit(i).name,GaussFit(i).startp,[],GaussFitAll.fp,xx);
                end
                GaussFit(i).p=[1,0];
                GaussFit(i).Identity={'L','U'};
                'Ploop set to 1.'
                cc=13;
            elseif cc=='n' %Only an unlooped state ("_n_one looped")
                GaussFit(i).p=[0,0,1];
                GaussFit(i).Identity={'B','M','U'};
                'Ploop set to 0.'
                cc=13;
                
            end
        end
    end
end




