%makeallreptracesV3(dirname,parentdir,allsets,filenames,nolacallsets,...
%    nolacnames,corres3,GaussFit,fps,datatype)
%
%Makes figures that have nolac traces on the top and lac tracces on the bottom.
%The lac trace has a dashed line corresponding to the nolac length.  Both
%also have histograms to the right, the lac histo with the fit from
%GaussFit.
%
%Inputs are:
%dirname: name of a new directory where it will put the eps figure files.
%
%Stephanie Johnson 2/11

function makeallreptracesV3(dirname,parentdir,allsets,filenames,nolacallsets,...
    nolacnames,corres3,GaussFit,fps,datatype)

%Set up the paramters that will control how the histogram is binned and
%displayed, taken from FitIndivBds:
if strcmpi(datatype,'SL')
	binvect=[50:2:200];
	TraceRange=[50,180];
elseif strcmpi(datatype,'SB')
	binvect=[50:2:200];
	TraceRange=[50,180];
elseif strcmpi(datatype,'PUC')
	binvect=[100:2:360];
	xx=[100:2:360];
	TraceRange=[100,360];
elseif strcmpi(datatype,'WT') || strcmpi(datatype,'WTKO')
	binvect=[80:2:240];
	xx=[80:2:240];
	TraceRange=[80,240];
elseif strcmpi(datatype,'O') %Should do some error-handling here ...
	xx = input('Enter x range (= start:increment:end) :');
    binvect = xx;
	startp= input('Enter start parameters (six-element vector): ');
	TraceRange = input('Enter TraceRange (= [start end]) :');
end

fits=1;

numsets=length(allsets); %If there were multiple sets taken on the same area, some of these
    %sets will actually be empty
setnum=1;
bdnum=1; %Indexes the GaussFit array, which has all the beads together in one structure

mkdir(fullfile(parentdir,dirname));
savedir = fullfile(parentdir,dirname);

for i=1:numsets
    thisset=allsets{i};
    if ~isempty(thisset)
        thisnolacset=nolacallsets{i};
        thiscorres=corres3{i};
        L1=size(thisnolacset,2);
        
        [numbds L2]=size(thisset);
        name=strrep(filenames(i),'_','\_');
        nolacname=strrep(nolacnames(i),'_','\_');

        %For the nolac trace:
        secs1=L1/fps;%To plot in terms of time, not frame index
        x1=(1/fps):(1/fps):secs1;%first corresponds to 1/30 secs b/c 30fps, last is L/30
        %For the lac trace:
        secs2=L2/fps;%To plot in terms of time, not frame index
        x2=(1/fps):(1/fps):secs2;%first corresponds to 1/30 secs b/c 30fps, last is L/30

        for b=1:numbds
            %bd1 will be the nolac, bd2 the lac
            if ~isempty(nolacallsets{i})
               bd1=thisnolacset(thiscorres(b),:);
               [n, xout] = hist(bd1, binvect);
               prob = n./(sum(n(2:end))); %normalize by the total number of counts
               prob(1) = 0; %Because screenbeads sets 'bad' data to zero, the first bin will contain all these points
               opts = fitoptions('gauss1','Algorithm','Trust-Region');
               try
                    unloopedfit = fit(xout(2:end)',prob(2:end)','gauss1',opts);
                    unlooped = unloopedfit.b1; %This is where the unlooped state should be for this particular bead
               catch
                   unlooped = 0;
                   disp(strcat('Bd ',int2str(b),': Unable to fit, user-input contour length plotted instead.'))
               end
            else
               bd1=zeros(1,L1);
            end
            
            bd2=thisset(b,:);
   
            figure('Position',[200 200 900 600])
                subplot('Position', [.1,.6,.6,.3])
                title(strcat(nolacname,'\_Bead',int2str(b)),'Fontsize',16);
                hold on
                ylim(TraceRange)
                xlabel('Time (sec)','Fontsize',16);
                ylabel('<R> (nm)','Fontsize',16)
                set(gca,'FontSize',14)
                
                plot(x1,bd1)
                plot([0 L1],[unlooped unlooped],'k--')%Plotting "guides for the eyes"
                xlim([0 secs1])
                hold off

            subplot('Position', [.71,.6,.2,.3])
                barh(xout, prob)
                ylim(TraceRange)
                set(gca,'YTick',[]);
                set(gca,'XTick',linspace(0,max(prob)+0.01,3));
                xlim([0 max(prob)+0.01])
                xlabel('Probability','Fontsize',16)
                set(gca,'FontSize',14)
            
            %Plot the lac trace:
            subplot('Position', [.1,.13,.6,.3])
                title(strcat(name,'\_Bead',int2str(b)),'Fontsize',16);
                hold on
                ylim(TraceRange)
                xlabel('Time (sec)','Fontsize',16)
                ylabel('<R> (nm)','Fontsize',16)
                set(gca,'FontSize',14)

                plot(x2,bd2)
                plot([0 L2],[unlooped unlooped],'k--')
                xlim([0 secs2])
                hold off

            subplot('Position', [.71,.13,.2,.3])
                [n2, xout2] = hist(bd2, binvect);
                prob2 = n2./(sum(n2(2:end))); %normalize by the total number of counts
                prob2(1) = 0; %Because screenbeads sets 'bad' data to zero, the first bin will contain all these points

                if fits==1
                    hold on
                    if GaussFit(bdnum).Approved==1
                        if isequal(GaussFit(bdnum).p, [0 0 1])
                            barh(xout2,prob2)
                            title('p_{loop}=0','Fontsize',14)
                            ylim(TraceRange)
                            set(gca,'YTick',[]);
                            set(gca,'XTick',linspace(0,max(prob2)+0.01,3));
                            xlim([0 max(prob2)+0.01])
                            xlabel('Probability','Fontsize',16)
                            set(gca,'FontSize',14)
                        elseif isequal(GaussFit(bdnum).p,[1 0])
                            barh(xout2,prob2)
                            title('p_{loop}=1','Fontsize',14)
                            ylim(TraceRange)
                            set(gca,'YTick',[]);
                            set(gca,'XTick',linspace(0,max(prob2)+0.01,3));
                            xlim([0 max(prob2)+0.01])
                            xlabel('Probability','Fontsize',16)
                            set(gca,'FontSize',14)
                        else
                            bar(xout2,prob2)
                            if length(GaussFit(bdnum).p)==3
                                Gauss1=GaussFit(bdnum).fp.a1*...
                                    exp(-((xout2-GaussFit(bdnum).fp.b1)./GaussFit(bdnum).fp.c1).^2);
                                Gauss2=GaussFit(bdnum).fp.a2*...
                                    exp(-((xout2-GaussFit(bdnum).fp.b2)./GaussFit(bdnum).fp.c2).^2);
                                Gauss3=GaussFit(bdnum).fp.a3*...
                                    exp(-((xout2-GaussFit(bdnum).fp.b3)./GaussFit(bdnum).fp.c3).^2);
                                plot(xout2,Gauss1,'-r',...
                                    xout2,Gauss2,'-r',...
                                    xout2,Gauss3,'-r',...
                                    'LineWidth',1.5);
                            else 
                                Gauss1=GaussFit(bdnum).fp.a1*...
                                    exp(-((xout2-GaussFit(bdnum).fp.b1)./GaussFit(bdnum).fp.c1).^2);
                                Gauss2=GaussFit(bdnum).fp.a2*...
                                    exp(-((xout2-GaussFit(bdnum).fp.b2)./GaussFit(bdnum).fp.c2).^2);
                                plot(xout2,Gauss1,'-r',...
                                    xout2,Gauss2,'-r',...
                                    'LineWidth',1.5);
                            end
                            view(90,-90)
                            xlim(TraceRange)
                            set(gca,'XTick',[]);
                            set(gca,'YTick',linspace(0,max(prob2)+0.01,3));
                            ylim([0 max(prob2)+0.01])
                            ylabel('Probability','Fontsize',16)
                            set(gca,'FontSize',14)
                        end
                    else
                        barh(xout2,prob2,'r')
                        ylim(TraceRange)
                        set(gca,'YTick',[]);
                        set(gca,'XTick',linspace(0,max(prob2)+0.01,3));
                        xlim([0 max(prob2)+0.01])
                        xlabel('Probability','Fontsize',16)
                        set(gca,'FontSize',14)
                    end
                    hold off
                else
                    barh(xout2, prob2)
                    ylim(TraceRange)
                    set(gca,'YTick',[]);
                    set(gca,'XTick',linspace(0,max(prob2)+0.01,3));
                    xlim([0 max(prob2)+0.01])
                    xlabel('Probability','Fontsize',16)
                    set(gca,'FontSize',14)
                end
            bdnum=bdnum+1;
            
            saveasname=strcat('Set',int2str(setnum),'Bead',int2str(b));
            print('-depsc',fullfile(savedir,saveasname))
            saveas(gca,fullfile(savedir,saveasname),'fig') %Also save as a matlab fig file so can edit later
            clear saveasname

            %pause
            close
        
        end
        setnum=setnum+1;
    end
    
    strcat('Set',int2str(i),'/',int2str(numsets),' Finished')  
    
end
