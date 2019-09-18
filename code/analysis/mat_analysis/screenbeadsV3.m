%function [screened_r,record,cutpts] = screenbeadsV3(nolacr,lacr,fps,nolacname, lacname, ...
%    path,savedirname,contourlength)
%
%Allows user to screen TPM traces: beads can be kept, discarded, or
%discarded after a certain point.  User can watch the beads' movies for any
%part of the trace.
%
%Inputs are:
%nolacr: a beads-by-frames matrix of nolac data
%lacr: a beads-by-frames matrix of lac data; same number of rows as nolacr
%fps: camera frame rate
%nolacname: set name for no lac loaded from masterscriptV2params
%lacname: set name for lac loaded from masterscriptV2params
%path: full path to where the .pxl movie files are
%savedirname: where to store the images of the traces
%contourlength: contour length of the unlooped DNA in bp; for guides to the
%   eye
%
%Stephanie Johnson 1/11

function [screened_r,record,cutpts] = screenbeadsV3(nolacr,lacr,fpsL,fpsNL, ...
    nolacname, lacname,path,savedirname,contourlength,wind,clA)

%%%% Parameters for making the plots nice %%%%
%Plot titles: Replace any underscores in nolacname and lacname with \_, so they appear
%nicely in the plot titles
lactitle=strrep(lacname,'_','\_');
nolactitle=strrep(nolacname,'_','\_');
%For plotting guides to the eyes, based on Soichi's calibration curves
% Resurgam
unlooped = -2.00*10^(-5)*contourlength^2+0.143*contourlength+92;
% Fluxus
%unlooped = -1.85*10^(-5)*contourlength^2+0.136*contourlength+92;
%{
States for truncated DNA
singly = 262.8;
doubly = 258.2;
looped = 210;
%}
%bent = 150;

%{
States for 2940 bp 12/23 RSS DNA, 40nM HMGB1
unlooped_hmgb1 = 305;
mid = 285;
looped = 260;
%}

%States for 2940bp 12/23RSS DNA, 80nM HMGB1
unlooped_hmgb1 = 318;
mid = 285;
looped = 255;
%}

%{
States for 2940bp 12/23RSS DNA, 80nM HMGB1, Ca2+
unlooped_hmgb1 = 300;
mid = 265;
looped = 250;
%}

%{
States for 2940 bp 12/23 RSS DNA, 160nM HMGB1
unlooped_hmgb1 = 270;
mid = 250;
looped = 230;
%}

%mid=310;
%looped=280;
%lppped = 197;
%DL = [500 1000 1500 2000 2500 3000];

%Get x-axes for plots: *1 is for nolac, *2 is for lac
[bds L2]=size(lacr);
L1=size(nolacr,2);
secs1=L1/fpsNL; %Total seconds of nolac trace
x1=(1/fpsNL):(1/fpsNL):secs1; %x-axis, in seconds, of nolac trace
secs2=L2/fpsL; %likewise for lac trace
x2=(1/fpsL):(1/fpsL):secs2;
%binvect = 130:1:210;
binvect = 155:1:380; %This is for the histogram
%binvect = 0:1:380;
%To keep track of user screening:
record=zeros(1,bds); %for compatibility with old code: If a bead is discarded (completely), record at that index will be 0, otherwise 1
cutpts = zeros(1,bds);%for each bead, is 0 if the entire bead is discarded, is size(lacr,2) if a bead is
    %kept completely, and is the frame number after which the bead was
    %discarded if part of it is cut

b = 1; %this indexes which bead is being examined.
    
while b <= bds %Loops until user has gone through all beads in this set

    bd1=nolacr(b,:);
    bd2=lacr(b,:);

    %Plot a trace for the user
    figure('Position',[200 200 900 600])
        %Plot the nolac trace
        subplot('Position', [.1,.6,.6,.3])
        title(strcat(nolactitle,'\_Bead',int2str(b)),'Fontsize',16);
        hold on
        %ylim([110 binvect(end)])
        ylim([binvect(1) binvect(end)]) %sets y axis limits for all plots
        % Truncated DNA length
        %ylim([0 binvect(end)])
        xlabel('Time (sec)','Fontsize',20);
        %ylabel('DNA length (bp)','Fontsize',20)
        ylabel('<RMS> (nm)','Fontsize',20)
        set(gca,'FontSize',14)
        plot(x1,bd1)
        plot([0 secs1],[unlooped unlooped],'k-')
        %set(gca,'YTick',[160 170 180 190]);
        %set(gca,'YTick',[binvect(1) 210 257 295 324 345 binvect(end)]);
        %set(gca, 'YTick',[200 looped doubly singly unlooped 280]);
        set(gca, 'YTick', [binvect(1) unlooped binvect(end)]);
        %set(gca,'YTickLabel',DL)
        set(gca,'FontSize',20)
        
        %plot([0 secs1],[unlop unlop],'k-')%Plotting "guides for the eyes"
        xlim([0 secs1])
        hold off
        %Plot the nolac histogram
        subplot('Position', [.71,.6,.2,.3])
        [n, xout] = hist(bd1, binvect);
        prob = n./(sum(n)); 
        prob(1) = 0; 
        barh(xout, prob)
        %ylim([110 binvect(end)])
        ylim([binvect(1) binvect(end)])
        % Truncated DNA screen
        %ylim([200 binvect(end)])
        set(gca,'YTick',[]);
        %set(gca,'YTickLabel',DL)
        set(gca,'XTick',[0 0.1 0.2]);
        xlim([0 0.2])
        xlabel('Probability','Fontsize',20)
        set(gca,'FontSize',20)
        
    	%Plotting the RAG and HMGB1 contributed trace:
        h=subplot('Position', [.1,.13,.6,.3]);
        title(strcat(lactitle,'\_Bead',int2str(b)),'Fontsize',16);
        hold on
        %ylim([110 binvect(end)])
        %ylim([155 binvect(end)])
        ylim([binvect(1) binvect(end)])
        xlabel('Time (sec)','Fontsize',20)
        %ylabel('DNA length (bp)','Fontsize',20)
        ylabel('<RMS> (nm)','Fontsize',20)
        set(gca,'FontSize',20)
        plot(x2,bd2)
        %
        plot([0 secs2],[unlooped_hmgb1 unlooped_hmgb1],'k--')
        plot([0 secs2],[mid mid], 'b-')
        plot([0 secs2],[looped looped],'r')
        %plot([0 secs2],[unlooped unlooped],'k-')%Plotting "guides for the eyes"
        %plot([0 secs2],[bent bent],'r')
        %{
        plot([0 secs2],[singly singly],'b')
        plot([0 secs2],[doubly doubly],'k--')
        %}
        %set(gca,'YTick',[160 170 180 190]);
        % Plot the 1 hour mark
        plot([3600 3600],[0 400],'r:')
        set(gca,'YTick',[binvect(1) binvect(end)]);
        % Truncated DNA
        set(gca, 'YTick',[200 looped mid unlooped_hmgb1 380]);
        %set(gca, 'YTick', [200 looped mid unlooped 380]);
        %set(gca, 'YTick',[bent unlooped]);
        %set(gca,'YTickLabel',DL)
        xlim([0 secs2])
        hold off
        %Plotting the lac histogram
        subplot('Position', [.71,.13,.2,.3])
        [n2, xout2] = hist(bd2, binvect); 
        prob2 = n2./(sum(n2)); 
        prob2(1) = 0; %Because screenbeads sets 'bad' data to zero, the first bin will contain all these points
        barh(xout2, prob2)
        %ylim([110 binvect(end)])
        %ylim([155 binvect(end)])
        ylim([binvect(1) binvect(end)])
        set(gca,'YTick',[]);
        set(gca,'XTick',[0 0.1 0.2]);
        xlim([0 0.2])
        xlabel('Probability','Fontsize',20)
        set(gca,'FontSize',20)
        %[1 11 14 22]  [1 5 7 17 26]  [4 6 12 16 19 21 29 30]  [5 7 8 9 16 19 21 22] 
        %[2 5 8 10 20 21 28]  [2 3 7 11 12 19 23]  [6 7 8 13 19]  [1 5 9 14 18]
    %Screen data
    examine_data = 1;

    while examine_data == 1
        pxl_index = input('Enter index to watch PXL file, 0 to decide, -1 to find breakpoint, -2 for prev. bd: ');
        if pxl_index == 0
            examine_data = 0;
            reply  = input('Keep all data (enter), discard bead (1), discard data (2): ');
    
            if reply == 1; %Note that the default for record and cutpts is discarded, since it's made of zeros
                title(h,strcat(lactitle, '\_Bead ',int2str(b),' DISCARDED '),'Fontsize',16)
            elseif reply == 2;
                ditch = 0;
                while ditch == 0
                    ditchind = input('From what index would you like to discard? '); %set number at which to start ignoring data
                    if ditchind*fpsL > size(lacr,2) || ditchind < 1
                        disp('Index out of bounds.')
                    else
                        ditch = 1;
                    end
                end
                cutpts(b)=floor(ditchind*fpsL);%User will enter index to discard from in seconds, need to convert to frames
                title(h,strcat(lactitle, '\_Bead ',int2str(b),' DISCARDED AFTER ',int2str(ditchind),' Seconds'),'Fontsize',16)
                record(b)=1;
            else
                cutpts(b)=size(lacr,2);
                record(b)=1;
            end 
            
            %SAVE FIGURE
            saveasname=fullfile(savedirname,strcat(lacname,'Bead',int2str(b)));
            print('-depsc',saveasname)
            clear saveasname
            close   

        elseif pxl_index*fpsL > size(lacr,2)
            disp('Requested index out of bounds.')
        elseif isempty(pxl_index) %If the user hits enter by accident, just go back to the beginning of the loop
            examine_data = 1;
        elseif pxl_index == -1
            try 
                [breakpoint,xvalues,yvalues]=find_bead_breakV3(path,lacname,b,fpsL); %find_bead_break returns breakpoint as a frame number
                usrinput = input('Discard data after this breakpoint? (y/n)','s');
                if strcmpi(usrinput,'y')
                    examine_data=0;
                    cutpts(b)=floor(breakpoint);
                    title(h,strcat(lactitle, '\_Bead ',int2str(b),' DISCARDED AFTER ',int2str(breakpoint/fpsL),' Seconds'),'Fontsize',16)
                    record(b)=1;

                    %SAVE FIGURE
                    saveasname=fullfile(savedirname,strcat(lacname,'Bead',int2str(b)));
                    print('-depsc',saveasname)
                    clear saveasname
                    close
                end
            catch
                disp('No breakpoint found.')
            end
        elseif pxl_index == -2
            b = b-2; %because there's a b+1 command at the end of the while loop
            examine_data = 0;
            close
        elseif isnumeric(pxl_index)
            pxl_index = pxl_index*fpsL;%Traces are plotted vs. seconds, not index--convert back so can watch movie
            pxl_index=ceil(pxl_index/wind);
            figure, open_pxlV3(path,strcat(lacname,'_',int2str(pxl_index),'.pxl'),clA);
            close
        else
            disp('Invalid input.')
        end
    end
   
   if b<1
       b=1;
   else
       b=b+1;
   end
   
end

%After the user has finished analyzing all the data, create the screened_r
%matrix that will be returned

kept = 1; %this indexes the screened_r matrix

for i = 1:bds 
    if cutpts(i) == size(lacr,2) %Keep all the data for this bead
        screened_r(kept, :) = lacr(i, :);
        kept = kept+1;
    elseif cutpts(i) ~= 0 %discard after a certain point for these bead
        screened_r(kept, :) = lacr(i, :);
        screened_r(kept,cutpts(i):end) = 0; %Sets discarded data to 0
        kept = kept+1;
    end
end

