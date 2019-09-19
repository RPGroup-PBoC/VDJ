function ragrates_hairpin_2
%{
Processes trajectories of beads to determine looped and unlooped
states as well as times spent in each state.
%}
clear all
f=275;
path = '/Users/soichihirokawa/Documents/vdj_recombination/analysis';
specific_path = '/Analyzed Data/new_analysis/coding_flank/';

all_mutations = dir(fullfile(path,specific_path));

for nn=1:length(all_mutations)
    
    if all_mutations(nn).name(1) ~= '1' | all_mutations(nn).name(end-3:end)=='.pdf'
        continue
    end
    all_mutations(nn).name
    
    underscore_index = findstr(all_mutations(nn).name,'_');
    mutant_name = all_mutations(nn).name(underscore_index(1)+1:underscore_index(2)-1);
    
    names={{all_mutations(nn).name}};
    for fff=1:length(names)
    currexp=names{fff};
    currexp=currexp{end};
    %load data file
    load( fullfile(path, specific_path, currexp,...
        '/lacdata.mat'),'alllacnames','lacANAL_RMS','record')
    load( fullfile(path, specific_path, currexp,...
        '/nolacdata.mat'),'nolacANAL_RMS','record')

    n = 0;
    for i = 1:length(lacANAL_RMS)
        n = n + size(lacANAL_RMS{i},1);
    end

    %initialization
    statediagram=[];
    gaus_signal=[];
    b=[];
    temp=[];
    numbeads=0;
    globalmax=0;
    t=0;

    threshold=f;
    
    %{80nM HMGB1
    trackthresh=350;
    stickthresh=240;

    lowstateaverage=260;
    highstateaverage=316;
    %}
    %{
    40nM HMGB1
    trackthresh=350;
    stickthresh=220;

    lowstateaverage=265;
    highstateaverage=320;
    %}
    %A gaussian filter
    sigma = 240;
    deadfilter=2*sigma*1.34;
    filtsize = 4*sigma;
    x = linspace(-filtsize / 2, filtsize / 2, filtsize);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize

    %zero sample specific stuff
            onstate_temp=[];
            offstate_temp=[];
            onstate_nolac_temp=[];
            offstate_nolac_temp=[];

            offstate_comp=[];
            onstate_comp=[];
            offstate_nolac_comp=[];
            onstate_nolac_comp=[];        

            allstates_comp=[];
            allstates_nolac_comp=[];
            statetrace_comp=[];
            statetrace_nolac_comp=[];
            trace_comp=[];
            trace_nolac_comp=[];

            % Include length of time that bead moves in the experiment
            % While beads that get stuck and become unstuck during course of
            % the experiment will count in the experiment, we will remove beads
            % that become stick for the remainder of the experiment or float
            % away either from passive bead loss or RAG cleavage.
            time_comp = [];

    %loop over each cell in the data file (for each position)
    for l=1:length(lacANAL_RMS)
        %take the rms data
        a=lacANAL_RMS{l};
        b=nolacANAL_RMS{l};
        r=find(record{l}==1);

          % clear things
            trace=[];
            segmented=[];
            nolactrace=[];
            nolacsegmented=[];
        %gaussian filter
        for ff=1:size(a,1)
            t=t+1;
            trace=real(a(ff,3*sigma:end-3*sigma));
            l;
            ff;
            [statetrace,allstates]=segmenttrace_hairpins(trace,threshold,trackthresh,stickthresh,deadfilter);
            segmented=statetrace;
            segmented(segmented==3)=highstateaverage;
            segmented(segmented==2)=lowstateaverage;

            nolactrace=real(b(r(ff),1:end));
            [nolacstatetrace,nolacallstates]=segmenttrace_hairpins(nolactrace,threshold,trackthresh,stickthresh,deadfilter);
            nolacsegmented=nolacstatetrace;
            nolacsegmented(nolacsegmented==3)=highstateaverage;
            nolacsegmented(nolacsegmented==2)=lowstateaverage;

            %deconvolve states final
            onstate_temp=zeros(length(allstates),1);
            offstate_temp=zeros(length(allstates),1);
            onstate_nolac_temp=zeros(length(nolacallstates),1);
            offstate_nolac_temp=zeros(length(nolacallstates),1);

            onstate_temp(statetrace==2)=1;
            offstate_temp(statetrace==3)=1;
            onstate_nolac_temp(nolacstatetrace==2)=1;
            offstate_nolac_temp(nolacstatetrace==3)=1;

            offstate_comp{l}(ff,:)=offstate_temp;
            onstate_comp{l}(ff,:)=onstate_temp;
            offstate_nolac_comp{l}(ff,:)=offstate_nolac_temp;
            onstate_nolac_comp{l}(ff,:)=onstate_nolac_temp;        

            allstates_comp{l}(ff,:)=allstates;
            allstates_nolac_comp{l}(ff,:)=nolacallstates;
            statetrace_comp{l}(ff,:)=statetrace;
            statetrace_nolac_comp{l}(ff,:)=nolacstatetrace;
            trace_comp{l}(ff,:)=trace;
            trace_nolac_comp{l}(ff,:)=nolactrace;

            % Take lacANAL_RMS and segment to 0s and 1s where the threshold is
            % set to the stickthresh RMS level.
            rag_RMS_thresholded = (a(ff,:) > stickthresh);
            % Find last 1 in list and note index. This is the time when the
            % bead is declared stuck or floats away for the remainder of the 
            % experiment. We then divide by the frame rate to get the time in
            % seconds.
            if isempty(find(rag_RMS_thresholded,1,'last'))
                time_comp{l}(ff,:) = 0;
            else
                time_comp{l}(ff,:) = find(rag_RMS_thresholded,1,'last')./30;
            end


%{
             figure(1)
             subplot(2,1,2);
             plot(segmented,'k--')
             hold
             plot(trace,'r-')
             plot([1 length(trace)+1],[f f],'b--')
             set(gca, 'xscale','linear', 'yscale','linear','ylim',[200 350],'xscale','linear','xlim',[1 length(trace)]);
             subplot(2,1,1);
             plot(nolacsegmented,'k--')
             hold
             plot(nolactrace,'r-')
             plot([1 length(nolactrace)+1],[f f],'b--')
             set(gca, 'xscale','linear', 'yscale','linear','ylim',[200 350],'xscale','linear','xlim',[1 length(nolactrace)]);
             leg = legend(['l=',num2str(l),' and f=',num2str(ff)],'Location','NorthOutside');
             set(leg,'Box','on','fontsize', 28);
             %pause(1)
             %clf

             waitforbuttonpress
    %}       

        end
    end    

    ontime_comp=[];
    ontime_nolac_comp=[];
    offtime_comp=[];
    offtime_nolac_comp=[];
    tot_trace=0;
    tot_nolac_trace=0;
    ontime_accept=[];
    ontime_nolac_accept=[];
    offtime_accept=[];
    offtime_nolac_accept=[];
    %make sure there are no transitions in the calibration data
    accepted=zeros(1,length(onstate_nolac_comp));
    accepted_nolac=zeros(1,length(onstate_nolac_comp));
    ontime_accept=[];
    ontime_nolac_accept=[];
    offtime_accept=[];
    offtime_nolac_accept=[];
    temp_on=[];
    temp_nolac_on=[];
    temp_off=[];
    temp_nolac_off=[];
    tmp=[];
    tmp_nolac=[];
    st=[];
    st_nolac=[];
    for l=1:size(onstate_nolac_comp,2)
        ontime_comp{l} = 0;
        for m=1:size(onstate_nolac_comp{l},1)
            %First do lac transitions
            acc=0;
            temp_on=onstate_comp{l}(m,:);
            temp_off=offstate_comp{l}(m,:);  
            tmp=trace_comp{l}(m,:);
            st=statetrace_comp{l}(m,:);
            %This gets rid of data from trajectories that ended early.  Otherwise
            %state continues until the end of the array
        %     temp_on(tmp==0)=0;
        %     temp_off(tmp==0)=0;
        %     st(tmp==0)=NaN;
            temp_on_label=bwlabel(temp_on);
            temp_off_label=bwlabel(temp_off);
            ontime_temp=regionprops(temp_on_label);
            offtime_temp=regionprops(temp_off_label);
            offstate_comp{l}(m,:)=temp_off;
            onstate_comp{l}(m,:)=temp_on;
            statetrace_comp{l}(m,:)=st;

            %make sure there are more than 2 states
        %     if (max(temp_on_label)+max(temp_off_label))<3
        %         acc=1;
        %     end
        %     % make sure it didn't transition to looped in calibration 
        %     if max(onstate_nolac_comp{t})==1
        %         acc=2;
        %     end
        %     % Check if both flags got flipped
        %     if (max(temp_on_label)+max(temp_off_label))<3 && max(onstate_nolac_comp{t})==1
        %         acc=3;
        %     end
            tot_trace=tot_trace+1;
            centroids=cat(2,[ontime_temp.Centroid],[offtime_temp.Centroid]); %compile centroids, takes both x and y component
            centroidsnew=[];
            for kk=2:2:length(centroids)
                centroidsnew(kk/2)=centroids(kk); %remove the extra dimension... why is this so hard.
            end
            centroids=centroidsnew;

            firststate=min(centroids);
            endstate=max(centroids);
            for i=1:length(ontime_temp)
                %if ontime_temp(i).Centroid(2)<endstate && ontime_temp(i).Centroid(2)>firststate
                ontime_comp{l}(m+1,i)=ontime_temp(i).Area;
                if ontime_temp(i).Centroid(2)==firststate
                ontime_accept=cat(1,ontime_accept,1);
                elseif ontime_temp(i).Centroid(2)==endstate
                ontime_accept=cat(1,ontime_accept,1);
                else
                ontime_accept=cat(1,ontime_accept,acc);
                end        
                %end
            end
            for i=1:length(offtime_temp)
                %if offtime_temp(i).Centroid(2)<endstate && offtime_temp(i).Centroid(2)>firststate
                offtime_comp{l}(m,i)=offtime_temp(i).Area;
                if offtime_temp(i).Centroid(2)==firststate
                offtime_accept=cat(1,offtime_accept,1);
                elseif offtime_temp(i).Centroid(2)==endstate
                offtime_accept=cat(1,offtime_accept,1);
                else
                offtime_accept=cat(1,offtime_accept,acc);
                end        
                %end
            end
            %      end
            % end
            %accepted(t)=acc;

            %Now do the whole thing with no lac
            acc_nolac=0;
            temp_nolac_on=onstate_nolac_comp{l}(m,:);
            temp_nolac_off=offstate_nolac_comp{l}(m,:);  
            tmp_nolac=trace_nolac_comp{l}(m,:);
            st_nolac=statetrace_nolac_comp{l}(m,:);
            %This gets rid of data from trajectories that ended early.  Otherwise
            %state continues until the end of the array
        %     temp_on(tmp==0)=0;
        %     temp_off(tmp==0)=0;
        %     st(tmp==0)=NaN;
            temp_nolac_on_label=bwlabel(temp_nolac_on);
            temp_nolac_off_label=bwlabel(temp_nolac_off);
            ontime_nolac_temp=regionprops(temp_nolac_on_label);
            offtime_nolac_temp=regionprops(temp_nolac_off_label);
            offstate_nolac_comp{l}(m,:)=temp_nolac_off;
            onstate_nolac_comp{l}(m,:)=temp_nolac_on;
            statetrace_nolac_comp{l}(m,:)=st_nolac;

            %make sure there are more than 2 states
        %     if (max(temp_on_label)+max(temp_off_label))<3
        %         acc_nolac=1;
        %     end
        %     % make sure it didn't transition to looped in calibration 
        %     if max(onstate_nolac_comp{t})==1
        %         acc_nolac=2;
        %     end
        %     % Check if both flags got flipped
        %     if (max(temp_on_label)+max(temp_off_label))<3 && max(onstate_nolac_comp{t})==1
        %         acc_nolac=3;
        %     end
            tot_nolac_trace=tot_nolac_trace+1;
            centroids_nolac=cat(2,[ontime_nolac_temp.Centroid],[offtime_nolac_temp.Centroid]); %compile centroids, takes both x and y component
            centroidsnew_nolac=[];
            for kk=2:2:length(centroids_nolac)
                centroidsnew_nolac(kk/2)=centroids_nolac(kk); %remove the extra dimension... why is this so hard.
            end
            centroids_nolac=centroidsnew_nolac;

            firststate_nolac=min(centroids_nolac);
            endstate_nolac=max(centroids_nolac);
            for i=1:length(ontime_nolac_temp)
                %if ontime_temp(i).Centroid(2)<endstate && ontime_temp(i).Centroid(2)>firststate
                ontime_nolac_comp{l}(m,i)=ontime_nolac_temp(i).Area;
                if ontime_nolac_temp(i).Centroid(2)==firststate_nolac
                ontime_nolac_accept=cat(1,ontime_nolac_accept,1);
                elseif ontime_nolac_temp(i).Centroid(2)==endstate_nolac
                ontime_nolac_accept=cat(1,ontime_nolac_accept,1);
                else
                ontime_nolac_accept=cat(1,ontime_nolac_accept,acc_nolac);
                end        
                %end
            end
            for i=1:length(offtime_nolac_temp)
                %if offtime_temp(i).Centroid(2)<endstate && offtime_temp(i).Centroid(2)>firststate
                offtime_nolac_comp{l}(m,i)=offtime_nolac_temp(i).Area;
                if offtime_nolac_temp(i).Centroid(2)==firststate_nolac
                offtime_nolac_accept=cat(1,offtime_nolac_accept,1);
                elseif offtime_nolac_temp(i).Centroid(2)==endstate_nolac
                offtime_nolac_accept=cat(1,offtime_nolac_accept,1);
                else
                offtime_nolac_accept=cat(1,offtime_nolac_accept,acc_nolac);
                end        
                %end
            end
            %      end
            % end
            %accepted_nolac(t)=acc_nolac;
        end
        
        num = size(ontime_comp{l},1);

        if num < size(time_comp{l},1)
            size_diff = size(time_comp{l},1) - num;
            ontime_comp{l}(end+1:end+size_diff+1,:) = 0;
            ontime_comp{l} = ontime_comp{l}(2:end,:);
        elseif num > size(time_comp{l},1)
            ontime_comp{l} = ontime_comp{l}(2:end,:);
        end
    end
    
    analysisfile=fullfile(path, specific_path,...
        currexp,[mutant_name,'_analysis_',num2str(threshold),'_analyzed.mat']);
    save(analysisfile,'threshold','ontime_accept','ontime_nolac_accept',...
        'offtime_accept','offtime_nolac_accept','offtime_comp',...
        'offtime_nolac_comp','ontime_comp','ontime_nolac_comp',...
        'offstate_comp','offstate_nolac_comp','onstate_comp',...
        'onstate_nolac_comp','allstates_comp','allstates_nolac_comp',...
        'statetrace_comp','statetrace_nolac_comp','trace_comp',...
        'trace_nolac_comp','tot_trace','tot_nolac_trace','n','time_comp',...
        'alllacnames');        

    end
end