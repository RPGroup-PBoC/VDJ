function [statetrace,allstates]=segmenttrace(trace,threshold,trackthresh,stickthresh,deadfilter)
offstates_label=[];
onstates_label=[];
trackerror_label=[];
stuck_label=[];
trackerror=[];
onstates=[];
offstates=[];
stuck=[];

spurfirst=0;%1: resolve spurious states first then short states, 0 resolve short states first

onstates=trace>threshold & trace<=trackthresh;
offstates=trace<=threshold & trace> stickthresh;
trackerror=trace>trackthresh;
stuck=trace<=stickthresh;

%each trace labeled by state occurence
trackerror_label=bwlabel(trackerror);
onstates_label=bwlabel(onstates);
offstates_label=bwlabel(offstates);
stuck_label=bwlabel(stuck);

%Hard coded to allow for 100 occurences of a state per trace
trackerror_label(trackerror_label>0)=trackerror_label(trackerror_label>0)+300;
onstates_label(onstates_label>0)=onstates_label(onstates_label>0)+200;
offstates_label(offstates_label>0)=offstates_label(offstates_label>0)+100;

%combine them and relabel
allstates=trackerror_label+onstates_label+offstates_label+stuck_label;
[allstates]=relabeltrace(allstates);
%also need the segmented trace with stick=1, offstate = 2, onstate=3,
%trackerror=4
statetrace=double(stuck)+2*double(offstates)+3*double(onstates)+4*double(trackerror);
%statetrace_short_spur=statetrace;
%properties

if spurfirst==1
    [allstates,statetrace]=removespur(allstates,statetrace);
    allstates_prop=regionprops(allstates);
    %while (sum([allstates_prop.Area]<2*deadfilter))>0
    [allstates,statetrace]=removeshort(allstates,statetrace,deadfilter);
    allstates_prop=regionprops(allstates);
    %end
else
    allstates_prop=regionprops(allstates);
    sss=0;
    while (sum([allstates_prop.Area]<2*deadfilter))>0
        sss=sss+1
    [allstates,statetrace]=removeshort(allstates,statetrace,deadfilter);   
    allstates_prop=regionprops(allstates);
    end
    
    %Let's remove ending spurious states now
    
    
    [allstates,statetrace]=removespur(allstates,statetrace);
end


%remove spurious states

%first remove short states

end
%
function [allstates,statetrace]=removespur(allstates,statetrace)
allstates_prop=regionprops(allstates);
laststate=max([allstates_prop.Centroid]);
centlist=[allstates_prop.Centroid];
laststateindex=(find(centlist==laststate)+1)/2;
%This is an attempt to eliminate data with trailing garbage states because
%often the data at the end is just pure garbage and needs removal
allstates_prop(laststateindex).Centroid(1)
while statetrace(round(allstates_prop(laststateindex).Centroid(1)))==4 || statetrace(round(allstates_prop(laststateindex).Centroid(1)))==1
        allstates_prop(laststateindex).Centroid(1)
        cent = allstates_prop(laststateindex).Centroid(1);
        halfwidth = allstates_prop(laststateindex).Area/2;
        prevstate=statetrace(ceil(cent-halfwidth)-1);
        statetrace=statetrace(1:ceil(cent-halfwidth)-1);
        allstates=allstates(1:ceil(cent-halfwidth)-1);
        [allstates,statetrace]=decomposerecompose(statetrace);
        [allstates]=relabeltrace(allstates);
        allstates_prop=regionprops(allstates);
        
        allstates_prop=regionprops(allstates);
        laststate=max([allstates_prop.Centroid]);
        centlist=[allstates_prop.Centroid];
        laststateindex=(find(centlist==laststate)+1)/2;
end
        
        
    
firststate=10^100;
for q=1:length(centlist)/2
    firststate=min([firststate centlist(2*q-1)]);
end



for ii=1:size(allstates_prop,1)
    %check for state 4 and not an endpoint
    if statetrace(round(allstates_prop(ii).Centroid(1)))==4 && allstates_prop(ii).Centroid(1)<laststate && allstates_prop(ii).Centroid(1)>firststate
        %divide the spurious state evenly among the two next door
        %neighbor states (nextstate and prevstate)
        cent = allstates_prop(ii).Centroid(1);
        halfwidth = allstates_prop(ii).Area/2;
        nextstate=statetrace(floor(cent+halfwidth)+1);
        prevstate=statetrace(ceil(cent-halfwidth)-1);
        statetrace(ceil(cent-halfwidth):ceil(cent-1))=prevstate;
        statetrace(ceil(cent):floor(cent+halfwidth))=nextstate;
        allstates(ceil(cent-halfwidth):ceil(cent-1))=allstates(ceil(cent-halfwidth)-1);
        allstates(ceil(cent):floor(cent+halfwidth))=allstates(floor(cent+halfwidth)+1);
        
        %check for state 1 and not endpoint
    elseif statetrace(round(allstates_prop(ii).Centroid(1)))==1 && allstates_prop(ii).Centroid(1)<laststate && allstates_prop(ii).Centroid(1)>firststate
        %divide the spurious state evenly among the two next door
        %neighbor states (nextstate and prevstate)
        cent = allstates_prop(ii).Centroid(1);
        halfwidth = allstates_prop(ii).Area/2;
        nextstate=statetrace(floor(cent+halfwidth)+1);
        prevstate=statetrace(ceil(cent-halfwidth)-1);
        statetrace(ceil(cent-halfwidth):ceil(cent-1))=prevstate;
        statetrace(ceil(cent):floor(cent+halfwidth))=nextstate;
        allstates(ceil(cent-halfwidth):ceil(cent-1))=allstates(ceil(cent-halfwidth)-1);
        allstates(ceil(cent):floor(cent+halfwidth))=allstates(floor(cent+halfwidth)+1);
        
        %check for last state as 4
    elseif statetrace(round(allstates_prop(ii).Centroid(1)))==4 && allstates_prop(ii).Centroid(1)==laststate
        %if the last state is spurious put the whole thing in the second
        %to last state
        cent = allstates_prop(ii).Centroid(1);
        halfwidth = allstates_prop(ii).Area/2;
        prevstate=statetrace(ceil(cent-halfwidth)-1);
        statetrace(ceil(cent-halfwidth):ceil(cent-1))=prevstate;
        statetrace(ceil(cent):floor(cent+halfwidth))=prevstate;
        allstates(ceil(cent-halfwidth):floor(cent+halfwidth))=allstates(ceil(cent-halfwidth)-1);
        
        %check for last state as 1
    elseif statetrace(round(allstates_prop(ii).Centroid(1)))==1 && allstates_prop(ii).Centroid(1)==laststate
        %if the last state is spurious put the whole thing in the second
        %to last state
        cent = allstates_prop(ii).Centroid(1);
        halfwidth = allstates_prop(ii).Area/2;
        prevstate=statetrace(ceil(cent-halfwidth)-1);
        statetrace(ceil(cent-halfwidth):ceil(cent-1))=prevstate;
        statetrace(ceil(cent):floor(cent+halfwidth))=prevstate;
        allstates(ceil(cent-halfwidth):floor(cent+halfwidth))=allstates(ceil(cent-halfwidth)-1);
        
        %check for first state as 4
    elseif statetrace(round(allstates_prop(ii).Centroid(1)))==4 && allstates_prop(ii).Centroid(1)==firststate
        %if the first state is short put the whole thing in the second
        %state
        cent = allstates_prop(ii).Centroid(1);
        halfwidth = allstates_prop(ii).Area/2;
        nextstate=statetrace(floor(cent+halfwidth)+1);
        statetrace(ceil(cent-halfwidth):ceil(cent-1))=nextstate;
        statetrace(ceil(cent):floor(cent+halfwidth))=nextstate;
        allstates(ceil(cent-halfwidth):floor(cent+halfwidth))=allstates(floor(cent+halfwidth)+1);
        
    elseif statetrace(round(allstates_prop(ii).Centroid(1)))==1 && allstates_prop(ii).Centroid(1)==firststate
        %if the first state is short put the whole thing in the second
        %state
        figure(133);
        plot(allstates);
        cent = allstates_prop(ii).Centroid(1);
        halfwidth = allstates_prop(ii).Area/2;
        nextstate=statetrace(floor(cent+halfwidth)+1);
        cent
        halfwidth
        ceil(cent-halfwidth)
        ceil(cent-1)
        nextstate;
        statetrace(ceil(cent-halfwidth):ceil(cent-1))
        statetrace(ceil(cent-halfwidth):ceil(cent-1))=nextstate;
        statetrace(ceil(cent):floor(cent+halfwidth))=nextstate;
        allstates(ceil(cent-halfwidth):floor(cent+halfwidth))=allstates(floor(cent+halfwidth)+1);
        close(133)
    end
end
%combine connected states (the reason for this is subtle.  If a short state
%is in between the same state nothing (prior to this) ever tells those two now connected
%states that they ar ethe same).
%relabel states
[allstates,statetrace]=decomposerecompose(statetrace);
[allstates]=relabeltrace(allstates);
allstates_prop=regionprops(allstates);
end

function [allstates,statetrace]=removeshort(allstates,statetrace,deadfilter)
allstates_prop=regionprops(allstates);
laststate=max([allstates_prop.Centroid]);
centlist=[allstates_prop.Centroid];
firststate=10^100;
for q=1:length(centlist)/2;
    firststate=min([firststate centlist(2*q-1)]);
end
for i=1:size(allstates_prop,1)
    if allstates_prop(i).Area<2*deadfilter && allstates_prop(i).Centroid(1)<laststate && allstates_prop(i).Centroid(1)>firststate
        %divide the too short state evenly among the two next door
        %neighbor states (nextstate and prevstate)
        cent = allstates_prop(i).Centroid(1);
        halfwidth = allstates_prop(i).Area/2;
        nextstate=statetrace(floor(cent+halfwidth)+1);
        prevstate=statetrace(ceil(cent-halfwidth)-1);
        statetrace(ceil(cent-halfwidth):ceil(cent-1))=prevstate;
        statetrace(ceil(cent):floor(cent+halfwidth))=nextstate;
        allstates(ceil(cent-halfwidth):ceil(cent-1))=allstates(ceil(cent-halfwidth)-1);
        allstates(ceil(cent):floor(cent+halfwidth))=allstates(floor(cent+halfwidth)+1);
    elseif allstates_prop(i).Area<2*deadfilter && allstates_prop(i).Centroid(1)==laststate
        %if the last state is short put the whole thing in the second
        %to last state
        cent = allstates_prop(i).Centroid(1);
        halfwidth = allstates_prop(i).Area/2;
        prevstate=statetrace(ceil(cent-halfwidth)-1);
        statetrace(ceil(cent-halfwidth):ceil(cent-1))=prevstate;
        statetrace(ceil(cent):floor(cent+halfwidth))=prevstate;
        allstates(ceil(cent-halfwidth):ceil(cent-1))=allstates(ceil(cent-halfwidth)-1);
        allstates(ceil(cent):floor(cent+halfwidth))=allstates(ceil(cent-halfwidth)-1);
    elseif allstates_prop(i).Area<2*deadfilter && allstates_prop(i).Centroid(1)==firststate
        %if the first state is short put the whole thing in the second
        %state
        cent = allstates_prop(i).Centroid(1);
        halfwidth = allstates_prop(i).Area/2;
        nextstate=statetrace(floor(cent+halfwidth)+1);
        statetrace(ceil(cent-halfwidth):ceil(cent-1))=nextstate;
        statetrace(ceil(cent):floor(cent+halfwidth))=nextstate;
        allstates(ceil(cent-halfwidth):ceil(cent-1))=allstates(floor(cent+halfwidth)+1);
        allstates(ceil(cent):floor(cent+halfwidth))=allstates(floor(cent+halfwidth)+1);        
    end
end
%combine connected states (the reason for this is subtle.  If a short state
%is in between the same state nothing (prior to this) ever tells those two now connected
%states that they ar ethe same).
%relabel states
[allstates,statetrace]=decomposerecompose(statetrace);
[allstates]=relabeltrace(allstates);
allstates_prop=regionprops(allstates);
end

% decompose and recompose statetraces
function [allstates,statetrace]=decomposerecompose(statetrace)
stuckstate=zeros(1,length(statetrace));
onstate=zeros(1,length(statetrace));
offstate=zeros(1,length(statetrace));
trackstate=zeros(1,length(statetrace));

stuckstate(statetrace==1)=1;
onstate(statetrace==2)=1;
offstate(statetrace==3)=1;
trackstate(statetrace==4)=1;

stuckstate_label=bwlabel(stuckstate);
onstate_label=bwlabel(onstate);
offstate_label=bwlabel(offstate);
trackstate_label=bwlabel(trackstate);

stuckstate_label(stuckstate_label>0)=stuckstate_label(stuckstate_label>0)+300;
onstate_label(onstate_label>0)=onstate_label(onstate_label>0)+200;
offstate_label(offstate_label>0)=offstate_label(offstate_label>0)+100;

allstates=stuckstate_label+onstate_label+offstate_label+trackstate_label;
end

