%function [newlacANAL_RMS,newlacnames,unloopedlengths]=bdconcatV3(lacANAL_RMS,nolacANAL_RMS,corres,corres2,corres3,alllacnames)
%
%Concatenate data on the same beads from different data sets.  Called by
%masterscriptV3; see that function for input definitions.
%
%Stephanie Johnson 1/11

function [newlacANAL_RMS,newlacnames]=bdconcatV3(lacANAL_RMS,corres2,alllacnames)

newlacnames = cell(length(lacANAL_RMS)/2,1);
newlacANAL_RMS = cell(length(lacANAL_RMS)/2,1);

%Duplicate sets are the second half of the lacANA_RMS array
for i=(1:length(lacANAL_RMS)/2)
    newlacnames{i}=alllacnames{i};
    thisset=lacANAL_RMS{i};
    set2=lacANAL_RMS{length(lacANAL_RMS)/2+i}; %This is the second set on the area in question
    if ~isempty(set2) 
        thiscorres2=corres2{i};
        k=1; %Increments through the corres2 matrix and set2, which should be the same length
        j=1; %Increments through thisset
        newset=[];
        newsetlength=size(thisset,2)+size(set2,2);
        while k<=length(thiscorres2) || j<=size(thisset,1)
            if k<=length(thiscorres2) && thiscorres2(k)==j 
                newrow=[thisset(j,:),set2(k,:)];
                newset=[newset;newrow];
                k=k+1;
                j=j+1;
            elseif k<=length(thiscorres2) && thiscorres2(k)==0%Means data was taken on this bead during the second set but not the first.
                newrow=zeros(1,newsetlength);
                newrow(1,(size(thisset,2)+1):end)=set2(k,:);
                newset=[newset;newrow];
                k=k+1;
            elseif j<=size(thisset,1)
                newrow=zeros(1,newsetlength);
                newrow(1,1:size(thisset,2))=thisset(j,:);
                newset=[newset;newrow];
                j=j+1;
            else
                strcat('i=',int2str(i))
                strcat('j=',int2str(j))
                strcat('k=',int2str(k))
            end
        end
    else
        newset=thisset; %Means lacnames2 for this set was empty
    end
    newlacANAL_RMS{i}=newset;  
end