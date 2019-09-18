%[] = createCorres2and3()
%
%Companion function to masterscriptV3 that creates the corres2 and corres3
%matrices.  Corres2 relates the beads numbers for multiple sets on the same area and is used
%primarily for bdconcat; corres3 relates the bead numbers in lacdataconcat
%to those in the saved nolac data file.
%
%Steph 2/11

function [corres2, corres3] = createCorres2and3(lacnames,alllacnames,record,nolacANAL_RMS)

%Corres2: will have empty matrices for any areas on which there was only
%one set; no second sets results in a length (numareas) array with all empty matrices;
%otherwise corres2 is constructed as described in the comments in
%bdconcat.m

corres2{1}=[];

if ~(isequal(length(lacnames),length(alllacnames))) %Means there's at least one second set
    for i=1:length(alllacnames)/2
        if ~strcmpi(alllacnames{i+length(alllacnames)/2},'') %This area had a second set
            %There's got to be a faster way to do this ...
            firstsetrec=record{i};
            secondsetrec=record{i+length(alllacnames)/2};
            tempcorres2=[];
            firstindex=1; %Indexes the first set bead number
            for j=1:length(firstsetrec)
                if firstsetrec(j)==1 && secondsetrec(j)==1 %Data kept on both sets
                    tempcorres2=[tempcorres2 firstindex];
                    firstindex=firstindex+1;
                elseif firstsetrec(j)==0 && secondsetrec(j)==1 %Data kept on second set only
                    tempcorres2=[tempcorres2 0];
                elseif firstsetrec(j)==1 && secondsetrec(j)==0 %Data kept on first set only
                    firstindex=firstindex+1;
                %Do nothing if both have 0's
                end
            end
            corres2{i}=tempcorres2;
        else
                corres2{i}=[];
        end
    end
else
    %Again there's got to be a faster way than this ...
    for i=1:length(alllacnames)
        corres2{i}=[];
    end
end

%Corres3: Length is numbds in the lacconcat set, value of each element is
%index of bd in corresponding nolac matrix.
if ~(isequal(length(lacnames),length(alllacnames))) %Means there's at least one second set
    for i=1:length(alllacnames)/2
        temp=1:size(nolacANAL_RMS{i},2);
        if ~strcmpi(alllacnames{i+length(alllacnames)/2},'') %This area had a second set
            summedrec=record{i}+record{i+length(alllacnames)/2}; %If data were kept during either the first or second sets,
                %this will have a 1 at that index; if data were kept in
                %both sets, a 2; in neither, a 0.
        else
            summedrec=record{i};
        end
        corres3{i}=temp(logical(summedrec));%Logical turns all nonzero elements to 1
        clear summedrec
    end
else
    for i=1:length(alllacnames)
        temp=1:size(nolacANAL_RMS{i},2);
        corres3{i}=temp(logical(record{i}));
    end
end