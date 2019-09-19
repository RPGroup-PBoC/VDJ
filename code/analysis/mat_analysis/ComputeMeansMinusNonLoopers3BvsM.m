%[alldistribssub0,Meanssub0,SEssub0,alldistribssub0B,Meanssub0B,SEssub0B,...
%    alldistribssub0M,Meanssub0M,SEssub0M] = ComputeMeansMinusNonLoopers3BvsM(alldistribs,concrange,alldistribsB,alldistribsM)
%
%Same as ComputeMeansMinusNonLoopers2BvsM but updated 7/2011 to accept
%discrete elements across which to calculate the nonloopers, instead of
%accepting a range.  (This is mainly for the length series.)  Subtracts nonloopers from
%alldistribsB and alldistribsM for cases with two looped states.
%
%Calculates the nonloopers as: instead of subtracting an average number of nonloopers, for
%the elements in concrange (which the user should have chosen such
%that the rest of the distribution is well separated from the zero bin),
%the whole zero bin is excluded from the calculation of the mean.  The
%average fraction of nonloopers is used for those concentrations where the
%distribution includes the zero bin.
%
%Inputs are: one of the alldistribs created by pLoopsvslengthDistrib in
%ComparePLoopCalcMethods; and a concrange vector that gives indices
%across which to calculate the average number of nonloopers.
%
%Steph 3/11

function [alldistribssub0,Meanssub0,SEssub0,alldistribssub0B,Meanssub0B,SEssub0B...
    alldistribssub0M,Meanssub0M,SEssub0M] = ComputeMeansMinusNonLoopers3BvsM(alldistribs,concrange,alldistribsB,alldistribsM)

alldistribssub0=cell(length(alldistribs),1);
alldistribssub0B=cell(length(alldistribs),1);
alldistribssub0M=cell(length(alldistribs),1);

fracnonloops = zeros(length(alldistribs),5); %Only works for 1000 to 5000 seconds

%Find how many nonloopers at each concentration, for each time.  This is
%done only on alldistribs, not alldistribsB/M
for i=1:length(alldistribs)
    tempdistribs=alldistribs{i};
    for j=1:5 %Just doing >1000 sec, >2000 sec ... >5000sec
        thisdistrib = tempdistribs{j}; 
        thisdistrib = sortrows(thisdistrib,2);
        firstnonzero = find(thisdistrib(:,2),1);
        fracnonloops(i,j) = (firstnonzero-1)/size(thisdistrib,1);
        clear thisdistrib
    end
    clear tempdistribs
end

%Find the average number of nonloopers
%Update 7/2011: concrange can be more than 2 elements now
% for c = 1:length(concrange)
%     touse(c,:) = fracnonloops(concrange(c),:);
% end
%faster way:
% touse = fracnonloops(concrange,:);
% avgnonloops = mean(touse,1)
%should be equivalent to:
avgnonloops = mean(fracnonloops(concrange,:),1)
%pause

%For the concentrations in concrange, subtract all nonloopers.  Otherwise
%subtract this average; and calcuate the new mean and SE
Meanssub0 = cell(length(alldistribs),1);
SEssub0 = cell(length(alldistribs),1);
Meanssub0B = cell(length(alldistribs),1);
SEssub0B = cell(length(alldistribs),1);
Meanssub0M = cell(length(alldistribs),1);
SEssub0M = cell(length(alldistribs),1);
for i=1:length(alldistribs)
    tempdistribs=alldistribs{i};
    newdistrib = cell(5,1);  %Again only goes 1000sec ... 5000 sec
    means=zeros(5,1);
    SEs = zeros(5,1);
    tempdistribsB=alldistribsB{i};
    newdistribB = cell(5,1);  %Again only goes 1000sec ... 5000 sec
    meansB=zeros(5,1);
    SEsB = zeros(5,1);
    tempdistribsM=alldistribsM{i};
    newdistribM = cell(5,1);  %Again only goes 1000sec ... 5000 sec
    meansM=zeros(5,1);
    SEsM = zeros(5,1);
    for j=1:5
        thisdistrib = tempdistribs{j};
        thisdistribB = tempdistribsB{j};
        thisdistribM = tempdistribsM{j};
        thisdistrib = sortrows(thisdistrib,2);
        thisdistribB = sortrows(thisdistribB,2);
        thisdistribM = sortrows(thisdistribM,2);
        firstnonzero = find(thisdistrib(:,2),1); %Need to subtract the same number of nonloopers from thisdistrib, thisdistribB, and thisdistribM
        %if i<concrange(1) || i>concrange(2)
        if ~ismember(i,concrange)
            if avgnonloops(j)*size(thisdistrib,1) <= firstnonzero-1 %There are some valid nonloopers
                %There will never be fewer pLoops = 0 for B or M than for the tot
                newdistrib{j} = thisdistrib(floor(avgnonloops(j)*size(thisdistrib,1))+1:end,:);
                newdistribB{j} = thisdistribB(floor(avgnonloops(j)*size(thisdistrib,1))+1:end,:); %and thisdistribB,M have the same length as thisdistrib
                newdistribM{j} = thisdistribM(floor(avgnonloops(j)*size(thisdistrib,1))+1:end,:);
            else %there were fewer than average nonloopers: discard all zero's
                newdistrib{j} = thisdistrib(firstnonzero:end,:);
                newdistribB{j} = thisdistribB(firstnonzero:end,:);
                newdistribM{j} = thisdistribM(firstnonzero:end,:);
            end
        else %discard all nonloopers in concrange
            newdistrib{j} = thisdistrib(firstnonzero:end,:);
            newdistribB{j} = thisdistribB(firstnonzero:end,:);
            newdistribM{j} = thisdistribM(firstnonzero:end,:);
        end
        
        means(j) = mean(newdistrib{j}(:,2));
        SEs(j) = std(newdistrib{j}(:,2))/sqrt(length(newdistrib{j}(:,2))-1);
        meansB(j) = mean(newdistribB{j}(:,2));
        SEsB(j) = std(newdistribB{j}(:,2))/sqrt(length(newdistribB{j}(:,2))-1);
        meansM(j) = mean(newdistribM{j}(:,2));
        SEsM(j) = std(newdistribM{j}(:,2))/sqrt(length(newdistribM{j}(:,2))-1);
    end
    Meanssub0{i} = means;
    SEssub0{i} = SEs;
    alldistribssub0{i} = newdistrib;
    Meanssub0B{i} = meansB;
    SEssub0B{i} = SEsB;
    alldistribssub0B{i} = newdistribB;
    Meanssub0M{i} = meansM;
    SEssub0M{i} = SEsM;
    alldistribssub0M{i} = newdistribM;
    clear means SEs newdistrib meansB SEsB newdistribB meansM SEsM newdistribM
end