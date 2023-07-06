
function [meanie, stderr] = SingleBdAnalysisNoRSS(corres3, nolacANAL_RMS)
    
num=[];
data = [];

for i = 1:length(corres3)
    for j = 1:length(corres3{1,i})
        data{1,end+1}(1,:) = nolacANAL_RMS{1,i}(corres3{1,i}(j),:);
    end
end

for i=1:length(data);
    set=data{i};
    for j=1:size(set,1)
        %[3 7 20 23]  [1 9 27 32]  [13 17]
         bd = set(j,:);
        [n, xout] = hist(bd, 0:1:400);
        prob = n./(sum(n(2:end))); %normalize by the total number of counts
        prob(1) = 0; %Because screenbeads sets 'bad' data to zero, the first bin will contain all these points
        opts = fitoptions('gauss1','Algorithm','Trust-Region');
        fitresult = fit(xout(2:end)',prob(2:end)','gauss1',opts);
        num(end+1) = fitresult.b1;     
    end
end
meanie =mean(num);
stderr = std(num)/sqrt(length(num)-1);