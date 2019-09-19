
function [meanie] = SingleBdAnalysisNoProtein(nolacANAL_RMS)
 
num=[];
    i=1:length(nolacANAL_RMS);
    set=nolacANAL_RMS{2};
    for j =[1 2 3 4 6 11 12 13]
       %[1 5 6 7 9 10]  [1 7 8 11 14 15]}
        
         bd = set(j,:);
        [n, xout] = hist(bd, 0:1:400);
        prob = n./(sum(n(2:end))); %normalize by the total number of counts
        prob(1) = 0; %Because screenbeads sets 'bad' data to zero, the first bin will contain all these points
        opts = fitoptions('gauss1','Algorithm','Trust-Region');
        fitresult = fit(xout(2:end)',prob(2:end)','gauss1',opts);
        num(end+1) = fitresult.b1;     
    end
meanie =(num);
%stderr = std(num)/sqrt(length(num)-1);

%{[1 2 15 17]  [3 4 7 14 27]  [18 20 21 28 29]c

%[1 5 6 10 11 14]  [1 6 7 8 9 10 11 13]}