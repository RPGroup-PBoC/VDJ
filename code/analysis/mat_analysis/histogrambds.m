
function [xout, prob] = histogrambds(nolacANAL_RMS)
 
num=[];
bd=[];
for i=1:3
    set=nolacANAL_RMS{1:3};
    for j={[1 2 15 17 21 29],[2 3 4 7 13 14 27 ],[18 20 21 28 29]};
        %[2 3 4 5 7 8 10 11 13 15 16 20]  [1 3 4 5 6 7 8 11 12 14 15];
        %1:size(set,1)
        %[3 9 12 13 16 17 19];[3 4 5 7 13];[1 4 7 9 11 16 19 22 23 24 26];
        bd = [bd, set(j,:)];
        [n, xout] = hist(bd, 50:1:350);
        prob = n./(sum(n(2:end))); %normalize by the total number of counts
        prob(1) = 0; %Because screenbeads sets 'bad' data to zero, the first bin will contain all these points
    end
end
 %Because screenbeads sets 'bad' data to zero, the first bin will contain all these points

bar(xout, prob);


