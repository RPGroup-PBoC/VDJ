%function [unloopedlengths,varagout] = calc_unlooped_lengths(nolacANAL_RMS,corres3)
%
%Called by masterscriptV3.  Returns a vector of unlooped lengths for
%the single bead analysis, and optionally a vector of the c parameter of
%the Gaussian fit.
%
%Stephanie Johnson 1/11

function [unloopedlengths,varargout] = calc_unlooped_lengths(nolacallsets,corres3)

unloopedlengths=[];
unloopedwidths = [];
for i=1:length(corres3)
   thisnolacset=nolacallsets{i};
   
   if size(thisnolacset,1) > size(thisnolacset,2)
       thisnolacset = transpose(thisnolacset); %MasterscriptV2 and V3 save these in different ways :(
   end
   
   thiscorres=corres3{i};
   for b=1:length(thiscorres)
       bd1=thisnolacset(thiscorres(b),:);
       [n, xout] = hist(bd1, 0:1:300);
       prob = n./(sum(n(2:end))); %normalize by the total number of counts
       prob(1) = 0; %Because screenbeads sets 'bad' data to zero, the first bin will contain all these points
       opts = fitoptions('gauss1','Algorithm','Trust-Region');
       try
            unloopedfit = fit(xout(2:end)',prob(2:end)','gauss1',opts);
            unloopedlengths(end+1) = unloopedfit.b1;
            unloopedwidths(end+1) = unloopedfit.c1;
       catch
           unloopedlengths(end+1) = 0;
           unloopedwidths(end+1) = 0;
       end
       
   end 
end
varargout{1} = unloopedwidths;