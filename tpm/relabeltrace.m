function [as]=relabeltrace(as)
%relabel the labeled matrix as without  combining states
statenum=1;
        for i=1:max(as)
            if ismember(i,as)
                as(as==i)=statenum;
                statenum=statenum+1;
            end
        end