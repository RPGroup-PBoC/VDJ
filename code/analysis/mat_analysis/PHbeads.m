function PHbeads


figure('Position',[200 200 900 600])
        %Plot the nolac trace
        subplot('Position', [.1,.6,.6,.3])
%        title(strcat(nolactitle,'\_Bead',int2str(b)),'Fontsize',16);
        hold on
        ylim([130 170]) %sets y axis limits for all plots
        xlabel('Time (sec)','Fontsize',16);
        ylabel('<R> (nm)','Fontsize',16)
        set(gca,'FontSize',14)
        plot(x1,bd)
        %plot([0 secs1],[unlop unlop],'k-')%Plotting "guides for the eyes"
        xlim([0 500])
        hold on
        %Plot the nolac histogram
        subplot('Position', [.71,.6,.2,.3])
        [n, xout] = hist(bd, 1);
        prob = n./(sum(n)); 
        prob(1) = 0; 
        barh(xout, prob)
        ylim([130 170])
        set(gca,'YTick',[]);
        set(gca,'XTick',[0.04 0.08 0.12]);
        xlim([0 0.15])
        xlabel('Probability','Fontsize',16)
        set(gca,'FontSize',14)