function ragANAL_RMS = LACtoRAG(lacANAL_RMS)

ragANAL_RMS = [];

for i = 1:length(lacANAL_RMS)
    for j = 1:size(lacANAL_RMS{1,i},1)
        if lacANAL_RMS{1,i}(j, end) > 320 || lacANAL_RMS{1,i}(j, end) == 0
            k = find(lacANAL_RMS{1,i}(j,:),1,'last');
            [~, peaks] = findpeaks(-lacANAL_RMS{1,i}(j,1:k));
            if isempty(peaks)
                ragANAL_RMS{1, end+1}(1,:) = lacANAL_RMS{1,i}(j,1:k);
                continue
            else
                last_peak = peaks(end);
            end
            if last_peak + 400 < length(lacANAL_RMS{1,i}(j,:))
                ragANAL_RMS{1, end+1}(1,:) = lacANAL_RMS{1,i}(j,1:last_peak+400);
            else
                ragANAL_RMS{1, end+1}(1,:) = lacANAL_RMS{1,i}(j,:);
            end
        else
            ragANAL_RMS{1, end+1}(1,:) = lacANAL_RMS{1,i}(j,:);
        end
    end
    
end

sample = [];

for i=1:length(ragANAL_RMS)
    sample{1,i}(1,:) = ragANAL_RMS{i}(1,:);
end

m = ceil(sqrt(length(sample)));
unlooped_hmgb1 = 318;
mid = 285;
looped = 255;
for r = 1:4
    num = [9,10,12,15];
    j = num(r);
    subplot(2,2,r);
    len = length(sample{j});
    x = linspace(0,len/120,len);
    plot(x,sample{j}, 'LineWidth',2)
    
    hold on
    plot([0 len/120],[unlooped_hmgb1 unlooped_hmgb1],'g','LineWidth',2)
    plot([0 len/120],[mid mid],'--k','LineWidth',2)
    plot([0 len/120],[looped looped],'r','LineWidth',2)
    %}
    xlim([0 len/120])
    ylim([210 350])
    if rem(r,2)==1
        ylabel('Bead RMS (nm)', 'FontSize',16)
    end
    if r==3
        xlabel('Time (sec)', 'FontSize',16)
    end
end