function sample_trace_plots(n, corres3, lacANAL_RMS, nolacANAL_RMS)

% n corresponds to the nth trajectory in lacANAL_RMS

% Length of time of no-protein control
t_np = size(nolacANAL_RMS{1}(corres3{1}(n),:),2);
time_np = linspace(0, t_np./30, t_np);

% Length of time of experiment
t_p = size(lacANAL_RMS{1}(n,:),2);
time_p = linspace(0, t_p./30, t_p);

np_rms = 339.5;
unlooped_hmgb1 = 318;
looped = 255;

subplot(2,1,1)
plot(time_np, nolacANAL_RMS{1}(corres3{1}(n),:),'LineWidth',2)
ylim([0 400]);
hold on
plot([0 t_np./30], [np_rms np_rms],'--r','LineWidth',2)
title('No Protein Control','FontSize',18);

ylabel('RMS (nm)','FontSize',18)

subplot(2,1,2)
plot(time_p, lacANAL_RMS{1}(n,:),'LineWidth',2)
ylim([0 400])
hold on
plot([0 t_p./30], [unlooped_hmgb1 unlooped_hmgb1],'--r', 'LineWidth', 2)
plot([0 t_p./30], [looped looped],'--g', 'LineWidth', 2);
title('8pM protoRAG, 80nM HMGB1','FontSize',18);

xlabel('Time (seconds)','FontSize',18)
ylabel('RMS (nm)','FontSize',18)