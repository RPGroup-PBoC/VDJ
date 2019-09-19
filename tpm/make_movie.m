function make_movie

path = '/Volumes/Soichi Backup/Volumes/Seagate Expansion Drive/tpm_data';
directory = '12CodC6A';
file_directory =...
    '170810_12CodC6A_cons23rss_2900bp_ISD1200bp_8nMcMR1cMR2_80nMHMGB1_NoProtein3';

t = 1;
n = 4;

file_number = ceil(t/n);

filename = [file_directory '_' num2str(file_number) '.pxl'];

d = open_pxlV3(fullfile(path,directory,file_directory),filename);

for i=1:120
    imagesc(d{i})
    colormap gray
    hold on
    plot([205 215], [200 200], '-w', 'LineWidth', 3)
    hold off
    text(210, 207, '1 micron', 'HorizontalAlignment', 'center', 'Color', 'w')
    print(strcat('tpm_vid/tpm_',int2str(i)),'-djpeg')
end

end