%Luke Breuer 12/09, modified to work with masterscriptV3 1/11 by Steph
%
% This is our algorithm:
% 1. start in the middle of the data
% 2. set the "increment" to be 1/4 the range
% 3. test if the max value (intensity of light) is above the threshold
%    a. if it is    , then move to the right by one increment
%    b. if it is not, then move to the left  by one increment
% 4. decrease the increment by half
% 5. while the increment is above 1, goto step 3.
%
% The desired result is that we bounce back and forth on the sides of the
% break, that is, when the maximum value in each frame drops below the
% threshold value.  This will work, as long as the series of max values are
% monotonic.  If they aren't, here are a few possible workarounds:
%   - instead of looking at just one frame, look at several and perform
%     some kind of averaging, or discard the outliers
%   - whenever a value below the threshold is detected:
%     - disallow increasing the current frame past that value and/or
%     - don't decrease the increment value
function [break_frame, x_values, y_values] = find_bead_breakV3(path,filename,beadno,fps)
%     if ismac
%         %path = '/Volumes/drobo/stephj/TPM data';
%         path = '~/Desktop/temp data+code'
%     else
%         path = 'T:\stephj\TPM data';
%     end
% 
%     filename = '091111_OidE85O2+11_slide2_500fMOLDlac_area1';
    
    file_count = length(dir(fullfile(path, strcat(filename, '*.pxl'))));    
    first_imgstemp = PXLtomatrix(path, filename, 1);
    [x_size, y_size, nroi, frames_per_file] = size(first_imgstemp);
    first_imgs=first_imgstemp(:,:,beadno,:);
    clear first_imgstemp

    % data_read is essentially a cache of the files our algorithm decided
    % to read in; get_frame will check this cache to see whether it needs
    % to load a new file (and cache it), or simply return a cached one
    data_read = cell(file_count);
    data_read{1} = first_imgs;
    loaded_files = 0;
    
    function img = get_frame(frame)
        file_num = floor(frame / frames_per_file) + 1; % 0-based to 1-based
        
        if isempty(data_read{file_num})
            data_readtemp = PXLtomatrix(path, filename, file_num);
            data_read{file_num}=data_readtemp(:,:,beadno,:);
            clear data_readtemp
            loaded_files = loaded_files + 1;
        end
        
        img = data_read{file_num}(:, :, 1, mod(frame, frames_per_file) + 1);
    end
    
    % this is our main algorithm
    %threshold = 50; %Steph 8/10 changed this to be 30% of the max of the
    %first frame:
    threshold = 0.3*max(max(data_read{1}(:,:,1,1)));
    current_frame = file_count * frames_per_file / 2;
    increment = file_count * frames_per_file / 4;
    
    while (increment > 1)
        img = get_frame(current_frame);
        value = max(img(:));
        
        if value > threshold
            current_frame = current_frame + increment;
        else
            current_frame = current_frame - increment;
        end
        
        increment = floor(increment / 2);
        
        disp(sprintf('current_frame: %i  increment: %i  value: %i', current_frame, increment, value));
    end
    
    % the search above didn't necessarily give us the exact frame; this 
    % will both do that and force the frames surrounding the break to be
    % loaded, which is useful for eyeballing the analysis
    for i = current_frame - 200 : current_frame + 200
        img = get_frame(i);
        value = max(img(:));
        
        if value < threshold
            break_frame = i;
            disp(sprintf('break frame: %i, which is %4.2f seconds', i, i / fps));
            break;
        end
    end
    
    %********Steph commented out the stuff between asterisks to get it to
    %run faster
%     % load some more frames for better eyeballing
%     for i = 1 : frames_per_file * 100 : frames_per_file * file_count
%         %disp(int2str(i));
%         img = get_frame(i);
%     end
%     
    %we can be lazy and 
    x_values = [];
    y_values = [];
%     
%     function the_max = get_max(matrix_4)
%         the_max = zeros(size(matrix_4, 4), 1);
%         
%         for j = 1 : size(matrix_4, 4)
%             img = matrix_4(:, :, 1, j);
%             the_max(j) = max(img(:));
%         end
%     end
%     
%     for i = 1 : file_count
%         if ~isempty(data_read{i})
%             x_values = [x_values ((i - 1) * frames_per_file + 1):(i * frames_per_file)];
%             y_values = [y_values get_max(data_read{i}(:, :, beadno, :))];
%         end
%     end
% 
%     % for some reason, x_values is 1-d but y_values is 2-d; this fixes that
%     y_values = y_values(:);
%     figure
%     plot(x_values, y_values);
%     pause
%     close
    %***********
    
    %Show the movie for the 4 seconds around the break point
    play_movie(PXLtomatrix(path, filename, ceil(break_frame/frames_per_file)),x_size,y_size,nroi,frames_per_file)
end