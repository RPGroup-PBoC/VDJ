%[corrx,corry,corrR]= meandriftcorrfunctionV3(raw_x,raw_y,fc)
%
% Apply a butterworth filter to bead position-vs.-time traces, on a bead-by-bead
% basis.  Inputs are raw x and y data (use loadrawV3 to load the data from POS files),
% the frame rate of the camera (fps, in frames per second), and fc, the cutoff 
%frequency of the Butterworth in Hz.  We usually use fc = 0.05 Hz.  Output
% is a cell array of drift-corrected x and y data, and also sqrt(corr_x^2+corr_y^2).
%
% Basic code based on butterfilt.m by Lin Han; modifications to fit this
% code with masterscriptV* by Kate Craig and Stephanie Johnson.
%
%In contrast to meandriftcorrfunctionV2, this code simply applies the
%Butterworth filter to the raw data and subtracts the results.
%meandriftcorrfunctionV2, on the other hand, applies the filter to the data
%averaged over 4 seconds, then subtracts this filtered average from the raw
%data.  It also computes the cutoff frequency differently.
%
%Stephanie Johnson 9/10  

function [corrx, corry, corrR]= meandriftcorrfunctionV3(raw_x,raw_y,fps,fc)
 
bds = size(raw_x, 1);
numframes = size(raw_x,2);

%Define and apply a first order Butterworth filter
n = 1; % Order of the filter; larger means a sharper transition at the cutoff frequency

if rem(numframes,2)==0 %Gives a frequency axis from numframes/2 to -numframes/2 
   freq=[(0:numframes/2),(-numframes/2+1:1:-1)];
else
   freq=[(0:numframes/2),((-numframes-1)/2+1:1:-1)];
end

cf=floor(numframes/(fps*fc^-1)); %A unit conversion, since the freq axis is unitless but fc is in Hz.
    %Note that this works out because fft returns a vector the same length
    %as the input, in frequencies that go from 0 to fps/2.  (see the fft
    %help file and the plot in one of the examples).
B=1./(1+(freq/cf).^(2*n)); %The Butterworth filter

%Fourier transform the data. If the input is a matrix, then the fft is done
%column-wise, so need to transpose the inputs.
ft_x = fft(raw_x'); %these are now numframes by beads
ft_y = fft(raw_y');

filt_x = zeros(numframes,bds); %so make these numframes by beads
filt_y = zeros(numframes,bds);

for i = 1:bds
    filt_x(:,i) = ft_x(:,i).*B';  %Apply the filter (need to transpose it to match dimensions)
    filt_y(:,i) = ft_y(:,i).*B';
end

filtered_xdata = transpose(real(ifft(filt_x))); %Transform back to time domain.  This contains the drift.
        %Transpose so it's the same dimensions as the original raw x and y data.
filtered_ydata = transpose(real(ifft(filt_y)));
    
corrx = zeros(bds,numframes);
corry = zeros(bds,numframes);
corrR = zeros(bds,numframes);

for j = 1:bds
    corrx(j,:) = raw_x(j,:) - filtered_xdata(j,:); %Subtract the drift
    corry(j,:) = raw_y(j,:) - filtered_ydata(j,:);
    corrR(j,:) = (corrx(j,:).^2 + corry(j,:).^2).^.5; %Compute non-averaged RMS motion
end

%For debugging:
%Plot the raw data and the drift-corrected data, both as time series and as
%x-vs-y blobs
% for r=1:bds
%     figure
%     subplot(3,1,1)
%         plot([1:numframes]./fps, raw_x(r,:),'-b')
%         hold on
%         plot([1:numframes]./fps, corrx(r,:),'-r')
%         legend('Raw data','Drift Corrected')
%         xlabel('Time (seconds)')
%         ylabel('x-coordinate (nm)')
%         ylim([-500 500])
%     subplot(3,1,2)
%         plot([1:numframes]./fps, raw_y(r,:),'-b')
%         hold on
%         plot([1:numframes]./fps, corry(r,:),'-r')
%         legend('Raw data','Drift Corrected')
%         xlabel('Time (seconds)')
%         ylabel('y-coordinate (nm)')
%         ylim([-500 500])
%     subplot(3,1,3)
%         plot([1:numframes]./fps, (raw_x(r,:).^2 + raw_y(r,:).^2).^.5,'-b')
%         hold on
%         plot([1:numframes]./fps, corrR(r,:),'-r')
%         legend('Raw data','Drift Corrected')
%         xlabel('Time (seconds)')
%         ylabel('RMS (nm)')
%         ylim([0 400])
% 
%     figure
%     subplot(1,2,1)
%         plot(raw_x(r,:),raw_y(r,:),'ob')
%         xlabel('x-coordinate (nm)')
%         ylabel('y-coordinate (nm)')
%         legend('Raw data')
%     subplot(1,2,2)
%         plot(corrx(r,:),corry(r,:),'or')
%         xlabel('x-coordinate (nm)')
%         ylabel('y-coordinate (nm)')
%         legend('Drift Corrected')
%     
%     disp(int2str(r))
%         
%     pause
%     close
%     close
% end
