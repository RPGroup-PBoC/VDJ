%function filt_vect = gaussfiltV3(inputR, fc,fps)
%
%Apply a Gaussian filter with cutoff frequency fc to TPM data.  Input R is
%a beads-by-frames matrix (or column vector to filter just one trace).
%NOTE: THIS SHOULD BE (x^2+y^2)!  That is, to compute the RMS, we filter
%R^2 and then take the square root.  Fps is the frame rate in Hz.
%
%Returns a beads-by-frames matrix of filtered data--this will be R^2 if inputR
%follows our convention of being R^2.
%
%Stephanie Johnson 1/11, based off of code from Lin Han (see gaussfiltV2
%and previous versions listed there)

function filt_vect = gaussfiltV3(inputR, fc,fps)

numframes = size(inputR, 2);
bds = size(inputR, 1);

fft_filt_vect = zeros(numframes, bds); %this will store the product of the data and the filter

%Set up the independent variable for the filter from -numframes/2 to numframes/2
if rem(numframes,2)==0
    freq = [ (0:numframes/2),(-numframes/2+1:1:-1)];
else
    freq = [ (0:numframes/2),((-numframes-1)/2+1:1:-1)];
end

%Fc is given in Hz, but there needs to be a unit conversion for what goes
%in the filter itself, because freq is unitless but the output of the fft
%is not
cf=floor(numframes/(fps*fc^-1));

%Define the filter
G = exp(-0.3466*(freq/cf).^2); %the 0.3466 is chosen to give 3dB
                                    %of attenuation at the cutoff freq

%Fourier transform the data.  If the input is a matrix, then the fft is done
%column-wise, so need to transpose inputR so that each column is a bead.
fftdata = fft(inputR');

%Apply the filter
for i = 1:bds 
    fft_filt_vect(:,i) = fftdata(:,i).*G'; %multiply by the filter
end

%Transform back to time.
filt_vect = transpose(real(ifft(fft_filt_vect)));

           