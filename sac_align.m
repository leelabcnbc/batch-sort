function sac_align(wave_length)



global WaveformInfo
global FileInfo

Threshold_location = FileInfo(1).ThresholdLocation;

%initialization
ToAlign = WaveformInfo.Waveforms;
[n, d] = size(ToAlign);% n is number of waveforms, d is length

peaks = zeros(1, n);
x = 1:d;

% x1 defines region of interest in the interpolated waveform (*10), where minima may occur 
x1 = ((Threshold_location - 1)*10 + 1):((Threshold_location + 9)*10); 
resamp = zeros(n, d+4);
resamp(:, 3:d+2) = ToAlign;
[r1, r2] = size(resamp);
To_resamp = reshape(resamp', n*r2, 1);
upsamp = reshape(interp(To_resamp, 10, 6, 0.8), r2*10, r1)';%perform 10x interpolation using the function interp (lowpass interpolation)
upsamp = upsamp(:, 20:end-20);% discard edges

%find shifts
[mi, I1] = max((diag(sparse(sign(ToAlign(:,  Threshold_location + 1))))*upsamp(:, x1))'); %find minima

goodpeaks = find((I1>1)&(I1<length(x1)));%reject waveforms that don't have a minimum in the region
ToAlign = ToAlign(goodpeaks, :);
peaks = (Threshold_location - 1) + (I1-1)/10;
peaks(setdiff(1:n, goodpeaks)) =  mean(peaks(goodpeaks));%if peak not in range,  don't shift

if wave_length == 0
    minpeaks = min(peaks);
    maxpeaks = max(peaks);
    
    % 2 peaks align version
    %     array = [-(Threshold_location - 1):d/2-maxpeaks   -(Threshold_location - 1) + d/2:d-maxpeaks ];%data not aligned yet
    %    'data not aligned yet'
    
    array = [-(Threshold_location - 1):d - maxpeaks]; 
    
else
    array = [1:wave_length]- Threshold_location;% length determined by previous alignment
    peaks = min(peaks, d+Threshold_location-wave_length);
end

WaveformInfo.Waveforms = zeros(n, length(array));
for i1 = 1:n
    %     whos n
    WaveformInfo.Waveforms(i1, :) = upsamp(i1, round((peaks(i1)+array)*10));
end   
