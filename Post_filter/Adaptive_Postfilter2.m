clc;
clear;
%% read input audio from path
file_name = "1cut1";
[sig, fs] = audioread(strcat(file_name,'.wav'));
t = length(sig)/fs; % signal length
fprintf('Signal duration= %f secs\n',t);
fprintf('Sampling frequency= %d Hz\n',fs);
%% Parameters
% Long-term:
Cz = 0.4;
Cp= 0.1;
Uth = 0.6;
% Short-term:
alpha = 0.5;
beta = 0.9;
mu = 0.5;
% Gain estimate params
delta = 0.9;
% LPC order - Frame parameter
lpc_ord = 10;
frame = 0.02; % ms second
overlap = 0.5; % overlap 50%

frame_length = floor(frame * fs); % frame in samples
overlap_by_samples = floor(frame_length * overlap);
fprintf('Frame width: %f (ms) = %d samples - Overlap: %.2f %c\n', frame*1e3, frame_length, overlap*100, 37);

%% find Pitch
%f0 - pitch frequencies
[f0, idx] = pitch(sig, fs, ...
    'Method', 'SRH', ...
    'WindowLength', frame_length, ...
    'OverlapLength', overlap_by_samples ...
    );

%% pitch indicator: pitch period as number of samples
pitch_indicator = round(fs ./ f0);
n_frames = length(f0);
%% ----------------    Long -term postfilter ---------------------
% Voicing Indicator
voice_indicator = zeros(n_frames, 1);
voice_indicator(1) = 0;
voice_indicator(2) = 0;
phi_matrix = zeros(3,3); % phi matrix for solving beta coeffecients
phi_coeff = zeros(3,1);
for ii = 3:n_frames
    
    st = 1 + (ii-2) * frame_length - (ii-1) * overlap_by_samples; % start of frame
    en = ii * frame_length - (ii-1) * overlap_by_samples; % end of frame
    double_frame = sig(st:en);
    % Case 1: 3-tap pitch predictor
    
    M = pitch_indicator(ii);
    % Create Phi matrix for solving Phi_matrix * beta = Phi_coeff
    phi_matrix(1, 1) = Phi(double_frame, M-1, M-1);
    phi_matrix(1, 2) = Phi(double_frame, M-1, M);
    phi_matrix(1, 3) = Phi(double_frame, M-1, M+1);
    
    phi_matrix(2, 1) = Phi(double_frame, M, M-1);
    phi_matrix(2, 2) = Phi(double_frame, M, M);
    phi_matrix(2, 3) = Phi(double_frame, M, M+1);
    
    phi_matrix(3, 1) = Phi(double_frame, M+1, M-1);
    phi_matrix(3, 2) = Phi(double_frame, M+1, M);
    phi_matrix(3, 3) = Phi(double_frame, M+1, M+1);
    %---------------------------------------------- Phi_coeff
    
    phi_coeff(1) = Phi(double_frame, 0, M-1);
    phi_coeff(2) = Phi(double_frame, 0, M);
    phi_coeff(3) = Phi(double_frame, 0, M+1);
    %---------------------------------------------- beta
    if det(phi_matrix) == 0
        voice_indicator(ii) = 0;
        continue;
    end
    beta_coeff = phi_matrix\phi_coeff;
    voice_indicator(ii) = sum(beta_coeff);
    %     Case 2: 1-tap pitch predictor
%          voice_indicator(ii) = Phi(double_frame, M, M)\Phi(double_frame, 0, M);
end

% param in use
gamma = Cz * f(voice_indicator, Uth);
lambda = Cp * f(voice_indicator, Uth);
G_long = ones(size(gamma)); % longterm scale factor
for ii = 1:n_frames
    if voice_indicator(ii) ~= 0
        G_long(ii) = (voice_indicator(ii)-lambda(ii)) / ...
            (voice_indicator(ii) + gamma(ii));
    end
end

%% ------------  Short term coeffecients --------------------
short_tempnum = zeros(n_frames, lpc_ord+1);
short_denum = zeros(n_frames, lpc_ord+1);
temp1 = ones(lpc_ord+1, 1);
temp2 = ones(lpc_ord+1, 1);
for ii = 1:lpc_ord
    temp1(ii+1) = -beta^(ii+1);
    temp2(ii+1) = -alpha^(ii+1);
end
for ii = 1: n_frames
    st = 1 + (ii-1) * (frame_length - overlap_by_samples);
    en = ii * frame_length - (ii-1) * overlap_by_samples;
    mid = (1 + st + en) /2;
    
    sig_frame = sig(st :en) .* hann(frame_length); 
    temp3 = lpc(sig_frame, lpc_ord);
    short_tempnum(ii, :) = temp3' .* temp1;
    short_denum(ii, :) = temp3' .* temp2;
end
%% Long term filter


%% signal goes through long-term filter

long_out = zeros(size(sig));  % signal after long-term filter
short_out = zeros(size(sig)); % signal after 2 filter (unscale)
scaled_out = zeros(size(sig));

for ii = 1:n_frames
    st = 1 + (ii-1) * (frame_length - overlap_by_samples);
    en = ii * frame_length - (ii-1) * overlap_by_samples;
    mid = (1 + st + en) /2;
    sig_frame = sig(st :en) .* hann(frame_length);
    
    %% long-term
    ord = pitch_indicator(ii); % order of long-term
    long_num = zeros(ord+1, 1);
    long_num(1) = 1;
    long_num(ord+1) = gamma(ii);
    
    long_denum = zeros(ord+1, 1);
    long_denum(1) = 1;
    long_denum(ord+1) = -lambda(ii);
    
    long_frame_o = filter(G_long(ii)*long_num, long_denum, sig_frame);
    
    long_out(st: mid-1) = long_frame_o(1: frame_length/2) + ...
        long_out(st: mid-1);
    long_out(mid: en) = long_frame_o(frame_length/2+1 :frame_length);
    
    %% combine filter coeffecients
    % short-term
    short_num = conv(short_tempnum(ii, :), [1 -mu]);
    combine_num = conv(short_num, G_long(ii)*long_num);
    combine_denum = conv(short_denum(ii, :), long_denum);
    % combine
    % unscaled output after short term filter
    temp4 = filter(combine_num, combine_denum, sig_frame);
   
    short_out(st: mid-1) = temp4(1: frame_length/2) + ...
        short_out(st: mid-1);
    short_out(mid: en) = temp4(frame_length/2+1 :frame_length);
    
    %% Automatic Gain Scale
    g1 = zeros(frame_length, 1);    % gain estimator of unfiltered signal
    g2 = zeros(frame_length, 1);    % gain estimator of filtered signal
    g = ones(frame_length, 1);
    for ii2 = 2:frame_length
        g1(ii2) = delta * g1(ii2-1) + (1-delta)*abs(sig_frame(ii2));
        g2(ii2) = delta * g2(ii2-1) + (1-delta)*abs(temp4(ii2));
        g(ii2) = g1(ii2) / g2(ii2);
    end
    
    % rescale temp4 to get output
    temp5 = temp4 .* g;
    if isnan(temp5)
        fprintf('%d nan\n',ii);
    end
    scaled_out(st: mid-1) = temp5(1: frame_length/2) + ...
        scaled_out(st: mid-1);
    if(isnan(scaled_out(st: mid-1)))
        fprintf('%d\n',ii);
    end
    if(isinf(scaled_out(st: mid-1)))
        fprintf('%d is Inf\n',ii);
    end
    scaled_out(mid: en) = temp5(frame_length/2+1 :frame_length);
    if(isnan(scaled_out(mid: en)))
        fprintf('%d is NaN\n',ii);
    end
    if(isinf(scaled_out(mid: en)))
        fprintf('%d is Inf\n',ii);
    end
end
scaled_out(isnan(scaled_out))=0;
% disp(find(isnan(scaled_out)))
%%
% audiowrite(strcat(file_name,'_long.wav'), long_out, fs);
% audiowrite(strcat(file_name, '_short.wav'), short_out, fs);
audiowrite(strcat(file_name, '_scaled.wav'), scaled_out, fs);
%%

figure(1);
spectrogram(sig,hann(frame_length),80,1024,'yaxis');
colormap jet;
title('Unfiltered Signal');

figure(2);
spectrogram(long_out,hann(frame_length),80,1024,'yaxis');
title('Long-term Filtered Signal');
colormap jet;

figure(3);
spectrogram(short_out,hann(frame_length),80,1024,'yaxis');
title('Short-term Filtered Signal');
colormap jet;

figure(4);
spectrogram(scaled_out,hann(frame_length),80,1024,'yaxis');
colormap jet
title('Combine Filtered Signal');

figure(5);
plot(f0);
xlabel('Sample');
ylabel('Frequency');
title('Pitch Frequency');