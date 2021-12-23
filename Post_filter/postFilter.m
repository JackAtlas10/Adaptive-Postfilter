clear;
close all;
%% read input audio from path
%clean_sig = audioread('audio_transmit_male.wav');
file_name = "sample_edabk2";
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
alpha = 0.9;
beta = 0.5;
mu = 0.5;
% Gain estimate params
delta = 0.99;
% LPC order - Frame parameter
lpc_ord = 10;
frame = 0.02; % ms second
overlap = 0.5; % overlap 50%

frame_length = floor(frame * fs); % frame in samples
overlap_by_samples = floor(frame_length * overlap);
fprintf('Frame width: %f (ms) = %d samples - Overlap: %.2f %c\n',...
    frame*1e3, frame_length, overlap*100, 37);

%%
n_frames = floor((length(sig) - overlap_by_samples) / ...
    (frame_length - overlap_by_samples));
% output

long_out = zeros(size(sig)); % output after long term filter
short_out = zeros(size(sig)); % output after short_term filter
scaled_out = zeros(size(sig)); % output after rescaled of combine filter
% inner paramters
pitch_freq = pitch(sig, fs, ...
        'Method', 'SRH', ...
        'WindowLength', frame_length, ...
        'OverlapLength', overlap_by_samples);
pitch_i = zeros(n_frames, 1); % pitch indicator
voiced_i = zeros(n_frames, 1); % voiced indicator

g_long = ones(n_frames, 1); % gain rescale factor of long-term filterter
% long term scale factor
temp1 = ones(lpc_ord+1, 1);
temp2 = ones(lpc_ord+1, 1);
for ii = 1:lpc_ord
    temp1(ii+1) = beta^(ii);
    temp2(ii+1) = alpha^(ii);
end
%% Apply filters
for ii=1:n_frames
    if (ii<=2)
        continue;
    end
    prev_st = 1 + (ii-2) * frame_length - (ii-1) * overlap_by_samples;
    st = 1 + (ii-1) * frame_length - (ii-1) * overlap_by_samples; % start of frame
    en = ii * frame_length - (ii-1) * overlap_by_samples; % end of frame
    mid = (st+en+1)/2;
    double_frame = sig(prev_st: en); % double frame without windowing
    frame = sig(st: en);
    w_frame = frame .* hann(frame_length); % windowed frame signal
    % Find pitch
    
    temp_pitch = round(fs / pitch_freq(ii));
    pitch_i(ii) = temp_pitch;
    % voiced indicator
    
    % Case 1: 3-tap pitch predictor
    
    M = temp_pitch;
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
        voiced_i(ii) = 0;
        temp_voiced = 0;
    else
        beta_coeff = phi_matrix\phi_coeff';
        temp_beta_coeff = sum(beta_coeff);
        % Case 2: 1-tap pitch predictor
        %     voice_indicator(ii) = Phi(double_frame, M, M)\Phi(double_frame, 0, M);
        temp_voiced = f(temp_beta_coeff, Uth);
        voiced_i(ii) = temp_voiced;
    end
    gamma = Cz * temp_voiced;
    lambda = Cp * temp_voiced;
    if (temp_voiced ~= 0)
        g_long(ii) = (temp_voiced-lambda) / ...
            (temp_voiced + gamma);
    end
    
    %% Long-term
    
    ord = temp_pitch; % order of long-term = pitch indicator
    long_num = zeros(ord+1, 1);
    long_num(1) = 1;
    long_num(ord+1) = gamma;
    
    long_denum = zeros(ord+1, 1);
    long_denum(1) = 1;
    long_denum(ord+1) = -lambda;
    
    long_frame_temp = filter(g_long(ii)*long_num, long_denum, w_frame);
    
    long_out(st: mid-1) = long_frame_temp(1: frame_length/2) + ...
        long_out(st: mid-1);
    long_out(mid: en) = long_frame_temp(frame_length/2+1 :frame_length);
    %% Short-term
    
    if(sum(abs(w_frame))~=0)
        temp3 = lpc(w_frame, lpc_ord);
    else
        temp3 = zeros(1, lpc_ord+1);
        temp3(1) = 1;
    end
    short_temp1 = temp3' .* temp1;
    short_num = conv(short_temp1, [1 -mu]);
    short_denum = temp3' .* temp2;
    
    combine_num = conv(short_num, g_long(ii)*long_num);
    combine_denum = conv(short_denum, long_denum);
    
    % unscaled output after short term filter
    temp4 = filter(combine_num, combine_denum, w_frame);
    
    short_out(st: mid-1) = temp4(1: frame_length/2) + ...
        short_out(st: mid-1);
    short_out(mid: en) = temp4(frame_length/2+1 :frame_length);
    %% Automatic Gain Scale
    g1 = zeros(frame_length, 1);    % gain estimator of unfiltered signal
    g2 = zeros(frame_length, 1);    % gain estimator of filtered signal
    g = ones(frame_length, 1);
    for ii1 = 2:frame_length
        g1(ii1) = sqrt(delta * g1(ii1-1)^2 + (1-delta)*w_frame(ii1)^2);
        g2(ii1) = sqrt(delta * g2(ii1-1)^2 + (1-delta)*temp4(ii1)^2);
        
        if(g2(ii1)==0)
            g(ii1) = 1;
        else
            g(ii1) = g1(ii1) / g2(ii1);
        end
    end
    
    % rescale temp4 to get output
    temp5 = temp4 .* g;
    
    scaled_out(st: mid-1) = temp5(1: frame_length/2) +...
        scaled_out(st: mid-1);
    scaled_out(mid: en) = temp5(frame_length/2+1 :frame_length);
end
%% Write to wav file
% audiowrite(strcat(file_name,'_long.wav'), long_out, fs);
% audiowrite(strcat(file_name, '_short.wav'), short_out, fs);
audiowrite(strcat(file_name, '_scaled.wav'), scaled_out, fs);

% %%
% figure(1);
% spectrogram(sig,hann(frame_length),overlap_by_samples,1024,'yaxis');
% colormap jet;
% title('Unfiltered Signal');
% %%
% figure(2);
% spectrogram(long_out,hann(frame_length),overlap_by_samples,1024,'yaxis');
% title('Long-term Filtered Signal');
% colormap jet;
% %%
% figure(3);
% spectrogram(short_out,hann(frame_length),overlap_by_samples,1024,'yaxis');
% title('Short-term Filtered Signal');
% colormap jet;
% %%
% figure(4);
% spectrogram(scaled_out,hann(frame_length),overlap_by_samples,1024,'yaxis');
% colormap jet
% title('Combine Filtered Signal');
% %%
% figure(5);
% plot(pitch_freq);
% xlabel('Sample');
% ylabel('Frequency');
% title('Pitch Frequency');
% %%
% figure(6);
% spectrogram(clean_sig,hann(frame_length),overlap_by_samples,1024,'yaxis');
% colormap jet
% title('Clean Signal');