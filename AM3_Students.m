%%%% Amplitude Modulation Simulation %%%%

clear 
close all
clc


%% Generate sampled 5000Hz Carrier 

% Set carrier frequency
carrier_freq = 5000;

% Set sampling frequency (must be > 2X highest frequency component in
% simulation)
sampling_frequency = 20*carrier_freq;
sampling_period = 1 / sampling_frequency;

% Set duration of simulation
num_samples = 10000;
total_time = sampling_period * num_samples;
time = linspace( 0, total_time, num_samples );

% Create carrier array
carrier = cos( 2 * pi * carrier_freq * time );

%generate other carrier - to demonstrate quadrature multiplexing
alt_carrier = sin( 2 * pi * carrier_freq * time);

% Plot carrier
figure
subplot(211);
plot( time(1:500), carrier(1:500) )
xlabel('t (seconds)')
ylabel('c_1(t) (V 1\Omega)')
title('Carrier (5000Hz)')

%plot carrier that is 90 degrees out of phase
subplot(212);
plot(time(1:500),alt_carrier(1:500));
xlabel('t (seconds)');
ylabel( 'c_2(t) (V 1\Omega)');
title('Carrier (5000Hz) 90\circ out of phase');


%% Generate LPF 

% Can also use fdatool.m
fir_length = 501;
message_freq_max = 500; % Hz
normalized_freq = message_freq_max / (sampling_frequency/2);
h = fir1( fir_length, normalized_freq );

% Plot impulse response
figure
stem( 0 : length(h)-1, h )
xlabel('n');
ylabel('h(n)');
title('Discrete LPF Impulse Response');

% Plot discrete frequency response
H = abs( fftshift( fft( h, min(length(h),1024) ) ) );
figure
plot( linspace(-pi,pi,length(H)), H )
xlabel('rads/sample');
ylabel('|H(\omega)|');
title('Discrete LPF Frequency Response');

% Plot effectice continuous frequency response
figure
plot( linspace(-sampling_frequency/2,sampling_frequency/2,length(H)), H )
xlabel('Hz');
ylabel('|H(f)|');
title('Effective Continuous LPF Frequency Response');

%% Generate message bandlimited to 500 Hz

randn('seed',20);
message = randn( 1, num_samples );
message = conv2( message, h, 'same' );

% Optional DSB-LC
% message = message + .25;  

% % Sinusoidal message
% message = cos( 2 * pi * 100 * time );
% message = message + 1.2;


% plot message
figure
plot( time, message )
xlabel('t (seconds)')
ylabel('f(t) (V 1\Omega)')
title('Bandlimited Message Signal (< 500Hz)')

figure
periodogram( message, [], 8192, sampling_frequency );
title('Power Spectral Density of Bandlimited Message ');

%% Perform AM Modulation (DSB-SC) 

modulated_carrier = message .* carrier;

% display the modulated carrier
figure
plot( time, modulated_carrier, 'b', time, message, 'r-' )
xlabel('t (seconds)')
ylabel('x(t) (V 1\Omega)');
axis tight
legend('Modulated Carrier','Original Message')
title('DSB-SC Modulation')


% display PSD of modulated carrier
figure
periodogram( modulated_carrier, [], 8192, sampling_frequency );
title('Power Spectral Density of Modulated Carrier');


%% Perform coherent demodulation 

received_and_mixed = modulated_carrier .* carrier;
demodulated_signal = conv2( received_and_mixed, 2*h, 'same' );

% display the demodulated signal and compare with message 
figure
plot( time, demodulated_signal, 'b', time, message, 'r-' )
xlabel('t (seconds)')
ylabel('z(t) (V 1\Omega)');
axis tight
legend('Demodulated Signal','Original Message')
title('Coherent Demodulation')


