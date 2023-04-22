close all;
% Load the audio sample
[y2, Fs] = audioread('female.wav');
y2=y2(:,1);
[y3, Fs] = audioread('male.wav');
y3=y3(:,1);
%
if length(y3) > length(y2)
    y2 = [y2.', zeros(1, length(y3) - length(y2))];
    y2=y2.';
else
    y3 = [y3.', zeros(1, length(y2) - length(y3))];
    y3=y3.';
end

% Apply FFT
Y = fft(y2);

% Shift FFT
Yshifted2 = fftshift(Y);

% Plot FFT
figure;
f = linspace(-Fs/2, Fs/2, length(y2));
ff2= abs(Yshifted2)/max(abs(Yshifted2));
plot(f, abs(Yshifted2)/max(abs(Yshifted2)));
xlim([-1000 1000]);
title('Real Female');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Load the audio sample

% Apply FFT
Y = fft(y3);

% Shift FFT
Yshifted3 = fftshift(Y);

% Plot FFT
f = linspace(-Fs/2, Fs/2, length(y3));
figure;
ff3= abs(Yshifted3)/max(abs(Yshifted3));
plot(f, abs(Yshifted3)/max(abs(Yshifted3)));
xlim([-1000 1000]);
title('Real Male');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%Prominant Bands

% initialize the third array with zeros
C = zeros(size(ff2));

% loop through the arrays and set the corresponding value in C
for i = 1:length(ff3)
    if ff2(i) > ff3(i)
        C(i) = 1;
    elseif ff2(i) < ff3(i)
        C(i) = -1;
    end
end

figure;
stem(f, C);
ylim([-2 2]);
xlim([110 400]);
title('Prominance Band');
xlabel('Frequency (Hz)');
ylabel('-1->Male 1->Female');

%Female-Male 

figure;
stem(f,ff2-ff3);
xlim([150 400]);
title('Female - Male Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');