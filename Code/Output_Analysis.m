close all;
% Load the audio sample
[y1, Fs] = audioread('female_to_male.wav');
[y4, Fs] = audioread('male_to_female.wav');
[y2, Fs] = audioread('female.wav');
y2=y2(:,1);
[y3, Fs] = audioread('male.wav');
y3=y3(:,1);
if length(y3) > length(y2)
    y2 = [y2.', zeros(1, length(y3) - length(y2))];
    y2=y2.';
else
    y3 = [y3.', zeros(1, length(y2) - length(y3))];
    y3=y3.';
end

if length(y1) > length(y2)
    y2 = [y2.', zeros(1, length(y1) - length(y2))];
    y2=y2.';
else
    y1 = [y1.', zeros(1, length(y2) - length(y1))];
    y1=y1.';
end

if length(y1) > length(y4)
    y2 = [y2.', zeros(1, length(y1) - length(y4))];
    y4=y4.';
else
    y1 = [y1.', zeros(1, length(y4) - length(y1))];
    y1=y1.';
end


  
%female_to_male

% Apply FFT
Y = fft(y2);

% Shift FFT
Yshifted2 = fftshift(Y);
Fs
% Plot FFT
figure;
f = linspace(-Fs/2, Fs/2, length(y2));
ff2= abs(Yshifted2)/max(abs(Yshifted2));
plot(f, abs(Yshifted2)/max(abs(Yshifted2)));hold on;
xlim([-1000 1000]);

xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Apply FFT
Y = fft(y1);

% Shift FFT
Yshifted1 = fftshift(Y);

% Plot FFT
f = linspace(-Fs/2, Fs/2, length(y1));
ff1= abs(Yshifted1)/max(abs(Yshifted1));
plot(f, abs(Yshifted1)/max(abs(Yshifted1)));
xlim([-1000 1000]);
title('Converted to male');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Female','Converted to male')


%male_to_female

% Apply FFT
Y = fft(y3);

% Shift FFT
Yshifted3 = fftshift(Y);
Fs
% Plot FFT
figure;
f = linspace(-Fs/3, Fs/3, length(y3));
ff3= abs(Yshifted3)/max(abs(Yshifted3));
plot(f, abs(Yshifted3)/max(abs(Yshifted3)));hold on;
xlim([-1000 1000]);

xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Apply FFT
Y = fft(y4);

% Shift FFT
Yshifted4 = fftshift(Y);

% Plot FFT
f = linspace(-Fs/2, Fs/2, length(y1));
ff1= abs(Yshifted4)/max(abs(Yshifted4));
plot(f, abs(Yshifted4)/max(abs(Yshifted4)));
xlim([-1000 1000]);
title('Converted to Female');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Male','Converted to Female')







