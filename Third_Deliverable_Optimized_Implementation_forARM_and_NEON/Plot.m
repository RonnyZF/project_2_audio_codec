clc;
close all;
clear all;

[y] = audioread('sample_1-8kHz.wav');
%[y] = audioread('Test.wav');
[m,fs] = audioread('out.wav');
IFFT = csvread('Decoder.csv');

figure()
hold on;
plot(y);
plot(IFFT);
xlim([0 2500]);
legend('Original','Descomprimido');
hold off;

%sound(y);
%sound(IFFT);