clear; close all; clc;
import dsp.*;

Fs = 44100;

figure();

% negative coefficients
a_1_n = [-0.8; -0.6; -0.4; -0.2];

n_labels = strings(size(a_1_n));
for q = 1:size(a_1_n)
    str = strcat('a = ', num2str(a_1_n(q)));
    n_labels(q, 1) = str;
end

subplot(2,1,1, 'align');
xlabel('Frequency (Hz)');
ylabel('Group Delay (samples)');
title("Negative Coefficents");

for p = 1:size(a_1_n)
    g = dsp.AllpassFilter('AllpassCoefficients', a_1_n(p));
    [G, W] = g.grpdelay(8192, Fs);
    hold on
    semilogx(W, G);
end
legend(n_labels);

% positive coefficients
a_1_p = [0.2; 0.4; 0.6; 0.8];
p_labels = strings(size(a_1_p));
for q = 1:size(a_1_p)
    str = strcat('a = ', num2str(a_1_p(q)));
    p_labels(q, 1) = str;
end

subplot(2,1,2, 'align');
xlabel('Frequency (Hz)');
ylabel('Group Delay (samples)');
title("Positive Coefficents")
for p = 1:size(a_1_p)
    g = dsp.AllpassFilter('AllpassCoefficients', a_1_p(p));
    [G, W] = g.grpdelay(8192, Fs);
    hold on
    semilogx(W, G);
end
legend(p_labels);
