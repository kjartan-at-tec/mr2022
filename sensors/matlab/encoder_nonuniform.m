% Script for generating encoder signals sequence
clear all
close all

PPR = 4;
N = 2000;
t = linspace(0, N, N)/N;
wprofile = 2*pi*exp(-((t-0.5)*2).^2);
%w = 2*pi*sin(t/N*pi);

Texp = 20;
tt = linspace(0, Texp, N);

F0 = 0.12;
F1 = 0.18;
yA = chirp(tt,F0,Texp/2, F1, 'quadratic', 0, 'concave');
yB = chirp(tt,F0,Texp/2, F1, 'quadratic', 45, 'concave');

figure(1)
clf
subplot(211)
plot(yA)
subplot(212)
plot(yB)

A = (yA > sind(45)) + (-yA > sind(45));
B = (yB > sind(45)) + (-yB > sind(45));
figure(1)
subplot(211)
hold on
plot(A)
subplot(212)
hold on
plot(B)

% Find rising and falling edges
dA = gradient(A);                % Take Derivative
figure(2)
clf
plot(dA)

[dArpks,iAr] = findpeaks(dA, 'MinPeakDistance',7, 'MinPeakHeight',0.3);
[dAfpks,iAf] = findpeaks(-dA, 'MinPeakDistance',7, 'MinPeakHeight',0.3);
hold on
plot(iAr, dArpks, '^g', 'MarkerFaceColor','g')
plot(iAf, dAfpks, '^r', 'MarkerFaceColor','r')

% Find rising and falling edges
dB = gradient(B);                % Take Derivative
figure(2)
clf
plot(dB)

[dBrpks,iBr] = findpeaks(dB, 'MinPeakDistance',7, 'MinPeakHeight',0.3);
[dBfpks,iBf] = findpeaks(-dB, 'MinPeakDistance',7, 'MinPeakHeight',0.3);
hold on
plot(iBr, dBrpks, '^g', 'MarkerFaceColor','g')
plot(iBf, dBfpks, '^r', 'MarkerFaceColor','r')




Np = min([length(iAr), length(iAf), length(iBr), length(iBf)]);
iArv = cat(2, iAr(1:Np)', ones(Np,1));
iAfv = cat(2, iAf(1:Np)', zeros(Np,1));
iBrv = cat(2, iBr(1:Np)', ones(Np,1));
iBfv = cat(2, iBf(1:Np)', zeros(Np,1));
iAunsorted = cat(1, iArv, iAfv);
[iA, sortind] = sort(iAunsorted(:,1));
iAv = iAunsorted(sortind, :);
iBunsorted = cat(1, iBrv, iBfv);
[iB, sortind] = sort(iBunsorted(:,1));
iBv = iBunsorted(sortind, :);

iAB = sort(cat(1, iAv(:,1), iBv(:,1)));
iABAv = interp1(iAv(:,1), iAv(:,2), iAB, 'previous');
iABBv = interp1(iBv(:,1), iBv(:,2), iAB, 'previous');

counts = 0:(length(iAB)-1);

evtimes = tt(iAB)';
diffs = evtimes(2:end) - evtimes(1:end-1);
diffs = cat(1, diffs(1), diffs); % Repeat first element
freqs = 1./diffs;
data = cat(2, evtimes, iABAv, iABBv, counts', freqs);

dlmwrite('../../figures/chirpdata.dat', data, ',');



%% Velocity and acc

dts = [0.4, 2];
tend = tt(end);
figure(4)
figure(5)
plot(tt(iAB), counts+ 1)

filenames = {'../../figures/chirpdata_dt1.dat', '../../figures/chirpdata_dt2.dat'};
for i = 1:2
    dt = dts(i);
    ti = 0:dt:tend;
    ci = interp1(tt(iAB), counts, ti, 'previous');
    veli = gradient(ci)/dt * (2*pi/32)*1000;
    figure(4)
    subplot(2,1,i)
    plot(ti, veli)
    dlmwrite(filenames{i}, cat(2, ti(:), veli(:)), ',');
end


    

