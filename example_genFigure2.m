%%%%%%%%%%%%% INIT %%%%%%%%%%%%%

rng(1001)
% Build a homogenous and heterogenous population to show the outputs of the model
NumInputNeurons = 1000; % number of available neurons in the input population
NumStim    = 90; % stimulus space
Tuning = 'Orientation'; % Gaussian, Von Mises, Orientation (Von Mises wrapped at pi)
StimSupport = [-pi/2, pi/2]; % changing this can have bad effects if the tuning curve isn't matched

%%%%%%%%%%%%% construct INPOP and Readouts %%%%%%%%%%%%%

Ihom = InputPopulation('NumNeurons', NumInputNeurons, ...
    'Homogenous',       true, ...
    'NumStim',          NumStim, ...
    'StimSupport',      StimSupport, ...
    'Tuning',           Tuning, ...
    'bandwidthScale', 4, ...
    'maxNoiseCorr', 0.2);

Rhom = ReadoutPopulation(Ihom);

Ihet = InputPopulation('NumNeurons', NumInputNeurons, ...
    'Homogenous',       false, ...
    'NumStim',          NumStim, ...
    'StimSupport',      StimSupport, ...
    'Tuning',           Tuning, ...
    'bandwidthScale', 4, ...
    'maxNoiseCorr', 0.2);

Rhet = ReadoutPopulation(Ihet);

%%%%%%%%%%%%% PLOT RESULTS %%%%%%%%%%%%%

f = figure; clf

%%%% homogenous (simpliest case)
ax1 = subplot(3,2,1);
Ihom.plotTuningSparse
axis square
set(gca,'Xtick',-pi/2:pi/2:pi/2)
set(gca,'Xticklabel',-90:90:90)
set(gca,'TickDir','out')
box off
ylim([0 10])
xlim([(-pi/2)-pi/32 pi/2+pi/32])
xlabel('')
ylabel 'firing rate (Hz)'
title ''

ax3 = subplot(3,2,3);
plotHandles = Rhom.plotWeights(0);
plotHandles.Color = [.2 .2 .9];
set(gca,'Xtick',-pi/2:pi/2:pi/2)
set(gca,'Xticklabel',-90:90:90)
set(gca,'TickDir','out')
xlabel('')
ylim([-.2 .2])
box off
axis square
ylabel 'amplitude'

ax5 = subplot(3,2,5);
SomaPreferredTheta = 0;
plotHandles = Rhom.plotSomaTuning(SomaPreferredTheta);
plotHandles.Color = [.2 .2 .9];
ylim([0 .6])
set(gca,'Ytick',0:.3:.6)
set(gca,'Xtick',-pi/2:pi/2:pi/2)
set(gca,'Xticklabel',-90:90:90)
set(gca,'TickDir','out')
ylabel 'response'
xlabel('orientation (\theta)')
box off
title ''
axis square
%%%%


%%%% hetereogenous 
ax2 = subplot(3,2,2);
Ihet.plotTuningSparse
axis square
set(gca,'Xtick',-pi/2:pi/2:pi/2)
set(gca,'Xticklabel',-90:90:90)
set(gca,'TickDir','out')
box off
ylim([0 10])
xlim([(-pi/2)-pi/32 pi/2+pi/32])
xlabel('')
ylabel 'firing rate (Hz)'
title ''

ax4 = subplot(3,2,4);
plotHandles = Rhet.plotWeights(0);
plotHandles.Color = [.9 .2 .2];
set(gca,'Xtick',-pi/2:pi/2:pi/2)
set(gca,'Xticklabel',-90:90:90)
set(gca,'TickDir','out')
xlabel('')
ylim([-2 2])
title ''
box off
axis square
ylabel 'amplitude'

ax6 = subplot(3,2,6);
SomaPreferredTheta = 0;
plotHandles = Rhet.plotSomaTuning(SomaPreferredTheta);
plotHandles.Color = [.9 .2 .2];
ylim([0 .6])
set(gca,'Ytick',0:.3:.6)
set(gca,'Xtick',-pi/2:pi/2:pi/2)
set(gca,'Xticklabel',-90:90:90)
set(gca,'TickDir','out')
ylabel 'response'
title ''
xlabel('orientation (\theta)')
box off
axis square