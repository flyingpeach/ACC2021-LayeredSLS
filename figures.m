% Topology
fig1h = figure(1);
plot_graph(adjMtx, nodeCoords, 'k');

plot_special_vertex(plotNode, nodeCoords, 'r');

x = 10; y = 50; width = 580; height = 500;
set(fig1h, 'Position', [x y width height]); % display only

%% Save topology image
savepath            = '..\figures\';
fn                  = [savepath, 'topology'];
fig1h.PaperUnits    = 'centimeters';
fig1h.PaperPosition = [-0.69 -0.4 10 9]; % print properties
fig1h.PaperSize     = [9 8.5];
print(fn,'-dpdf','-r300');

%% Plot trajectory from specific node
fig2h = figure(2);
tCut  = 22;

green  = [0 0.9 0];
blue   = [0.2 0 1];
red    = [0.8 0 0];

cenlin = red;
cenmpc = blue;
loclayered = green;

subplot(4,1,1); title('Phase','fontsize', 10 ,'interpreter','latex'); hold on;
ylabel('s','fontsize', 9 ,'interpreter','latex');
box on;
p1=plot(timeVals(1:tCut), xsCenLin(phaseIdx, 1:tCut) + thetas(plotNode, 1:tCut), '.-', 'MarkerSize', 6);
p2=plot(timeVals, xsCenMPC(phaseIdx, :) + thetas(plotNode, :), '.-', 'MarkerSize', 6);
p3=plot(timeVals, xsLocLayered(phaseIdx, :) + thetas(plotNode, :), '.-', 'MarkerSize', 6);
p4=plot(timeVals, thetas(plotNode, :), '--');
ylim([-6 1]);

legend({'SatCenLin', 'CenMPC', 'LocLayered'}, ...
        'FontSize', 6);

set(p1,'Color', cenlin);
set(p2,'Color', cenmpc);
set(p3,'Color', loclayered);
set(p4,'Color', [0 0 0]);

subplot(4,1,2); title('Frequency','fontsize', 10 ,'interpreter','latex'); hold on;
ylabel('rad/s','fontsize', 9 ,'interpreter','latex');
box on;
f1=plot(timeVals(1:tCut-1), xsCenLin(freqIdx, 1:tCut-1), '.-', 'MarkerSize', 6);
f2=plot(timeVals, xsCenMPC(freqIdx, :), '.-', 'MarkerSize', 6);
f3=plot(timeVals, xsLocLayered(freqIdx, :), '.-', 'MarkerSize', 6);
ylim([-12 8]);

set(f1,'Color', cenlin);
set(f2,'Color', cenmpc);
set(f3,'Color', loclayered);

subplot(4,1,3); title('Actuation','fontsize', 10 ,'interpreter','latex'); hold on;
ylabel('N $\cdotp$ m','fontsize', 9 ,'interpreter','latex');
box on;
a1=plot(timeVals, usCenLin(plotNode, :), '.-', 'MarkerSize', 6);
a2=plot(timeVals, usCenMPC(plotNode, :), '.-', 'MarkerSize', 6);
a3=plot(timeVals, usLocLayered(plotNode, :), '.-', 'MarkerSize', 6);
ylim([-12 12]);
 
set(a1,'Color', cenlin);
set(a2,'Color', cenmpc);
set(a3,'Color', loclayered);

subplot(4,1,4); title('Disturbance','fontsize', 10 ,'interpreter','latex'); hold on;
ylabel('N $\cdotp$ m','fontsize', 9 ,'interpreter','latex');
box on;
pw=plot(timeVals, wsDist(freqIdx, :), '.-', 'MarkerSize', 6);
xlabel('Time (s)');

set(pw,'Color', [0.3 0.3 0.3]);

x = 10; y = 50; width = 550; height = 590;
set(fig2h, 'Position', [x y width height]); % display only

%% Save time trajectory image
savepath            = '..\figures\';
fn                  = [savepath, 'trajectory'];
fig2h.PaperUnits    = 'centimeters';
fig2h.PaperPosition = [-0.6 -0.5 12 13]; % print properties
fig2h.PaperSize     = [10.6 12];
print(fn,'-dpdf','-r300');