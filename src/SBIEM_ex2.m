% Example #2 for SBIEM:
% An asperity breaks and triggers rupture in a neighboring asperity

% get default parameters, with a few modified values:
% play with a range of segment sizes (L) and resolutions (NX)
% L allows control on the spurious periodicity effect
% NX controls the quality of the solution
% To explore efficiently the behavior of the asperities as a function
% of model parameters it is convenient (shorter computing time) to
% start with a small domain and low resolution.
% It can be useful then to verify selected simulations by running 
% them again with 2*L and/or 2*NX

% high resolution
pars = SBIEM([],'L',300e3,'NX',2048);

% low resolution
%pars = SBIEM([],'L',300e3,'NX',1024);

% high resolution, smaller domain:
%pars = SBIEM([],'L',100e3,'NX',1024);

pars.TMAX = 1.5 * pars.L/pars.CS;
pars.OX_IX = [1:pars.NX];
pars.OX_IDT = 1;
pars.OT_IX = pars.NX/2-[0 40 80 120];

% static and dynamic strength:
taus = pars.SIG0*pars.FRIC.MUs;
taud = pars.SIG0*pars.FRIC.MUd;

pars.TAU0 = zeros(pars.NX,1);

% set initial shear stress: two patches at the center of the fault 
% patch #1 and #2: size, location, stress
HET_L1 = 50e3 ; HET_X1=50e3; HET_A1 = (taus-taud)*0.2; 
HET_L2 = 50e3 ; HET_X2=-50e3 ; HET_A2 = (taus-taud)*0.4;
% set stress inside patches (two Heaviside functions)
pars.TAU0(abs(pars.X-HET_X1)<=HET_L1/2) = taud + HET_A1 ;
pars.TAU0(abs(pars.X-HET_X2)<=HET_L2/2) = taud + HET_A2 ;

% set stress between the two asperities slightly below dynamic strength
pars.TAU0(pars.X>HET_X2+HET_L2/2 & pars.X<HET_X1-HET_L1/2) = taud -0.0*(taus-taud);

% nucleation patch at the northern edge of asperity #1
%LNUC = 2e3; XNUC = HET_X1+HET_L1/2 -LNUC/2 ; SNUC = taud +(taus-taud)*1.01; 
% nucleation patch at the center of asperity #1
LNUC = 2e3; XNUC = HET_X1 ; SNUC = taud +(taus-taud)*1.01; 
pars.TAU0( abs(pars.X-XNUC)<=LNUC/2 ) = SNUC;

% nucleation patch at the northern edge of asperity #2
%LNUC = 2e3; XNUC = HET_X2+HET_L2/2 -LNUC/2 ; SNUC = taud +(taus-taud)*0.90; 
%pars.TAU0( abs(pars.X-XNUC)<=LNUC/2 ) = SNUC;

% set low initial stress far out to stop rupture
LBAR = 75e3 ;	% distance from center to barrier
SBAR = taud -0.4*(taus-taud) ; 	% stress in the barrier
pars.TAU0( abs(pars.X)>LBAR ) = SBAR;

% do the computations
[pars,xdat,tdat] = SBIEM(pars);

% plot the results
disp('Plotting ...')
figure(1)
clf

% space-time plot of slip rate 
subplot(5,5,[1 13])
contourf(xdat.Time,xdat.X/1e3,xdat.SlipRate,10,'LineStyle','none')
colorbar('E') % place colorbar INSIDE the contour plot to avoid axis
              % mis-aligned with the other subplots
title('Slip Rate')
ylabel('X (km)')
% show also, for reference, the S wave front emanating from
% the hypocenter and a front 3 times slower
hold on
plot([0 30]*1e3/pars.CS, -[0 30]+XNUC/1e3, 'y',...
     [0 30]*1e3/pars.CS*3, -[0 30]+XNUC/1e3, 'y--' )
hold off

% moment rate 
subplot(5,5,[16 18])
plot(tdat.Time,tdat.MeanSlipRate*pars.L)
ylabel('Moment rate /m')

% maximum slip rate, timeseries
subplot(5,5,[21 23])
plot(tdat.Time,tdat.MaxSlipRate);
xlabel('Time (s)')
ylabel('Max slip rate')

% slip snapshots
subplot(5,5,[4 14])
plot(xdat.Slip(:,1:100:end), xdat.X/1e3,'b')
xlabel('Slip')

% max slip rate, spatial distribution
subplot(5,5,[5 15])
plot(max(xdat.SlipRate'), xdat.X/1e3)
xlabel('Max slip rate')

figure(2)
% initial stress, strength
subplot(211)
xkm = [-1 1]*pars.L/2/1e3; % [xmin xmax] in km
plot(xdat.X/1e3,xdat.Stress(:,1)/1e6, ...
     xkm,[1 1]*pars.SIG0*pars.FRIC.MUs/1e6, ...
     xkm,[1 1]*pars.SIG0*pars.FRIC.MUd/1e6 )
ylabel('Initial stress (MPa)')
legend('Initial stress','Static strength','Dynamic strength', 'Location','NW')
% rupture speed
subplot(212)
vr = speed(xdat.X,xdat.RuptureTime);
plot(xdat.X/1e3,abs(vr)/1e3, ...
     xkm,[1 1]*pars.CS/1e3, '--',...
     xkm,[1 1]*pars.CS/1e3/3, ':')
ylabel('Rupture speed (km/s)')
xlabel('X (km)')

