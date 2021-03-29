% Example #1 for SBIEM:
% a bilateral rupture, nucleated on an overstressed patch and 
% stopped by a strong barrier (based on the SCEC benchmark #3)

% This example script illustrates the usage of SBIEM
% and the visualization of its outputs.

% get the default "input parameter structure"
pars = SBIEM

% The following sample case is based on the SCEC benchmark #3.
% Actually, most of the input parameters are set by default in SBIEM.m
% to the values of the SCEC benchmark #3, 
% so we only need to do few changes in the input parameter structure:

% set initial shear stress: an overstressed patch at the center of the fault 
% this is where rupture will start
HET_L = 3e3 ;		% patch size [3 km]
HET_A = 11.6e6;		% patch stress [11.6 MPa]
pars.TAU0 = pars.TAU0 + HET_A*(abs(pars.X)<=HET_L/2);

% set initial normal stress : a stopping barrier with high strength
BAR_L = pars.L/4 ;	% distance from center to barrier
BAR_A = 10 ; 		% relative strength of the barrier
pars.SIG0 = pars.SIG0 *(1+ (BAR_A-1)*(abs(pars.X)>=BAR_L) );

% do the computations ...
[pars,xdat,tdat] = SBIEM(pars);

% plot the results ...
figure(gcf)
clf

% mean slip rate (average over the whole fault) as a function of time
subplot(321)
plot(tdat.Time,tdat.MeanSlipRate)
ylabel('Mean slip rate')

% mean slip (average over the whole fault) as a function of time
subplot(323)
plot(tdat.Time,tdat.MeanSlip)
ylabel('Mean slip')

% rupture front plot: rupture time, tail of the process zone and healing time
% as a function of crack front position
% The healing front in this simple symmetric case 
% is obtained from the output 'CrackLength' 
% (plots the same as rupture time before the barriers are reached)
subplot(325)
plot(tdat.CrackLength/2,tdat.Time,xdat.X,xdat.RuptureTime,xdat.X,xdat.ProcessZone);
grid on
xlabel('X')
ylabel('Time')
legend('Half-crack length','Rupture front','Tail of process zone')

% slip rate at selected points
% their location was set by pars.OX_IX 
% their positions are tdat.X
subplot(322)
plot(tdat.Time,tdat.SlipRate)
ylabel('Slip rate')
xlabel('Time')
legend(num2str(tdat.X))

% overview of slip rate over half of the fault, for all times 
% space and time are coarsely sampled for this output
% time sampling was set by pars.OX_IDT and output locations by pars.OX_IX
subplot(3,2,[4 6])
surf(xdat.X,xdat.Time,xdat.SlipRate')
xlabel('X')
ylabel('Time')
zlabel('Slip rate')
shading interp
