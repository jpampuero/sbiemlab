% Example #0 for SBIEM:
% a stationary crack (fixed size) with a sudden uniform stress drop

% get default parameters
pars = SBIEM([],'L',8e3,'NX',1024);

pars.FRIC.MUs=0.6;
pars.FRIC.MUd=0.6;
pars.TAU0 = 70e6;
pars.SIG0 = 100e6;
pars.TMAX = 10e3/pars.CS;

% set initial normal stress : a stopping barrier
BAR_L = 1e3 /2 ;	% distance from center to barrier
BAR_A = 100 ; 	% relative strength of the barrier
pars.SIG0 = pars.SIG0 *(1+ (BAR_A-1)*(abs(pars.X)>=BAR_L) );

[cosa,pars.OT_IX] = min(abs(pars.X-0));
pars.OX_IX = find(abs(pars.X)<1e3);

% do the computations
[pars,xdat,tdat] = SBIEM(pars);

% plot the results
plot(tdat.Time,tdat.Stress)
ylabel('Stress')
