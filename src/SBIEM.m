% SBIEM a Spectral Boundary Integral Equation Method for 2D mode III rupture dynamics
%
% SYNTAX
%
%  Typical usages:
%
%    par = SBIEM()  
%    get the default input parameters structure
%
%    par = SBIEM(par,'Prop',Value,...)
%    modify some values in an existing input parameters structure
%
%    [par,ox,ot] = SBIEM(par)
%    solve a problem described by an existing input parameters structure
%
%  Less frequent but valid usages:
%
%    par = SBIEM([],'Property1',Value1,'Property2',Value2)  
%    modify some values in the default input parameters structure
%
%    [par,ox,ot] = SBIEM(par,'Prop',Value,...)
%    solve a problem described by an existing input parameters structure
%    with some modified values
%  
%
% INPUTS 	par	Input parameters structure containing the following fields:
%
%        		Physical properties:
%	                 L	fault length. Periodic boundary conditions are
%				assumed beyond the fault limits, so L must be 
%				large enough to avoid contamination by periodicity
%			 CS	S wave velocity
%			 RHO	density
%			 FRIC	friction law structure, each component can be a vector 
%		 		(see "help friction" for details)
%			 HEAL	friction healing (restrengthening):
%				0 = no healing
%				1 = immediate healing upon arrest
%
%        		Background stress:
%			 TAU0	initial shear stress
%			 SIG0	initial normal stress (positive compressive)
%			 TAUR	tectonic shear stressing rate
%
%			Numerical settings:
%			 NX	number of computational grid nodes
%			 TMAX	duration (final time)
%			 RTCUT	truncation time for mode 1 of the 
%				elastodynamic stress transfer kernel
%				normalized by L/CS (typically =2)
%			 QW	ratio of truncation times of Nyquist mode
%				to mode 1, divided by NX/2
%				(typically =4, uniform truncation if =NX/2)
%			 CFL	stability number = timestep*CS/DX
%				(typically =0.5)
%			 REFT	grid refinement ratio (expert)
%
%			Output settings:
%			 OX_FLAG (0 or 1) snapshot output at selected times
%				usually large outputs
%			 OX_IDT	time stride for snapshots
%			 OX_IX	node indices for snapshots (useful for
%				setting outputs only on a specified area
%				or for setting a coarse space sampling)
%			 OT_FLAG (0 or 1) output at every timestep
%				of quantities at selected nodes 
%				and of macroscopic source quantities
%			 OT_IX	node indices for timeseries outputs
%			 ECHO	output information level (2 = verbose)
%
%			Default values are assigned to any empty field in 'par',
%			or to the whole structure if it is empty.
%
%		'Prop' 	the name of a field in the input parameter structure to be reset
%		Value   new value of the above field, overrides defaults and input 'par'
%
% OUTPUTS 	par	structure containing input parameters, as above,
%			and the following additional fields
%			 MU	shear modulus (=RHO*CS^2)
%			 DX	grid spacing (=L/NX)
%			 DT	timestep (=CFL*DX/CS)
%			 NT	total number of timsteps (=ceil(TMAX/DT))
%			 LC	nucleation size in linear slip weakening
%				(Uenishi and Rice, 2003)
%			 SM	inverse nucleation timescale
%				(Ampuero et al, 2002)
%			These fields are reset internally each time SBIEM is called, 
%			and are available for convenient manipulation of 
%			output data or for resetting other input parameters, 
%			so it is useless to modify them externally
%
%		ox	snapshot output structure containing the following fields:
%
%			Basic information:
%			 X 	(size nx = size OX_IX)
%			 Time	(size nt = floor(NT/OX_IDT) )
%
%			Snapshots at selected times, size =[nx,nt]:
%			 Slip
%			 SlipRate
%			 Stress
%			 Strength
%
%			Rupture front information, size =nx
%			 RuptureTime (=inf at unruptured points)
%			 ProcessZone (time at which slip = Dc)
%
%		ot	timeseries output structure containing the following fields:
%
%			Basic information:
%			 X 	(size nx = size OT_IX)
%			 Time 	(size = NT)
%
%			Timeseries at selected fault nodes, size =[NT,nx]:
%			 Slip
%			 SlipRate
%			 Stress
%
%			Timeseries of macroscopic source quantities, size =[NT,1]:
%			 MeanSlip
%			 MaxSlip
%			 MeanSlipRate
%			 MaxSlipRate
%			 CrackLength (size of actively slipping area)
%			 MeanStress
%
%
% EXAMPLES	See SBIEM_ex*.m 
%
% REQUIRES	friction.m
%
% REFERENCES	J. W. Morrisey and P. H. Geubelle (1997),
%		"A numerical scheme for mode III dynamic fracture problems"
%		Int. J. Num. Meth. Eng., 40 (7), 1181-1196.
%
%		N. Lapusta et al. (2000), "Elastodynamic analysis for slow 
%		tectonic loading with spontaneous rupture episodes on faults 
%		with rate- and state-dependent friction"
%		J. Geoph. Res. 105 (B10), 23765-23789.
%
% AUTHOR	Jean-Paul Ampuero	ampuero@erdw.ethz.ch 	
% LAST MODIF	January 2008
% 

function [pars,oX,oT] = SBIEM(parsin,varargin)

%-----------------------------------------------
% parse inputs
%
% WARNING: no checks yet 
% Field names in parsin and Property's should match field names in pars
% or else the default value will be taken.

pars = SBIEM_get_defaults();

% 1. Override defaults by input structure parsin
if nargin && ~isempty(parsin)
  if ~isstruct(parsin)
    error('SBIEM:parsinNotaStructure','First argument must be a parameter structure')
  end
  fparsin = fieldnames(parsin);
  for k=1:length(fparsin),
    pars.(fparsin{k}) = parsin.(fparsin{k}) ;
  end
end

% 2. Override defaults and previously set parameters 
%    with the 'Property',Value pairs in VARARGIN
np2 = length(varargin);
if np2
  if np2/2 ~= floor(np2/2), error('Property/Value arguments must be in pairs'), end
  for k=1:2:np2,
    pars.(upper(varargin{k})) = varargin{k+1};
  end
end

%-----------------------------------------------
% INITIALIZE derived parameters

pars.MU = pars.RHO*pars.CS^2; 	% shear modulus
pars.DX = pars.L/pars.NX ;
pars.X = (-pars.NX/2+1:pars.NX/2).' *pars.DX ;
pars.DT = pars.CFL*pars.DX/pars.CS ; 
pars.NT = ceil(pars.TMAX/pars.DT);

% unfold the parameter structure fields into variables with the same name
fpars = fieldnames(pars);
for k=1:length(fpars),
  eval([fpars{k} '=pars.' fpars{k} ';']);
end

% some indicative quantities
impedance = MU/(2*CS) ;
W = SIG0.*(FRIC.MUs-FRIC.MUd)./ FRIC.Dc; % slip weakening rate
%if length(W)>1, W=W(NX/2); end % value at center
W = max(W);
LC = 2*0.57888694*MU/W ;	% nucleation size
SM = W / impedance ; 		% growth rate
pars.LC = LC;
pars.SM = SM;

% exit here without computations, when called as: pars = SBIEM(...)
if nargout<=1, return, end

%-- print some information abot numerical resolution
if ECHO>1 
  disp(sprintf( 'Nucleation length Lc = %0.4g',LC));
  disp(sprintf( 'Space resolution Lc/dx = %0.4g',LC/DX));
  disp(sprintf( 'Nucleation time scale Tnuc = %0.4g',1/SM) );
  disp(sprintf( 'Time resolution Tnuc/dt = %0.4g',1/(SM*DT)) );
  if pars.FRIC.kind~=1, disp('NOTE: the parameters above are designed for linear slip-weakening'); end
end

%-----------------------------------------------
% Initialize kernels
%
% Convolution integrated with semi-open second order rule
% (Numerical Recipes, combined eqs. 4.1.7 and 4.1.11)
% Use of a semi-open rule simplifies the second pass of the 
% time evolution scheme.
%
% Mode 0 has null contribution, don't need to store it,
% only strictly positive wavenumbers.
% For negative wavenumbers we will use symetries of real FFT
%
% Also implemented: mode-dependent KERNEL CUT-OFF (TCUT) as in 
% Lapusta et al (2000)

NK=NX/2; % number of strictly positive wavenumbers
NKnyq = NX/2+1 ; % Nyquist wavenumber position in classical FFT storage

ik = (1:NK)';
TCUT1 = floor(RTCUT*NX/CFL);
TCUT = ceil(  TCUT1*( 1+(QW-1)*(ik-1)/(NK-1) ) ./ik ); 

pathstr = fileparts(mfilename('fullpath'));
kernel_dir = 'kernels';
kernel_name = fullfile(pathstr,kernel_dir);
if exist(kernel_name,'file') ~= 7, mkdir(pathstr,kernel_dir); end
kernel_name = fullfile(kernel_name, sprintf('kernelV_%u',NX) );

Compute_Kernel = 1;
if exist([kernel_name '.mat'],'file')
  if ECHO > 1, disp('Reading the kernel'), end
  load(kernel_name);
 % --> kernel, KER_NX, KER_RTCUT, KER_QW, KER_CFL
  Compute_Kernel = any([KER_NX,KER_RTCUT,KER_QW,KER_CFL] ~= [NX,RTCUT,QW,CFL]);
end

if Compute_Kernel
%
% Elastodynamic kernel:
%   Kernel(k,t) = 0.5*MU*k * W_III( k*CS*t ) 
%
% Static kernel:
%   StKernel(k) = -0.5*MU*k
%
% Here we compute and save discretized, non-dimensional versions of the kernels :
%
%   kernel(ik).f(it) 	= Kernel(k,t)*DT/impedance 	   with k=2*pi*ik/L and t=it*DT
%			= k*CS*DT * W_III(k*CS*t)
%                       = CFL*knum * W_III( CFL*knum*it )  where knum=2*pi*ik/NX
%
%   stkernel(ik) = StKernel(k)/(MU/DX)
%		 = -0.5*knum
%
% Quadrature weights for semi-open rule (Numerical Recipes eq. 4.1.20)
%   it=1, t=0  : ff1 = 0
%   it=2, t=DT : ff1 = K_1*V_1
%   it=3       : ff1 = 3/2*K_1*V_2 + 1/2*K_2*V_1
%   it=4       : ff1 = 3/2*K_1*V_3 + K_2*V_2 + 1/2*K_3*V_1
%   it=n       : ff1 = 3/2*K_1*V_n + K_2*V_n-1 + ... + 1/2*K_n*V_1
%
  if ECHO > 1, disp('Computing the kernel'), end
  knum = 2*pi*(1:NK)'/NX;
  for ik=1:NK,
    kernel(ik).f = CFL*knum(ik)* W3( CFL*knum(ik)*(1:TCUT(ik)) ) ;
    kernel(ik).f(1) = 3/2*kernel(ik).f(1) ; % = 3/2*K_1 (see quadrature)
  end
  stkernel = - 0.5*knum ; % static kernel
  KER_NX = NX;
  KER_RTCUT = RTCUT;
  KER_QW = QW;
  KER_CFL = CFL;
  save(kernel_name,'kernel','stkernel','KER_NX','KER_RTCUT','KER_QW','KER_CFL')
end

% fix dimensions
for ik=1:NK,
  kernel(ik).f = impedance*kernel(ik).f;
end
stkernel = MU/DX*stkernel ;


%-----------------------------------------------
%-- initialize fields 
dload1 = TAUR*DT ;
load1 = TAU0 ;
stress1 = zeros(NX,1);
stress2 = zeros(NX,1);
stress2 = SIG0;
u = zeros(NX,1);
theta = zeros(NX,1);
v = zeros(NX,1);
for ik=1:NK,
  fv(ik).f = zeros(TCUT(ik),1);
end  
ff1 = zeros(NX,1);
ff2 = zeros(NX,1);

%-- intialize counters
ito = 1;
techo = ceil(NT/10);
ite = 1;
itcyc = ones(NK,1);

%-- initialize outputs
if OX_FLAG
  oX = new_ox(OX_IX, OX_IDT, X,NT,DT);
  if length(FRIC.Dc)==1
    dcpz = FRIC.Dc;
  else
    dcpz = FRIC.Dc(OX_IX);
  end
else
  oX = [];
end
if OT_FLAG
  oT = new_ot(OT_IX, X,NT,DT);
else
  oT = [];
end


%-----------------------------------------------
% BEGIN TIME LOOP

% ABOUT THE TIME EVOLUTION SCHEME:
% Equations are of the form
%
%   impedance*v = -Ffriction(u) + Fstatic(u) + Fdynamic(v) 
%
% where, owing to the use of a semi-open integration rule,
% Fdynamic(v) depends only on previous values of v.
% A model equation is:
%   v = f[u]
% Discrete version:
%   v(n+1) = f[u(n+1)]
%   u(n+1) = u(n) + 0.5*DT*[ v(n)+v(n+1) ]
% Implemented as a predictor-corrector explicit scheme (SECOND ORDER, Heun's method):
%   1. u_1(n+1) = u(n) + DT*v(n)
%   2. v_1(n+1) = f[u_1(n+1)]
%   3. u_2(n+1) = u(n) + 0.5*DT*[ v(n)+v_1(n+1) ]
%               = u_1(n+1) +0.5*DT*[ v_1(n+1)-v(n) ]
%   4. v_2(n+1) = f[u_2(n+1)]

if ECHO > 1, disp('Begin time loop'), end

for it = 1:NT ,

 if it==1
  % initial conditions, t=0
  % to guarantee second order when initial stress drop is abrupt, 
  % the initial velocity must be set here, analytically
  % impedance*v(0) = max( initial_stress - initial_strength, 0)
  stress1(:) = load1;
  strength = stress2.*friction(theta,FRIC);
  v = max( stress1-strength, 0 )/impedance;

 else

  u_old = u;
  v_old = v;
  theta_old = theta;
  
 % only needed for evaluation of rupture time:
  str_old = strength(OX_IX); 
  s1s_old = stress1(OX_IX); 

 %Update external load  
  load1 = load1 + dload1 ;

  for ipass=1:2,

   %Update slip
    du =  0.5*DT*(v+v_old); 
    u = u_old + du;
   
   %Update static stress term
    fftu = fft(u,NX); 
    ff2(1) = 0 ;
    ff2(2:NK+1) = stkernel.*fftu(2:NK+1) ;
    ff2(NKnyq+1:NX) = conj( ff2(NK:-1:2) );
    
   %Trial state: assuming no further slip
   % = load(n+1) +Fstatic(u) +Fdynamic(n+1)
    stress1 = real(ifft(ff1+ff2,NX)) ;
    stress1 = stress1 + load1 ;
  
   %Update strength Ffriction(u)
    theta = theta_old + du;
    frimu = friction(theta,FRIC);
    strength = stress2.*frimu ; % normal stress >0 is compressive
    
   %Solve friction -> update slip rate
    v = max( stress1-strength ,0)/impedance ;
    s1s = stress1(OX_IX); % stress_stick, only needed for evaluation of rupture time
    stress1 = min(stress1,strength) ;

  end
  if HEAL==1, theta(v==0) = 0; end % instantaneous healing upon slip arrest
  
 end

% Update slip rate spectrum and stress functionals for next timestep
% This loop is the most time-consuming section of the code
  fftv = fft(v,NX);
  for ik=1:NK , 

   % Update slip rate spectrum for NEXT time step
   % "fv" stores the slip rate spectrum history: fv(MODE).f(TIME)
   % The storage is CYCLIC in the time dimension:
   % don't need to keep slip history longer than TCUT(mode).
   % "itcyc" points, on entry, to the position of the current time step 
   %         in "fv(mode).f(:)"   = mod(it,TCUT(mode))+1
    itc = itcyc(ik);
    fv(ik).f(itc) = fftv(ik+1);

   % itclast => earliest velocity stored
    tcut = TCUT(ik);
    itclast = itc+1;
    if itclast==tcut, itclast=1; end

   % Update stress functional for NEXT time step
   % Some fixes are done here to match the quadrature weights
   % of the semi-open rule (Numerical Recipes eq. 4.1.20) :
   % it=1, t=0  : ff1 = 0
   % it=2, t=DT : ff1 = K_1*V_1
   % it=3       : ff1 = 3/2*K_1*V_2 + 1/2*K_2*V_1
   % it=4       : ff1 = 3/2*K_1*V_3 + K_2*V_2 + 1/2*K_3*V_1
   % it=n       : ff1 = 3/2*K_1*V_n + K_2*V_n-1 + ... + 1/2*K_n*V_1
    if it==1
      ff1(ik+1) = 2/3*kernel(ik).f(1)* fv(ik).f(1); %fix: kernel has 3/2*K_1
      fv(ik).f(1) = 0.5 *fv(ik).f(1); % fix: 1/2 in last term of quadrature
    elseif it <= tcut
      ff1(ik+1) = sum( kernel(ik).f(1:it) .* fv(ik).f(it:-1:1) );
    else
     % "ittab" points to the slip history from it+1 back to it-TCUT+2
      fv(ik).f(itclast) = 0.5*fv(ik).f(itclast); % fix: 1/2 in last term
      ittab = [ (itc:-1:1) (tcut:-1:itc+1) ];
      ff1(ik+1) = sum( kernel(ik).f .* fv(ik).f(ittab) );
    end

   % Update "itcyc": point to next time step
    itcyc(ik) = itclast;

  end

  ff1(1) = 0 ; % mode 0 is null
  ff1(NKnyq+1:NX) = conj( ff1(NK:-1:2) ); % real signal: complete the spectrum


%--------- outputs :
% Some progress info  
  if ECHO > 1 & (ite == techo | it==NT)
    ite = 1;
    disp(sprintf('Step %i/%i, Vmax = %0.4g, Umax = %0.4g' ...
                ,it,NT,max(v),max(u)))
  else
    ite = ite+1;
  end  

% STORE RESULTS
 % (x,t) fields, coarsely sampled
  if OX_FLAG 
    if it == ito*OX_IDT
      oX.Slip(:,ito)	= u(OX_IX) ;
      oX.SlipRate(:,ito)= v(OX_IX) ;
      oX.Stress(:,ito)  = stress1(OX_IX) ;
      oX.Strength(:,ito)= strength(OX_IX) ;
      ito = ito+1;
    end  
    if it>1 % first rupture time and tail of process zone
     % rupture time based on a velocity threshold (arbitrary definition)
%      oX.RuptureTime = zeroCrossingTime(oX.RuptureTime, ...
%                         v(OX_IX)-VDYN,v_old(OX_IX)-VDYN, it,DT);
     % rupture time based on stress threshold (natural definition)
     % = first zero crossing time of stress_stick and static_strength
      oX.RuptureTime = zeroCrossingTime(oX.RuptureTime,s1s-str_old,s1s_old-str_old, it,DT);
      oX.ProcessZone = zeroCrossingTime(oX.ProcessZone,u(OX_IX)-dcpz,u_old(OX_IX)-dcpz,it,DT);
    end  
  end  
 % full time sampled outputs
  if OT_FLAG
   % fields(t) at selected locations
    oT.SlipRate(it,:)	= v(OT_IX)';
    oT.Slip(it,:)	= u(OT_IX)';
    oT.Stress(it,:)	= stress1(OT_IX)';
   % macroscopic fields
    oT.MeanSlip(it)     = mean(u) ;
    oT.MaxSlip(it)      = max(u) ;
    oT.MeanSlipRate(it) = mean(v) ;
    oT.MaxSlipRate(it)  = max(v) ;
    oT.CrackLength(it)  = nnz(v>0)*DX; 
    oT.MeanStress(it)   = mean(stress1) ;
  end

end 
%----- END TIME LOOP


%====================================================================
%
% Default parameters
%
function pars = SBIEM_get_defaults()

% units
MPa = 1e6;
year = 356*24*60*60 ;

%-- physical properties --
pars.L = 20e3 ; 	% fault length (simulated period)
pars.CS = 3464.; 	% S wave velocity
pars.RHO = 2670.;	% density

%-- background stress --
pars.TAU0 = 70*MPa ;	% initial shear stress 
pars.SIG0 = 120 *MPa ;	% normal stress (+ compressive)
pars.TAUR = 0; %5e6*MPa/year; % tectonic shear loading rate

%-- numerical settings --
pars.REFT = 1; 			%[integer] grid refinement ratio.
pars.NX = 256*pars.REFT; 	% number of space elements, grid size [128]
pars.TMAX = pars.L/pars.CS;	% max time

%-- friction --
pars.FRIC.kind = 1;
pars.FRIC.Dc = 0.4;	% friction critical slip Dc
pars.FRIC.MUs = 0.677;	% static friction coefficient
pars.FRIC.MUd = 0.525;	% dynamic friction coefficient
pars.HEAL = 1;		% 0=no healing, 1=immediate healing upon arrest

%-- output settings --
%pars.VDYN = 1e-3;	% for the definition of the rupture front:
%			% velocity_threshold
pars.OX_FLAG = 1;	% snapshots at selected times
pars.OX_IDT = 2*pars.REFT;	% time stride for large outputs (snapshots)
pars.OX_IX = (pars.NX/2:2*pars.REFT:pars.NX);	% nodes for large outputs

pars.OT_FLAG = 1;	% timeseries at selected points
pars.OT_IX = round([0:4]*1e3*pars.NX/pars.L) +pars.NX/2;% nodes for full time outputs
pars.ECHO = 2;	% output info

%-- expert numerical settings --
pars.RTCUT = 2 ;% kernel_time-cut / fault_travel_time for mode_1, typically = 2
pars.QW = 4 ; 	% kernel_time-cut Nyquist_mode / mode_1 /(NX/2), typically = 4
                % NX/2 gives a non-truncated kernel
pars.CFL = 0.5 ;	% stability number (Courant), typically = 0.5


%====================================================================
%
% W3.m	Non-dimensional function W_{III} in mode III velocity kernel
%	= 1- Integral of C_{III} 
%	The integration is done numerically using the trapezoidal rule
%
% USES	C3.m
%
function WIII = W3(ARG)

  N = length(ARG);  % number of evaluation points
  nsub = 20; % number of sub-intervals for numerical integration
  nsub = max(nsub, nsub*(ARG(2)-ARG(1))/2/pi );
  Psub = nsub+1;

  if ARG(1)==0, error('FATAL ERROR in W3: ARG1 must be > 0'), end
  WIII = zeros(N,1);

 % first value, 
 % use same 'nsub' because usually ARG(1) = delay =~ ARG_step
  dsub = ARG(1)/nsub ;
  arg = (0:nsub).'*dsub ;
  CIII(1) = 0.5 ;  % = C3(0)
  CIII(2:Psub) = C3(arg(2:Psub)) ;
  WIII(1) = dsub*(0.5*CIII(1) + sum(CIII(2:nsub)) + 0.5*CIII(Psub) );

 % remaining values
  for i = 2:N
    dsub = (ARG(i)-ARG(i-1))/nsub ;
    arg  = ARG(i-1) + dsub*(0:nsub).' ;
    CIII = C3(arg) ;
    WIII(i) = WIII(i-1) ...
      + dsub*(0.5*CIII(1) + sum(CIII(2:nsub)) + 0.5*CIII(Psub) );
  end  

  WIII = 1 - WIII ;

%--------------------------------------------------------------------
%C3.m	Non dimensional kernel function for Mode III
%       C_{III} in Geubelle and Rice Eq.26
% WARNING: T=0 cannot be obtained with this formula
%	   (division by zero) and it is not implemented
%	   as special case for performance reasons
%	   Whenever needed: C3(0) = 0.5
function c = C3(T)
  c = besselj(1,T) ./ T ;

%--------------------------------------------------------------------
function ox = new_ox(IX,IDT,X,NT,DT)

  nt = floor(NT/IDT);
  nx = length(IX);

  ox.X 		= X(IX);
  ox.Time 	= (1:nt)*IDT*DT;
  ox.Slip	= zeros(nx,nt);
  ox.SlipRate 	= zeros(nx,nt);
  ox.Stress 	= zeros(nx,nt);
  ox.Strength 	= zeros(nx,nt);
  ox.RuptureTime = inf*ones(nx,1);
  ox.ProcessZone = inf*ones(nx,1);

%--------------------------------------------------------------------
function ot = new_ot(IX,X,NT,DT)

  nx = length(IX);

  ot.X 		= X(IX);
  ot.Time 	= (1:NT)*DT;
  ot.SlipRate 	= zeros(NT,nx);
  ot.Slip 	= zeros(NT,nx);
  ot.Stress 	= zeros(NT,nx);
  ot.MeanSlip 	= zeros(NT,1);
  ot.MaxSlip 	= zeros(NT,1);
  ot.MeanSlipRate = zeros(NT,1);
  ot.MaxSlipRate = zeros(NT,1);
  ot.CrackLength = zeros(NT,1);
  ot.MeanStress = zeros(NT,1);

%--------------------------------------------------------------------
% zero crossing time, with linear time interpolation
% tc must be initialized at = inf
function tc = zeroCrossingTime(tc, v,v_old, it,dt)

  ix = find( v_old<0 & v>=0 & tc==inf);
  tc(ix) = dt*(it- v(ix)./(v(ix)-v_old(ix)) ); 
