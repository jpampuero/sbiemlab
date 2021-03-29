% FRICTION slip-dependent friction coefficient
%
% mu = friction(u,f) gives the friction coefficient
% friction() plots the friction laws (demo)
%
% INPUT		u	slip (positive) or state variable (sort of slip variable
%			externally reset to zero at healing phases)
%   		f	friction data structure, contains the following fields:
%			MUs	static friction coefficient
%			MUd	dynamic friction coefficient
%			Dc	critical slip distance
%			p	exponent
%			kind	type of friction law:
%				1 linear slip-weakening
%				   mu(u) = MUs - (MUs-MUd)*u/Dc  if u<Dc
%                                        = MUd otherwise
%				2 non-linear law of Chambon et al. (JGR 2006b)
%				   mu(u) = MUd +(MUs-MUd)/(1+u/Dc)^p
%				3 initially linear slip-weakening law,
%				  then non-linear law of Abercrombie and Rice (GJI 2005)
%				   mu(u) = MUs - (MUs-MUd)*u/Dc  if u<Dc
%				         = MUs - (MUs-MUd)*(p-1+(u/Dc)^p)/p  otherwise
%				4 exponential slip-weakening
%				   mu(u) = MUd + (MUs-MUd)*exp(-u/Dc)
%
% OUTPUT	mu 	friction coefficient
%
% NOTE	Vector inputs: 'u' and any field of 'f' can be a vector, with same size.
%	'mu' is then a vector output, of same size.
%	'f' itself cannot be a vector.
%
function mu = friction(u,f)

% if no arguments: plot the friction laws
if nargin==0, mu = plotFriction(); return; end

% default = linear slip weakening
if ~isfield(f,'kind'), f.kind=1; end

switch f.kind

 case 1 % Linear slip weakening friction law
  W = (f.MUs-f.MUd) ./ f.Dc;
  mu = max(f.MUs - W.*u, f.MUd);

 case 2 %-- Chambon's non linear law --
% for large slip: strength ~ 1/u^p
  u = u ./(f.p*f.Dc);
  mu = f.MUd +(f.MUs-f.MUd) ./ (1+u).^f.p ; 
% With this definition:
% initial slope = (MUs-MUd)/Dc
% if p=0.4, strength drop at slip=2*Dc : 50% 
%	     			  20   : 80%
%	     			  126  : 90%
%            			  4e4  : 99%
% so if Dc=100e-6 m: 90% drop at slip = 1.2 cm
%                    99% ...          = 4 m

 case 3 %-- Abercrombie and Rice non-linear law --
% for large slip: strength_drop ~ u^p
  u = u ./ f.Dc;
  mu = f.MUs -(f.MUs-f.MUd).* ( u.*(u<=1) + (f.p-1+u.^f.p)./f.p .*(u>1) );
% With this definition: initial slope = (MUs-MUd)/Dc
% p =0.3 is usual (?)

 case 4 %-- exponential slip weakening friction law
  mu = f.MUd + (f.MUs-f.MUd).*exp(-u./f.Dc);

 otherwise
  error('friction','unknown kind')

end

%------- demo ----------
% plot the friction laws 
function f = plotFriction()

f.MUs = 0.6;
f.MUd = 0.5;
f.Dc = 1;

clf
f.kind = 1;
u=[0:2];
subplot(311)
plot(u,friction(u,f))
axis([0 2 0.4 0.7])
ylabel('Friction coefficient \mu')
title('Linear slip-weakening')

subplot(312)
u = logspace(-3,2,100);
f.kind=2;
for p=[0.1:0.1:0.4],
  f.p = p; 
  plot(u,friction(u,f))
  hold on
end
hold off
axis([0 inf 0.4 0.7])
legend('p=0.1','0.2','0.3','0.4')
title('Chambon''s law (1+slip)^{-p}')

subplot(313)
f.kind=3;
for p=[0.1:0.1:0.4],
  f.p = p; 
  plot(u,friction(u,f))
  hold on
end
hold off
axis([0 inf 0. 0.7])
xlabel('Slip')
legend('p=0.1','0.2','0.3','0.4')
title('Regularized Abercrombie and Rice''s law  a-b\timesslip^p')
