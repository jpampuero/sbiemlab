% FRICTION slip-dependent friction coefficient
%
% mu = friction(u,f) gives the friction coefficient
% friction() plots the friction laws (demo)
%
% INPUT     u       slip (positive) or state variable (sort of slip variable
%                   externally reset to zero at healing phases)
%           f       friction data structure, contains the following fields:
%           MUs     static friction coefficient
%           MUd     dynamic friction coefficient
%           Dc      critical slip distance
%           p       exponent
%           kind    type of friction law:
%                   1 linear slip-weakening
%                       mu(u) = MUs - (MUs-MUd)*u/Dc  if u<Dc
%                             = MUd otherwise
%                   2 non-linear (power-law) friction law of Chambon et al. (JGR 2006b)
%                       mu(u) = MUd +(MUs-MUd)/(1+u/(p*Dc))^p
%                   3 initially linear slip-weakening law,
%                     then non-linear law of Abercrombie and Rice (GJI 2005)
%                       mu(u) = MUs - (MUs-MUd)*u/Dc  if u<Dc
%                             = MUs - (MUs-MUd)*(p-1+(u/Dc)^p)/p  otherwise
%                   4 exponential slip-weakening
%                       mu(u) = MUd + (MUs-MUd)*exp(-u/Dc)
%                   5 non-linear (power-law) friction law of Chambon et al. (JGR 2006b)
%                       mu(u) = MUd + (MUs-MUd)/(1+u/Dc)^p;
%                       NOTE: the slope at u=0 depends on p with this normalization (u/Dc).
%                             Option 2 (with normalization u/(p*Dc)) has slope independent of p
%                   6 dual-scale linear slip-weakening
%                       mu(u) = MUs - (MUs-MUi) *    u   /  Dci     if       u <= Dci
%                       mu(u) = MUi - (MUi-MUd) * (u-Dci)/(Dc-Dci)  if Dci < u <= Dc
%                       mu(u) = MUd                                 if Dc  < u
%                       NOTE: the slope of the first stage must be steeper 
%                             than that of the second stage, i.e.:
%                             (MUs - MUi)/Dci > (MUi - MUd)/Dc 
%                             also, Dci < Dc
%
% OUTPUT    mu  friction coefficient
%
% NOTE  Vector inputs: 'u' and any field of 'f' can be a vector, with same size.
%   'mu' is then a vector output, of same size.
%   'f' itself cannot be a vector.
%
function mu = friction(u,f)

% if no arguments: plot the friction laws
if nargin==0
  mu = plotFriction(); 
  return; 
end

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
%                                 20   : 80%
%                                 126  : 90%
%                                 4e4  : 99%
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

 case 5 %-- power-law slip weakening friction law
  mu = f.MUd + (f.MUs-f.MUd)./(1+u./f.Dc).^(f.p);

 case 6 %-- dual-scale linear slip-weakening 
  W1 = (f.MUs - f.MUi)./f.Dci;
  W2 = (f.MUi - f.MUd)./(f.Dc-f.Dci);

  mu1 = f.MUs - W1.* u;
  mu2 = f.MUi - W2.*(u-f.Dci);

  mu = max(mu1,max(mu2,f.MUd));
  
 otherwise
  error('friction','unknown kind')

end

%------- demo ----------
% plot the friction laws 
function f = plotFriction()

f.MUs = 0.6;
f.MUd = 0.5;
f.Dc = 1;

figure;
f.kind = 1;
u=[0:2];
subplot(321)
plot(u,friction(u,f))
axis([0 2 0.4 0.7])
ylabel('Friction coefficient \mu');
xlabel('Slip [m]');
title('1: Linear slip-weakening')

subplot(323)
u = logspace(-3,2,100);
f.kind=2;
for p=[0.1:0.1:0.4]
  f.p = p; 
  plot(u,friction(u,f))
  hold on
end
hold off
axis([0 inf 0.4 0.7])
legend('p=0.1','0.2','0.3','0.4')
title('2: Chambon''s law (1+slip)^{-p}, normalization: slip/(p\timesDc)')
ylabel('Friction coefficient \mu');
xlabel('Slip [m]');

subplot(325)
f.kind=3;
for p=[0.1:0.1:0.4]
  f.p = p; 
  plot(u,friction(u,f))
  hold on
end
hold off
axis([0 inf 0. 0.7])
ylabel('Friction coefficient \mu');
xlabel('Slip [m]');
legend('p=0.1','0.2','0.3','0.4')
title('3: Regularized Abercrombie and Rice''s law  a-b\timesslip^p')

subplot(322)
f.kind=4;
plot(u,friction(u,f))
axis([0 2 0.4 0.7])
title('4: Exponential slip-weakening')
ylabel('Friction coefficient \mu');
ylabel('Friction coefficient \mu');
xlabel('Slip [m]');

subplot(324)
f.kind=5;
f.MUd = 0;
for p=[1/3,3]
  f.p = p;
  plot(u,friction(u,f));
  hold on
end
Lines = colormap(lines);
f.kind=2;
p=[1/3,3];
for i=1:2 
  f.p = p(i);
  plot(u,friction(u,f),':','Color',Lines(i,:));
  hold on
end
hold off
axis([0 10 0. 0.7])
legend('p = 1/3','3','Option 2, p = 1/3','3');
ylabel('Friction coefficient \mu');
xlabel('Slip [m]');
title('5 Chambon''s law (1+slip)^{-p}, normalization slip/Dc')

subplot(326)
f.kind=6;
f.MUs = 0.677;
f.MUi = 0.6;
f.MUd = 0.525;
f.Dci = 0.05;
f.Dc  = 0.4;
u = linspace(0,1.5*f.Dc,201);
plot(u,friction(u,f))
xlabel('Slip');
ylabel('Friction coefficient \mu');
title('6: Dual-scale LSW')
