close all
clear all

%% Inputs

tsv_diameter = 5e-6; % (m) Diameter of TSV metal core
l_tsv = 20e-6; % (m) TSV length
tox = 120e-9; % (m) Thickness of oxide liner

Na_cm = 2e15; % (cm^-3) Background substrate doping (p-type)
ni_cm = 1e10; % (cm^-3) Intrinsic carrier concentration in substrate (1e10 for Si)

epsr_si = 11.7; % (-) Relative permittivity of silicon
epsr_ox = 3.9; % (-) Relative permittivity of silicon dioxide

T = 300; % (K) Temperature


Qf_cm2 = 0; % (C/cm^2) Fixed charge at the Si/SiO2 interface
Qm_cm2 = 0; % (C/cm^2) Mobile charges in the oxide
Qot_cm2 = Qm_cm2 + Qf_cm2; % (C/cm^2) Total oxide charge

phi_cu = 4.7; % (eV) Work function of copper
chi_si = 4.05; % (eV) Electron affinity of silicon

Npts_v = 1e3;
Vtsv = linspace(-5, 5, Npts_v);
Ctsv = zeros(1, Npts_v);

%% Constants
eps0 = 8.854e-12; % (F/m) Vacuum permittivity
q = 1.602e-19; % (C) Electron charge
kb = 1.380e-23; % (J/K) Boltzmann constant

%% Conversions
eps_si = epsr_si * eps0; % (F/m)
eps_ox = epsr_ox * eps0; % (F/m)
Na = Na_cm * 100^3; % (m^-3) Background substrate doping (p-type)
ni = ni_cm * 100^3; % (m^-3) Intrinsic carrier concentration in substrate

E_bb = kb*T/q*log(Na/ni);
phi_si = chi_si + E_bb;

phi_ms = phi_cu - phi_si; % Work function difference between copper and silicon
Qot = Qot_cm2 * 100^2; % (C/m^2) Total oxide charge

Rmet = tsv_diameter/2; % (m) Radius of TSV metal core
Rox = Rmet + tox; % (m) Outer radius of TSV oxide liner


%% Useful quantities
a = q*Na/4/eps_si;
phi_f = kb*T/q*log(Na/ni);
beta = 2*phi_f/a - Rox^2;

%% Calculations

% Oxide capacitance
Cox = 2*pi*eps_ox*l_tsv/log(Rox/Rmet); % (F) oxide capacitance

% Approximate value for maximum depletion width
Rmax_guess = sqrt(beta)/sqrt(lambertw(0.367879*beta/Rox^2)); % guess from wolfram alpha solution of Rmax equation

% Function for finding more accurate Rmax
rmax_func = @(Rmax) abs( a*Rox^2 - 2*a*Rmax.^2.*log(Rox) + a*Rmax.^2.*(2*log(Rmax) - 1) - 2*phi_f );

% Find Rmax using either fminbnd or fminsearch
% Fminbnd may be faster, but fminsearch may avoid truncating search space
opts = optimset('TolX', 1e-16);
[Rmax, fval] = fminbnd( rmax_func, 0.1*Rmax_guess, 10*Rmax_guess, opts);
%[Rmax, fval] = fminsearch( rmax_func, Rmax_guess, opts);


Vfb = phi_ms - Rox*q*Qot*log(Rox/Rmet)/eps_ox;
%Vth = Vfb + 2*v_thermal*log(Na/ni) + q*Na*(Rmax^2-Rox^2)/2/eps_ox*log(Rox/Rmet); % Formula from Katti, possibly incorrect or inapplicable(Vth = Vfb - Qd/Ci + 2*phi_f)
Vth = Vfb + q*Na*(Rmax^2-Rox^2)/2/eps_si*log(Rmax/Rox); % Modified formula: Vth = Vfb - Qd/Cd

% Alt formulation: Vt = phi_ms + Qi/Ci + Qd/Ci + 2*phi_f
% Matches Above when Vtdep depends on Cdep_min
Cdep_min = 2*pi*eps_si*l_tsv/log(Rmax/Rox);
qdep_max = q*Na*pi*(Rmax^2-Rox^2)*l_tsv;
Vtdep = qdep_max/Cdep_min;
%Vtdep = qdep_max/Cox;
Vth2 = Vfb + Vtdep;% + 2*phi_f;


%% Find Rdep for appropriate voltages
Rdep = zeros(1, Npts_v );
vind_dep_start = find( Vtsv >= Vfb, 1, 'first');
vind_dep_stop = find(Vtsv <= Vth, 1, 'last');

for vind = vind_dep_start:vind_dep_stop
    vtsv = Vtsv(vind);
    
    rdep_func = @(Rdep) abs( 2*a*(Rdep.^2 - Rox^2).*log(Rdep/Rox) - (vtsv - Vfb) );
%     rdep_func = @(Rdep) abs( 2*a*(Rdep.^2 - Rox^2).*log(Rdep/Rox) - (vtsv - Vfb - 2*phi_f) );
%     rdep_func = @(Rdep) abs( q*Na/2/eps_ox*(Rdep.^2 - Rox^2).*log(Rox/Rmet) - (vtsv - Vfb - 2*phi_f) );
    rdep = fminsearch( rdep_func, Rox, opts);
    %[rdep, fval] = fminbnd( rdep_func, Rox, Rmax, opts);
    
    Rdep(vind) = rdep;
end

%% Capacitance
% Accumulation region (Vtsv < Vfb)
Ctsv(Vtsv <= Vfb) = Cox;

% Depletion ( Vfb <= Vtsv <= Vth)
Cdep = 2*pi*eps_si*l_tsv./log(Rdep/Rox);
in_depletion = ((Vtsv >= Vfb) & (Vtsv < Vth));
%Ctsv(in_depletion) = Cox*Cdep(in_depletion)./(Cox + Cdep(in_depletion));
Ctsv(in_depletion) = (1/Cox + 1./Cdep(in_depletion) ).^-1;

% Inversion (Vtsv >= Vth)
Cdep_min = 2*pi*eps_si*l_tsv/log(Rmax/Rox);
Ctsv_min = Cox*Cdep_min/(Cox+Cdep_min);
Ctsv(Vtsv >= Vth) = Cox*Cdep_min/(Cox+Cdep_min);

fprintf('Rmax: %.3g um\n', Rmax*1e6)
fprintf('Vth: %.3g V\n', Vth)
fprintf('Vfb: %.3g V\n', Vfb)
fprintf('Cox: %.3g fF\n', Cox*1e15)
fprintf('Cdep_min: %.3g fF\n', Cdep_min*1e15)
fprintf('Ctsv_min: %.3g fF\n', Ctsv_min*1e15)
fprintf('\n')

%% Plots

figure(1)
clf
hold on
plot(Vtsv, Ctsv*1e15)
plot(Vtsv(vind_dep_start:vind_dep_start+1), [min(Ctsv) max(Ctsv)]*1e15,'k:')
plot(Vtsv(vind_dep_stop-1:vind_dep_stop), [min(Ctsv) max(Ctsv)]*1e15,'k:')
xlabel('DC Bias (V)')
ylabel('TSV Capacitance (fF)')
%set(gca,'yscale','log')
fixfigs(1,3,14,12)

figure(4)
clf
hold on
plot(Vtsv(in_depletion), Ctsv(in_depletion)*1e15,'k')
plot(Vtsv(in_depletion), Cox*ones(1,length(Ctsv(in_depletion)))*1e15, 'b--')
plot(Vtsv(in_depletion), Cdep(in_depletion)*1e15, 'r')
plot(Vtsv(in_depletion), Cdep_min*ones(1,length(Ctsv(in_depletion)))*1e15, 'r--')
xlabel('DC Bias (V)')
ylabel('TSV Capacitance (fF)')
%set(gca,'yscale','log')
xlim( [min(Vtsv(in_depletion)) max(Vtsv(in_depletion)) ])
fixfigs(4,3,14,12)

figure(5)
clf
hold on
plot(Vtsv(in_depletion), Rox*ones(1,length(Ctsv(in_depletion)))*1e6, 'k')
plot(Vtsv(in_depletion), Rdep(in_depletion)*1e6,'b')
plot(Vtsv(in_depletion), Rmax*ones(1,length(Ctsv(in_depletion)))*1e6, 'r--')
xlabel('DC Bias (V)')
ylabel('Depletion Radius (um)')
%set(gca,'yscale','log')
xlim( [min(Vtsv(in_depletion)) max(Vtsv(in_depletion)) ])
fixfigs(5,3,14,12)


