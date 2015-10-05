close all
clear all

%% Inputs

Rmet = 5e-6; % (m) Radius of TSV metal core
l_tsv = 20e-6; % (m) TSV length
tox = 50e-9; % (m) Thickness of oxide liner

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
Vtsv = linspace(-10, 10, Npts_v);
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

Rox = Rmet + tox; % (m) Outer radius of TSV oxide liner


%% Useful quantities
alpha = q*Na/4/eps_si; % Prefactor in the Rmax calculation
v_thermal = kb*T/q; % (V) Thermal voltage
beta = 2*kb*T/q/alpha* log(Na/ni) - Rox^2;

%% Calculations

Cox = 2*pi*eps_ox*l_tsv/log(Rox/Rmet); % (F) oxide capacitance

Rmax = sqrt(beta)/sqrt(lambertw(0.367879*beta/Rox^2));

a = 2*kb*T/q/alpha * log(Na/ni);
b = Rox;
Rmax_alt = b*exp( 1/2*lambertw( (a-b^2)/exp(1)/b^2 ) + 1);
%Rmax = Rmax_alt;

Vfb = phi_ms - Rox*q*Qot_cm2*log(Rox/Rmet)/eps_ox;
%Vth = phi_ms - Vfb + 2*v_thermal*log(Na/ni) + q*Na*(Rmax^2-Rox^2)/2/eps_ox*log(Rox/Rmet);
Vth = phi_ms - Vfb + 2*v_thermal*log(Na/ni) + q*Na*(Rmax^2-Rox^2)/2/eps_si*log(Rmax/Rox);

% TSV voltage in depletion region (Vfb <= Vtsv <= Vth)
gamma = q*Na/2/eps_ox*log(Rox/Rmet); % prefactor
gamma_prime = q*Na/2/eps_si; % prefactor, corrected

%Rdep = sqrt( Rmax^2  - (Vth-Vtsv)/gamma);

%% Find Rdep for appropriate voltages
rr = linspace(Rox, Rmax, 1e2);
f4 = figure(4);
clf
hold on
f5 = figure(5);
clf
hold on

Rdep = zeros(1, Npts_v );
vind_dep_start = find( Vtsv >= Vfb, 1);
vind_dep_stop = find(Vtsv <= Vth, 1, 'last');
for vind = vind_dep_start:vind_dep_stop
    vtsv = Vtsv(vind);
    solve_func = @(Rdep) (Rdep.^2 - Rox^2).*log(Rdep/Rox) - 1/gamma_prime*(vtsv-Vfb);
    solve_func_real = @(Rdep) real(solve_func(Rdep));
    solve_func_imag = @(Rdep) imag(solve_func(Rdep));
    roots_real = fzero(solve_func_real, Rmax);
    roots_imag = fzero(solve_func_imag, Rmax);
    Rdep(vind) = roots_real;
    
    dd = solve_func(rr);
    figure(4)
    plot(rr, real(dd));
    figure(5)
    plot(rr, imag(dd));
end

%Rdep(vind_stop+1:end) = Rmax;



%% Capacitance
% Accumulation region (Vtsv < Vfb)
Ctsv(Vtsv < Vfb) = Cox;

% Depletion ( Vfb <= Vtsv <= Vth)
Cdep = 2*pi*eps_si*l_tsv./log(Rdep/Rox);
in_depletion = ((Vtsv >= Vfb) & (Vtsv <= Vth));
Ctsv(in_depletion) = Cox*Cdep(in_depletion)./(Cox + Cdep(in_depletion));

% Inversion (Vtsv >= Vth)
Cdep_min = 2*pi*eps_si*l_tsv/log(Rmax/Rox);
Ctsv_min = Cox*Cdep_min/(Cox+Cdep_min);
Ctsv(Vtsv >= Vth) = Cox*Cdep_min/(Cox+Cdep_min);



fprintf('Rmax: %.3g um\n', Rmax*1e6)
fprintf('Rmax_alt: %.3g um\n', Rmax_alt*1e6)
fprintf('Vth: %.3g V\n', Vth)
fprintf('Vfb: %.3g V\n', Vfb)
fprintf('Cox: %.3g fF\n', Cox*1e15)
fprintf('Cdep_min: %.3g fF\n', Cdep_min*1e15)
fprintf('Ctsv_min: %.3g fF\n', Ctsv_min*1e15)

%% Plots

figure(1)
clf
plot(Vtsv, Ctsv*1e15)
xlabel('DC Bias (V)')
ylabel('TSV Capacitance (fF)')
%set(gca,'yscale','log')
fixfigs(1,3,14,12)

figure(2)
clf
plot(Vtsv, imag(Ctsv))
xlabel('DC Bias (V)')
ylabel('TSV Capacitance')
fixfigs(1,3,14,12)

figure(3)
clf
plot(Vtsv(in_depletion), imag(Rdep(in_depletion)))
xlabel('DC Bias (V)')
ylabel('Depletion_width')
fixfigs(1,3,14,12)



