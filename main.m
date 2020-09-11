% Written by KyeoReh Lee, Ph.D., kyeo@kaist.ac.kr, 2020.09.03
% Original paper: TBD

%% set path
cd(fileparts(matlab.desktop.editor.getActiveFilename))

%% constants
h = 6.62607004e-34;         % Js, Plank constant
c = 2.998e8;                % m/s, speed of light
lambda_ex = 804.6e-9;       % m, excitation wavelength
lambda_em = 1064e-9;        % m, emmision wavelength
tau_2 = 230e-6;             % s, Effective upper level lifetime. 85us - 4% , 225us - 1%
sigma_st = 6.5e-23;         % m^2, Nd:YAG stimulated emission cross section @ 1064 nm
% Semwal, K. and S. Bhatt, Study of Nd3+ ion as a Dopant in YAG and Glass Laser. International Journal of Physics, 2013. 1(1): p. 15-21.

%% fiber parameters

% thoralbs, FP200ERT
coreDia =  200e-6;      % core diameter of fiber, m
cladDia =  225e-6;      % cladding diameter of fiber, m

NA_fiber = 0.5;         % numerical aperture of fiber
RI_fiber = 1.4496;      % RI of fused sillica @ 1064 nm

% no fiber, just an aperture
% coreDia =  10e-6;  % core diameter of fiber, m
% cladDia =  10e-6;  % cladding diameter of fiber, m
% NA_fiber = 1;     % numerical aperture of fiber
% RI_fiber = 1;

%% integrating sphere parameters
apertureDia = 250e-6;      %m, diameter of aperture.

mu_a_ex     = 20;          % m^-1, absorption coefficient of excitation light
mu_a_em     = 2.3e-1;      % m^-1, absorption coefficient of emission light

remu_s_ex   = 356e3; 	   % m^-1, (1-g_em)*mu_s_em, reduced scattering coefficient of excitation light
remu_s_em   = 317e3;       % m^-1, (1-g_em)*mu_s_em, reduced scattering coefficient of emission light

volfrac     = 0.75;        % volume fraction of Nd:YAG.
RI_medium   = 1;           % refractive index of medium (air) @ emission wavelength;
RI_NdYAG    = 1.8169;      % refractive index of Nd:YAG = 1.8169 @ 1064 nm; 1.8217 @ 808 nm
% Semwal, Kireet, and S. C. Bhatt. "Study of Nd3+ ion as a Dopant in YAG and Glass Laser." International Journal of Physics 1.1 (2013): 15-21.

couplingEff_fiber = (coreDia/apertureDia*NA_fiber)^2;   % out coupling efficiency of the channel.
effRI = (RI_NdYAG-RI_medium)*volfrac + RI_medium;       % estimated effective refractive index by the Nd:YAG volume fraction

asRatio_ex = mu_a_ex/remu_s_ex;                         % ratio btw absorption and reduced scattering coefficient of excitation light
asRatio_em = mu_a_em/remu_s_em;                         % ratio btw absorption and reduced scattering coefficient of emission light

Refl_oc = FresnelMeanRefl(1/RI_fiber);                  % the outcouping reflectivity estimated the fresnel coefficients.
Refl_NdYAG_specular = FresnelMeanRefl(1/RI_NdYAG);      % specular relection from the wall surface due to the RI mismatch)
effPenDepth = (1-Refl_NdYAG_specular)*(1/remu_s_em);    % m, Eq.(S10), effective wall penetration depth
tau_wall = 3/2*effRI/c * (1/remu_s_em)*(1+asRatio_em)/(1+sqrt(3*asRatio_em*(1+asRatio_em))); %s, Eq.(S10), mean dwell time for single reflection
% Patterson, M.S., B. Chance, and B.C. Wilson, Time resolved reflectance and transmittance for the noninvasive measurement of tissue optical properties. Applied Optics, 1989. 28(12): p. 2331-2336.

%% parameters
% costants
Refl_em = 0.989;                            % measured wall reflectivity for emission light   < 1
Refl_ex = 0.913;                            % measured wall reflectivity for excitation light < 1

% variables
Refl_noc = 0.1*5+0.01*-28;                  % effective reflectance of cylinderical channel surface + fiber cladding.
p_leak = 0.01*0+0.001*9;                    % unexpected out-coupling chance (area fraction)
eff_prac_ex = 1-(0.1*5+0.01*-4);            % practical efficiency of the pumping process  (mainly from the fiber loss)
eff_prac_em = 1-(0.1*4+0.01*2);             % practical efficiency of the emission process (mainly from the fiber loss)

% main
outParams = struct([]);
meanTransLib = load('meanTransportLength.mat');         % m, numerically calculated mean trasport length; rad, theta_fiber, Supplementary Figure 5.
cavDia =  apertureDia*linspace(1,7,1500);               % scattering cavity diameter as a main variable.
N = length(cavDia);

for kk = 1:1:N
    cavRadius   = cavDia(kk)/2;
    beta        = cavDia(kk)/apertureDia;
    theta_fiber = asin(1/beta);   
    
    % statistical probablities
    p_wall = 2./(3-cos(theta_fiber))-p_leak;               % Eq.(4)
    p_aper = (1-cos(theta_fiber))./(3-cos(theta_fiber));   % Eq.(4)   
    p_oc  = couplingEff_fiber*p_aper;                      % Eq.(5)
    p_noc = p_aper-p_oc;                                   % Eq.(5)

    % effective reflectances
    eff_Refl_em = p_wall*Refl_em + p_oc*Refl_oc + p_noc*Refl_noc;                       % Eq.(S1)
    eff_Refl_ex = p_wall*Refl_ex + p_oc*Refl_oc + p_noc*Refl_noc;                       % Eq.(S1)
    
    % mode volume
    V_phi = pi/24 * (cavDia(kk)+2*effPenDepth)^3 * (2-cos(theta_fiber)) * (1+cos(theta_fiber))^2; % m^3, Eq.(9)

    % mean transport time per wall reflection
    meanTransL = interp1(meanTransLib.theta_fiber,real(meanTransLib.meanTransporLength),theta_fiber,'pchip') * cavRadius;  % m, numerically calculated mean trasport length.
    tau_cav    = -(1/log(eff_Refl_em)) * (meanTransL/c+tau_wall); %s, Eq.(10), cavity lifetime, 
    
    % outcopuling and pumping efficincy
    eff_out  = p_oc  * (1-Refl_oc) / (1 - eff_Refl_em);             % Eq.(S2)
    eff_pump = p_wall * (1-Refl_ex)  / (1 - eff_Refl_ex);           % Eq.(S2)
        
    % pumping threshold and slope efficincy
    P_th = 1/eff_pump/eff_prac_ex * 1/tau_cav/tau_2 * ( h/lambda_ex/sigma_st*V_phi );   % W, Eq.(8)
    slopeEff = (eff_prac_em*eff_prac_ex) *eff_pump *eff_out *(lambda_ex/lambda_em);     % Eq.(1)
    
    % making a structure
    outParams(kk).cavDia        = cavDia(kk);
    outParams(kk).beta          = beta;    
    outParams(kk).tau_cav       = tau_cav;    
    outParams(kk).eff_Refl_em   = eff_Refl_em;
    outParams(kk).eff_Refl_ex   = eff_Refl_ex;
    outParams(kk).eff_pump      = eff_pump;
    outParams(kk).eff_out       = eff_out;    
    outParams(kk).P_th          = P_th;
    outParams(kk).slopeEff      = slopeEff;
end

%% plot
yyaxis left
plot([outParams.beta],[outParams.slopeEff])
xlabel('\beta (cavity diameter / channel diamter)');
ylabel('slope efficiency')

yyaxis right
plot([outParams.beta],[outParams.P_th])
ylabel('laser threshold (W)')
axis square

% theoretical max value calculation.
eff_Refl_aper = couplingEff_fiber * Refl_oc + (1 - couplingEff_fiber) * Refl_noc;

q           = 1 - p_leak;
gamma_ex    = sqrt(1 - q*Refl_ex);
gamma_em    = sqrt(1 - q*Refl_em);
gamma_a     = sqrt(1 - q*eff_Refl_aper);

beta_star = ...
    ( gamma_a^2 + (1-q)*gamma_ex*gamma_em ) /...
    ( 2 * sqrt( q*gamma_ex*gamma_em * ( gamma_a^2 + (1-2*q) * gamma_ex*gamma_em) ) );   % Eq. (S9)

eff_pump_eff_out_max = ...
    couplingEff_fiber *q^2 * (1 - Refl_oc) * (1 - Refl_ex) / ( gamma_a * ( gamma_ex + gamma_em ) )^2; % Eq. (S8)

slopeEff_max = (eff_prac_em*eff_prac_ex) *eff_pump_eff_out_max *(lambda_ex/lambda_em);  % Eq.(1)

yyaxis left, hold on
plot(beta_star,slopeEff_max,'*')
hold off

