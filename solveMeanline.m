function [meanline, span] = solveMeanline(inputs)

%inputs
P = inputs.P; %shaft power W
mdot = inputs.mdot; %mass flow rate kg/s
d_m = inputs.d_m/1000; %mean diameter m
omega = inputs.omega; %rotational speed rev/s
reaction = inputs.reaction;
beta_2 = inputs.beta_2; %deg
T_01 = inputs.T_01; %K
c_p = inputs.c_p; %specific heat at const P J/kgK
gamma = inputs.gamma; %ratio of specific heats
Mmol = inputs.Mmol; %kg/kmol
p_01 = inputs.p_01; %total pressure nozzle inlet Pa
degree_of_admission = inputs.degree_of_admission; 
delta_3 = inputs.delta_3;  %rotor inlet flow deviation deg
gamma_ri = inputs.gamma_ri; %rotor inlet wedge angle 0 < gamma_ri < 2*beta_3 deg
eta_m = inputs.eta_m; %mechanical efficiency
eta_n = inputs.eta_n; %nozzle efficiency
n = inputs.n; %vortex law exponent n
m = inputs.m; %vortex law exponent m
R_s = 8314.462618/Mmol;


%specific work
deltah_0 = P/(eta_m*mdot);

%blade speed 
u = pi*d_m*omega;

%flow velocities
c_3u = u*(1-reaction) + deltah_0/(2*u);
c_4u = u*(1-reaction) - deltah_0/(2*u);
w_3u = c_3u - u;
w_4u = c_4u - u;

%velocity triangles
%flow velocities
c_3 = c_3u/cosd(beta_2);
c_m = c_3*sind(beta_2); %constant meridional velocity assumption
c_4 = sqrt((c_m^2) + (c_4u^2));
w_3 = sqrt((c_m^2) + (w_3u^2));
w_4 = sqrt((c_m^2) + (w_4u^2));

%flow angles
alpha_3 = atan2d(c_3u, c_m);
alpha_3_r = atan2d(w_3u, c_m); %r for relative flow velocity
alpha_4_r = atan2d(w_4u, c_m);
alpha_4 = atan2d(c_4u, c_m);
alpha_e_r = alpha_3_r + delta_3 + (0.5*gamma_ri);

%metal angles
beta_3 = 90 - alpha_e_r + (0.5*gamma_ri);
beta_4 = 90 + alpha_4_r;

%total and static temps and pressures
T_03 = T_01;
T_3 = T_01 - c_3^2/(2*c_p);
T_3s = T_01 - c_3^2/(2*eta_n*c_p);
T_04 = T_01 - deltah_0/c_p;
T_4 = T_04 - c_4^2/(2*c_p);

p_3 = p_01*(T_3s/T_01)^(gamma/(gamma-1));
p_03 = p_3*(T_03/T_3)^(gamma/(gamma - 1));
p_4 = p_3;
p_04 = p_4*(T_04/T_4)^(gamma/(gamma - 1));
T_04s = T_01*(p_04/p_01)^((gamma-1)/gamma);
eta_tt = (T_01 - T_04)/(T_01 - T_04s); %turbine total to total efficiency

assert(p_03<=p_01, 'p_03 exceeds p_01. Nozzle efficiency may be unrealistic.');
assert(p_04<=p_03, 'p_04 exceeds p_03. Turbine efficiency may be unrealistic.');

%mach numbers
a_3 = sqrt(gamma*R_s*T_3);
a_4 = sqrt(gamma*R_s*T_4);
M_3 = c_3/a_3;
M_3_r = w_3/a_3;
M_4 = c_4/a_4;
M_4_r = w_4/a_4;

%densities
rho_3 = p_3/(R_s*T_3);
rho_2 = rho_3;
rho_4 = p_4/(R_s*T_4);

%blade heights m
H_2 = mdot/(degree_of_admission*rho_2*pi*d_m*c_m);
H_3 = mdot/(degree_of_admission*rho_3*pi*d_m*c_m);
assert(H_3 > 3.8e-3, 'H_3 is below min blade height of 3.8 mm.');
H_4 = mdot/(degree_of_admission*rho_4*pi*d_m*c_m);
assert(H_4 > 3.8e-3, 'H_4 is below min blade height of 3.8 mm.');


%span calculations

%radius
r_m = d_m/2;
r = linspace(r_m -H_3/2, r_m + H_3/2, 21);

%mean radius coefficients
R_m = reaction;
psi_m = deltah_0/u^2;

%speeds
u_span = u./(r_m./r);
c_3u_span = u*((1 - R_m)*(r_m./r).^n + (psi_m/2)*(r_m./r).^m);
c_4u_span = u*((1 - R_m)*(r_m./r).^n + (psi_m/2)*(r_m./r).^m - psi_m*(r_m./r));

%solving equations for radial equilibrium
dc_3m_squared_dr = -(2*u^2./r).*((1-R_m)*(r_m./r).^n + (psi_m/2)*(r_m./r).^m).*((1-R_m)*(1-n)*(r_m./r).^n + (psi_m/2)*(1-m)*(r_m./r).^m);
dc_4m_squared_dr = -(2*u^2./r).*((1-R_m)*(r_m./r).^n + (psi_m/2)*(r_m./r).^m - psi_m*(r_m./r)).*((1-R_m)*(1-n)*(r_m./r).^n+ (psi_m/2)*(1-m)*(r_m./r).^m);

integrated_c_3m_squared = cumtrapz(r, dc_3m_squared_dr);
integrated_c_4m_squared = cumtrapz(r, dc_4m_squared_dr);
mean_radius_i = (length(r) + 1)/2;

c_3m = sqrt(c_m^2 + integrated_c_3m_squared- integrated_c_3m_squared(mean_radius_i));
c_4m = sqrt(c_m^2 + integrated_c_4m_squared - integrated_c_4m_squared(mean_radius_i));

%relative whirl velocitie across span
w_3u_span = c_3u_span - u_span;
w_4u_span = c_4u_span - u_span;

%relative flow angles across span
alpha_3_r_span = atan2d(w_3u_span, c_3m);
alpha_4_r_span = atan2d(w_4u_span, c_4m);

alpha_e_r_span = alpha_3_r_span + delta_3 + 0.5*gamma_ri;
beta_3_span = 90 - alpha_e_r_span + (0.5*gamma_ri);
beta_4_span = 90 + alpha_4_r_span;

%absolute flow angle across span
alpha_3_span = atan2d(c_3u_span, c_3m);

beta_2_span = 90 - alpha_3_span;

deltah_span = u_span.*(c_3u_span - c_4u_span);

reaction_span = 1 - ((c_3m.^2 + c_3u_span.^2) - (c_4m.^2 + c_4u_span.^2))./(2*deltah_span);

%outputs
meanline.deltah_0 = deltah_0;
meanline.u = u;
meanline.c_3 = c_3;
meanline.c_4 = c_4;
meanline.c_m = c_m;
meanline.beta_3 = beta_3;
meanline.beta_4 = beta_4;
meanline.alpha_3 = alpha_3;
meanline.alpha_4 = alpha_4;
meanline.T_3 = T_3;
meanline.T_4 = T_4;
meanline.p_3 = p_3;
meanline.p_4 = p_4;
meanline.p_03 = p_03;
meanline.p_04 = p_04;
meanline.eta_tt = eta_tt;
meanline.M_3 = M_3;
meanline.M_4 = M_4;
meanline.M_3_r = M_3_r;
meanline.M_4_r = M_4_r;
meanline.H_2 = H_2;
meanline.H_3 = H_3;
meanline.H_4 = H_4;

span.r = r;
span.mean_radius_i = mean_radius_i;
span.beta_2_span = beta_2_span;
span.beta_3_span = beta_3_span;
span.beta_4_span = beta_4_span;
span.reaction = reaction_span;
end