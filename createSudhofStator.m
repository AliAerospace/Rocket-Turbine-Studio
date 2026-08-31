%%file createSudhofStator.m
function stator = createSudhofStator(inputs)
r_te = inputs.r_te;
chord = inputs.chord;
n_stagger = inputs.n_stagger;
beta_1 = inputs.beta_1;
wedge_angle_1 = inputs.wedge_angle_1;
beta_2 = inputs.beta_2;
gamma_no = inputs.gamma_no;
semi_minor_axis_le = inputs.semi_minor_axis_le;
eccentricity_le = inputs.eccentricity_le;
s_bar = inputs.s_bar;
u_bar = inputs.u_bar;
t_bar = inputs.t_bar;
w_45 = inputs.w_45;
Rc = inputs.Rc;
p_quartic = inputs.p_quartic;
b_quartic = inputs.b_quartic;
radius = inputs.radius;
N = inputs.N;
retraction = inputs.retraction;
contraction_ratio = inputs.contraction_ratio;

beta_2_suction = beta_2 - gamma_no/2;
beta_2_pressure = beta_2 + gamma_no/2;
suction_angle = beta_2_suction - 90;
pressure_angle = beta_2_pressure - 90;

%stagger
stagger_angle = calculateStaggerAngle(beta_1, beta_2, n_stagger);
axial_chord = chord*cosd(stagger_angle);
tangential_chord = chord*sind(stagger_angle);

te_tip = [axial_chord, 0];

q = linspace(0, 1, 101).';

%P1
P1 = te_tip - r_te*[cosd(suction_angle) sind(suction_angle)];

%P5
P5 = P1 + r_te*[cosd(pressure_angle - 90) sind(pressure_angle - 90)];

angle_1 = suction_angle + 180;
T1 = [cosd(angle_1) sind(angle_1)];

%P3 and P4 ellipse
theta_4 = beta_1 - 90 - wedge_angle_1;
a_le = semi_minor_axis_le/sqrt(1 - eccentricity_le^2);
lambda = 1/sqrt(cosd(theta_4)^2/a_le^2 + sind(theta_4)^2/semi_minor_axis_le^2);
phi_4 = mod(atan2(-lambda*cosd(theta_4)/a_le, lambda*sind(theta_4)/semi_minor_axis_le), 2*pi);

P3 = [0 tangential_chord];
P4 = [a_le + a_le*cos(phi_4), tangential_chord+semi_minor_axis_le*sin(phi_4)];

%P3 to P4 ellipse
phi = pi + (phi_4 - pi)*q;
P3_P4_curve = [a_le + a_le*cos(phi), tangential_chord + semi_minor_axis_le*sin(phi)];


%P4 to P5 conic
T4 = [cosd(theta_4) sind(theta_4)];
T5 = [cosd(pressure_angle) sind(pressure_angle)];
intersection = [T4(:) -T5(:)]\(P5 - P4).';
S1 = P4 + intersection(1)*T4;
P4_P5_curve = ((1 - q).^2*P4 + 2*w_45*q.*(1 - q)*S1 + q.^2*P5)./((1 - q).^2 + 2*w_45*q.*(1 - q) + q.^2);

%P2
pitch = 2*pi*radius/N;
t0 = pitch*sind(beta_2_suction);
throat = contraction_ratio*t0;

%arc length of pressure side P4 to P5
segment_length = vecnorm(diff(P4_P5_curve, 1, 1), 2, 2);
arc_length = [0; cumsum(segment_length)];

%retraction is from P5 towards P4
arc_target = (1 - retraction)*arc_length(end);

%move P2p on pressure side of the other blade
P2p_local = interp1(arc_length, P4_P5_curve, arc_target);
P2p = P2p_local + [0, pitch];

%P2p tangent on pressure side
pressure_dx = gradient(P4_P5_curve(:, 1), arc_length);
pressure_dy = gradient(P4_P5_curve(:, 2), arc_length);
pressure_tangent = interp1(arc_length, [pressure_dx, pressure_dy], arc_target);
pressure_tangent = pressure_tangent/norm(pressure_tangent);
throat_normal = [-pressure_tangent(2), pressure_tangent(1)];
P2 = P2p - throat*throat_normal;
T2 = -pressure_tangent;


%P1 to P2 non rational quartic
L = norm(P2 - P1);
G = L*sqrt(p_quartic*b_quartic);
H = L*sqrt(p_quartic*(1 - b_quartic));

Q0 = P1;
Q1 = P1 + G*T1;
Q3 = P2 - H*T2;
Q4 = P2;

n = 4;
k_0 = 0;
k_1 = 1/Rc;

Pax = Q1(1) - Q0(1);
Pay = Q1(2) - Q0(2);
Pbx = Q4(1) - Q3(1);
Pby = Q4(2) - Q3(2);

E1a = Pax*(Q0(2) - 2*Q1(2));
E2a = Pay*(Q0(1) - 2*Q1(1));
E1b = Pbx*(Q4(2) - 2*Q3(2));
E2b = Pby*(Q4(1) - 2*Q3(1));

La = G^3*n*k_0/(n - 1);
Lb = H^3*n*k_1/(n - 1);

Q2x = (Pbx*(La - E1a + E2a) - Pax*(Lb - E1b + E2b))/(Pax*Pby - Pay*Pbx);
Q2y = (Lb - E1b + E2b + Pby*Q2x)/Pbx;
Q2 = [Q2x Q2y];

P1_P2_curve = (1 - q).^4*Q0 + 4*(1 - q).^3.*q*Q1 + 6*(1 - q).^2.*q.^2*Q2+4*(1 - q).*q.^3*Q3 + q.^4*Q4;

%P2 to P3 paluszny
T3 = [0 -1];
curvature_P2 = 1/Rc;
curvature_P3 = a_le/semi_minor_axis_le^2;
intersection = [T2(:) -T3(:)]\(P3 - P2).';
M = P2 + intersection(1)*T2;
triangle_area = abs((M(1) - P2(1))*(P3(2) - P2(2)) - (M(2) - P2(2))*(P3(1) - P2(1)))/2;
K2 = curvature_P2*norm(M - P2)^3/(4*triangle_area);
K3 = curvature_P3*norm(P3-M)^3/(4*triangle_area);

X = -K3*t_bar*(K2*t_bar^2 - s_bar*u_bar)/(u_bar*(K3*t_bar^2 - s_bar*u_bar)- K3*t_bar*(K2*t_bar^2 - s_bar*u_bar));
Y = -K2*t_bar*(K3*t_bar^2 - s_bar*u_bar)/(s_bar*(K2*t_bar^2 - s_bar*u_bar)- K2*t_bar*(K3*t_bar^2 - s_bar*u_bar));

w1 = nthroot(X*Y^2/(K2^2*K3*(X - 1)^4*(Y - 1)^2), 3)/3;
w2 = nthroot(X^2*Y/(K2*K3^2*(X - 1)^2*(Y - 1)^4), 3)/3;

R0 = P2;
R1 = X*P2 + (1 - X)*M;
R2 = (1 - Y)*M + Y*P3;
R3 = P3;

B0 = (1 - q).^3;
B1 = 3*q.*(1 - q).^2;
B2 = 3*q.^2.*(1 - q);
B3 = q.^3;

P2_P3_curve = (B0*R0 + w1*B1*R1 + w2*B2*R2 + B3*R3)./(B0 + w1*B1 + w2*B2 + B3);

%P5 to tip curved
angle_te = linspace(pressure_angle - 90, suction_angle, 101).';
P5_tip_curve = P1 + r_te*[cosd(angle_te) sind(angle_te)];

%tip to P1 linear
tip_P1_curve = [te_tip;P1];


%outputs
stator.P1 = P1;
stator.P2 = P2;
stator.P3 = P3;
stator.P4 = P4;
stator.P5 = P5;
stator.P2p = P2p;
stator.throat = throat;
stator.t0 = t0;
stator.te_tip = te_tip;
stator.P1_P2_curve = P1_P2_curve;
stator.P2_P3_curve = P2_P3_curve;
stator.P3_P4_curve = P3_P4_curve;
stator.P4_P5_curve = P4_P5_curve;
stator.P5_tip_curve = P5_tip_curve;
stator.tip_P1_curve = tip_P1_curve;
stator.stagger_angle = stagger_angle;
stator.axial_chord = axial_chord;
stator.tangential_chord = tangential_chord;
stator.points = [P1; P2; P3; P4; P5];
stator.profile = [P1_P2_curve; P2_P3_curve(2:end, :); P3_P4_curve(2:end, :); P4_P5_curve(2:end, :); P5_tip_curve(2:end, :); tip_P1_curve(2:end, :)];
end