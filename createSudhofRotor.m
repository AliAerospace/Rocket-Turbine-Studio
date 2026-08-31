function rotor = createSudhofRotor(inputs)

%inputs
d_m = inputs.d_m;
chord = inputs.chord;
stagger_angle = inputs.stagger_angle;
axial_chord = chord*cosd(stagger_angle);
tangential_chord = chord*sind(stagger_angle);
gamma_turning_ri = inputs.gamma_turning_ri;
gamma_ri = inputs.gamma_ri;
beta_3 = inputs.beta_3;
semi_minor_axis_ri = inputs.semi_minor_axis_ri;
eccentricity_ri = inputs.eccentricity_ri;
gamma_turning_ro = inputs.gamma_turning_ro;
gamma_ro = inputs.gamma_ro;
beta_4 = inputs.beta_4;
r_ro = inputs.r_ro;
t_ri = inputs.t_ri;
curvature_P1 = inputs.curvature_P1;
curvature_P2 = inputs.curvature_P2;
quartic_pointedness = inputs.quartic_pointedness;
quartic_bias = inputs.quartic_bias;
lambda_5 = inputs.lambda_5;
lambda_6 = inputs.lambda_6;
s_bar = inputs.s_bar;
u_bar = inputs.u_bar;
channel_expansion_ratio = inputs.channel_expansion_ratio;
N = inputs.N;

pitch = pi*d_m/N;

%generating all points P

%finding P1 and P6
x_ro = 0;
z_ro = 0;

%trailing edge tangent angles
beta_4_suction = beta_4 - gamma_ro/2;
beta_4_pressure = beta_4 + gamma_ro/2;

%P1
%P1 suction trailing edge location of tangent
P1_x = x_ro + r_ro*sind(beta_4_suction);
P1_z = z_ro + r_ro*cosd(beta_4_suction);
P1 = [P1_x, P1_z];

%P6
%P6 pressure trailing edge location of tangent
P6_x = x_ro - r_ro*sind(beta_4_pressure);
P6_z = z_ro - r_ro*cosd(beta_4_pressure);
P6 = [P6_x, P6_z];

%finding P4 and P5
a_ri = semi_minor_axis_ri/sqrt(1 - eccentricity_ri^2);

x_ri = x_ro + tangential_chord;
z_ri = z_ro + r_ro - axial_chord + semi_minor_axis_ri;

%leading edge tangent angles
beta_3_suction = beta_3 - gamma_ri/2;
beta_3_pressure = beta_3 + gamma_ri/2;

theta_4 = atan2d(-semi_minor_axis_ri*cosd(beta_3_suction), a_ri*sind(beta_3_suction));
theta_5 = atan2d(-semi_minor_axis_ri*cosd(beta_3_pressure), a_ri*sind(beta_3_pressure)) + 180;

%P4
%P4 suction leading edge tangent point
P4_x = x_ri + a_ri*cosd(theta_4);
P4_z = z_ri + semi_minor_axis_ri*sind(theta_4);
P4 = [P4_x, P4_z];

%P5
%P5 pressure leading edge tangent point
P5_x = x_ri + a_ri*cosd(theta_5);
P5_z = z_ri + semi_minor_axis_ri*sind(theta_5);
P5 = [P5_x, P5_z];

%finding P2
t_ro = t_ri/channel_expansion_ratio;

P2_angle = beta_4 - gamma_ro/2 + gamma_turning_ro;

adjacent_trailing_edge_x = x_ro + pitch;
adjacent_trailing_edge_z = z_ro;

%P2
P2_x = adjacent_trailing_edge_x - (r_ro + t_ro)*sind(P2_angle);
P2_z = adjacent_trailing_edge_z - (r_ro + t_ro)*cosd(P2_angle);
P2 = [P2_x, P2_z];

%finding P3
P3_angle = beta_3 + gamma_turning_ri;

focal_distance_ri = a_ri*eccentricity_ri;

if cosd(P3_angle) >= 0
    P3_focus_x = x_ri + focal_distance_ri;
else
    P3_focus_x = x_ri - focal_distance_ri;
end

P3_focus_z = z_ri;
P3_x = cosd(P3_angle);
P3_z = sind(P3_angle);

A = P3_x^2/a_ri^2 + P3_z^2/semi_minor_axis_ri^2;
B = 2*((P3_focus_x  - x_ri)*P3_x/a_ri^2 + (P3_focus_z - z_ri)*P3_z/semi_minor_axis_ri^2);
C = (P3_focus_x - x_ri)^2/a_ri^2 + (P3_focus_z - z_ri)^2/semi_minor_axis_ri^2 - 1;
ellipse_distance = (-B + sqrt(B^2 - 4*A*C))/(2*A);

%P3 
P3_x = P3_focus_x + (ellipse_distance + t_ri)*P3_x;
P3_z = P3_focus_z + (ellipse_distance + t_ri)*P3_z;
P3 = [P3_x, P3_z];


%connecting all points P with curves based on sudhofs proposed methodology
q = linspace(0, 1, 101).';
%P1 to P6 circle 
theta_6 = atan2d(P6_z - z_ro, P6_x - x_ro);
theta_1 = atan2d(P1_z - z_ro, P1_x - x_ro);

delta_theta_61 = mod(theta_1 - theta_6 + 180, 360) - 180;
theta_61 = theta_6 + q*delta_theta_61;

P6_P1_x = x_ro + r_ro*cosd(theta_61);
P6_P1_z = z_ro + r_ro*sind(theta_61);
P6_P1_curve = [P6_P1_x, P6_P1_z];

%P1 to P2 non rational qaurtic (sudhof appendix B.3)
Q0 = P1;
Q4 = P2;

P1_tangent_x = sign(delta_theta_61)*(-sind(theta_1));
P1_tangent_z = sign(delta_theta_61)*( cosd(theta_1));

D0 = [P1_tangent_x, P1_tangent_z];
D4 = [cosd(P2_angle), -sind(P2_angle)];

g = sqrt(quartic_pointedness*quartic_bias);
h = sqrt(quartic_pointedness*(1 - quartic_bias));
G = g*norm(Q4 - Q0);
H = h*norm(Q4 - Q0);

Q1 = Q0 + G*D0;
Q3 = Q4 - H*D4;

Pax = Q1(1) - Q0(1);
Paz = Q1(2) - Q0(2);
Pbx = Q4(1) - Q3(1);
Pbz = Q4(2) - Q3(2);

E1a = Pax*(Q0(2) - 2*Q1(2));
E2a = Paz*(Q0(1) - 2*Q1(1));
E1b = Pbx*(Q4(2) - 2*Q3(2));
E2b = Pbz*(Q4(1) - 2*Q3(1));

La = G^3*4*curvature_P1/3;
Lb = H^3*4*curvature_P2/3;

Q2_x = (Pbx*(La - E1a + E2a) - Pax*(Lb - E1b + E2b))/(Pax*Pbz - Paz*Pbx);
Q2_z = (Lb - E1b + E2b + Pbz*Q2_x)/Pbx;
Q2 = [Q2_x, Q2_z];

%bezier
P1_P2_curve = (1 - q).^4*Q0 + 4*(1 - q).^3.*q*Q1+ 6*(1 - q).^2.*q.^2*Q2 + 4*(1 - q).*q.^3*Q3 + q.^4*Q4;


%P3 to P4 conic section
P3_angle = beta_3 + gamma_turning_ri;
P4_angle = beta_3 - gamma_ri/2;

P3_tangent_x = cosd(P3_angle);
P3_tangent_z = sind(P3_angle);
P4_tangent_x = cosd(P4_angle);
P4_tangent_z = sind(P4_angle);
P3_to_P4 = [P4_x - P3_x; P4_z - P3_z];

tangent_lengths = [P3_tangent_x, -P4_tangent_x; P3_tangent_z,-P4_tangent_z]\P3_to_P4;
S1_x = P3_x + tangent_lengths(1)*P3_tangent_x;
S1_z = P3_z + tangent_lengths(1)*P3_tangent_z;
S1 = [S1_x, S1_z];

%conic weight. 1 for parabola
w = 1;

%bezier eqn
P3_P4_curve = ((1 - q).^2*P3 + 2*w*q.*(1 - q)*S1 + q.^2*P4)./((1 - q).^2 + 2*w*q.*(1 - q) + q.^2);

%P4 to P5 ellipse
delta_theta_45 = mod(theta_5 - theta_4 + 180, 360) - 180;
theta_45 = theta_4 + q*delta_theta_45;

P4_P5_x = x_ri + a_ri * cosd(theta_45);
P4_P5_z = z_ri + semi_minor_axis_ri*sind(theta_45);
P4_P5_curve = [P4_P5_x, P4_P5_z];

%P5 to P6 non rational cubic
ellipse_direction = sign(delta_theta_45);

P5_tangent_length = sqrt((-a_ri*sind(theta_5))^2 + (semi_minor_axis_ri*cosd(theta_5))^2);

P5_tangent_x = ellipse_direction*(-a_ri*sind(theta_5))/P5_tangent_length;
P5_tangent_z = ellipse_direction*(semi_minor_axis_ri*cosd(theta_5))/P5_tangent_length;

P6_tangent_x = sign(delta_theta_61)*(-sind(theta_6));
P6_tangent_z = sign(delta_theta_61)*( cosd(theta_6));

V0 = P5;
V1_x = P5_x + lambda_5*P5_tangent_x;
V1_z = P5_z + lambda_5*P5_tangent_z;
V1 = [V1_x, V1_z];

V2_x = P6_x - lambda_6*P6_tangent_x;
V2_z = P6_z - lambda_6*P6_tangent_z;
V2 = [V2_x, V2_z];
V3 = P6;
%bezier eqn
P5_P6_curve = (1 - q).^3*V0 + 3*q.*(1 - q).^2*V1 + 3*q.^2.*(1 - q)*V2 + q.^3*V3;

%P2 to P3 paluszny rational cubic
t_bar = 1;

%P3 curvature from P3 to P4 rational conic
conic_first_derivative = 2*w*(S1 - P3);
conic_second_derivative = 2*(P3 - 2*w*S1 + P4) - 4*(1 - w)*P3 - 4*(w - 1)*2*w*(S1 - P3);
curvature_P3 = abs(conic_first_derivative(1)*conic_second_derivative(2) - conic_first_derivative(2)*conic_second_derivative(1))/ norm(conic_first_derivative)^3;

%tangent intersection
T2_line = D4;
T3_line = [P3_tangent_x, P3_tangent_z];

tangent_lengths_23 = [T2_line(:), -T3_line(:)]\(P3 - P2).';

M = P2 + tangent_lengths_23(1)*T2_line;

triangle_area = abs((M(1) - P2(1))*(P3(2) - P2(2))- (M(2) - P2(2))*(P3(1) - P2(1))) / 2;

%curvature parameters
K2 = abs(curvature_P2*norm(M - P2)^3/(4*triangle_area));
K3 = abs(curvature_P3*norm(P3 - M)^3/(4*triangle_area));

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


%outputs
rotor.P1 = P1;
rotor.P2 = P2;
rotor.P3 = P3;
rotor.P4 = P4;
rotor.P5 = P5;
rotor.P6 = P6;
rotor.P1_P2_curve = P1_P2_curve;
rotor.P2_P3_curve = P2_P3_curve;
rotor.P3_P4_curve = P3_P4_curve;
rotor.P4_P5_curve = P4_P5_curve;
rotor.P5_P6_curve = P5_P6_curve;
rotor.P6_P1_curve = P6_P1_curve;

rotor.points = [P1; P2; P3; P4; P5; P6];
end