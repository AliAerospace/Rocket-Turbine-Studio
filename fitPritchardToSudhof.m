function [sudhof_inputs, fit] = fitPritchardToSudhof(pritchard, t_ri, gamma_turning_ri)

%mapping parameters from pritchard 
radius = pritchard.inputs.radius;
axial_chord = pritchard.inputs.axial_chord;
tangential_chord = pritchard.inputs.tangential_chord;
throat = pritchard.inputs.throat;
unguided_turning = pritchard.inputs.unguided_turning;
beta_in = pritchard.inputs.beta_in;
epsilon_in = pritchard.inputs.epsilon_in;
r_le = pritchard.inputs.r_le;
beta_out = pritchard.inputs.beta_out;
r_te = pritchard.inputs.r_te;
N = pritchard.inputs.N;

%sudhof parameters
d_m = 2*radius;
gamma_ri = 2*epsilon_in;
beta_3 = 90 - beta_in;
semi_minor_axis_ri = r_le;
eccentricity_ri = 0;
gamma_turning_ro = unguided_turning;
gamma_ro = 2*pritchard.epsilon_out;
beta_4 = 90 + beta_out;
r_ro = r_te;
t_ro = throat;
channel_expansion_ratio = t_ri/t_ro;

%converting axial-tangential from pritchard to tangential-axial for sudhof
P1 = [pritchard.P1(2), pritchard.P1(1) - axial_chord + r_ro];
P2 = [pritchard.P2(2), pritchard.P2(1) - axial_chord + r_ro];
P4 = [pritchard.P3(2), pritchard.P3(1) - axial_chord + r_ro];
P5 = [pritchard.P4(2), pritchard.P4(1) - axial_chord + r_ro];
P6 = [pritchard.P5(2), pritchard.P5(1) - axial_chord + r_ro];

%curve targets
quartic_target = flipud([pritchard.unguided_curve(:, 2), pritchard.unguided_curve(:, 1) - axial_chord + r_ro]);

suction_target = flipud([pritchard.suction_curve(:, 2), pritchard.suction_curve(:, 1) - axial_chord + r_ro]);

pressure_target = flipud([pritchard.pressure_curve(:, 2), pritchard.pressure_curve(:, 1) - axial_chord + r_ro]);

%P3
x_ro = 0;
z_ro = 0;

a_ri = semi_minor_axis_ri/sqrt(1 - eccentricity_ri^2);

x_ri = x_ro + tangential_chord;
z_ri = z_ro + r_ro - axial_chord + semi_minor_axis_ri;

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
B = 2*((P3_focus_x - x_ri)*P3_x/a_ri^2 + (P3_focus_z - z_ri)*P3_z/semi_minor_axis_ri^2);
C = (P3_focus_x - x_ri)^2/a_ri^2 + (P3_focus_z - z_ri)^2/semi_minor_axis_ri^2 - 1;
ellipse_distance = (-B + sqrt(B^2 - 4*A*C))/(2*A);

P3_x = P3_focus_x + (ellipse_distance + t_ri)*P3_x;
P3_z = P3_focus_z + (ellipse_distance + t_ri)*P3_z;
P3 = [P3_x, P3_z];

%fitting P1 to P2 quartic
q = linspace(0, 1 ,size(quartic_target, 1)).';
%endpoint tangents
theta_6 = atan2d(P6(2), P6(1));
theta_1 = atan2d(P1(2), P1(1));
delta_theta_61 = mod(theta_1 - theta_6 + 180, 360) - 180;

D0 = sign(delta_theta_61)*[-sind(theta_1),cosd(theta_1)];

P2_angle = beta_4 - gamma_ro/2 + gamma_turning_ro;
D4 = [cosd(P2_angle), -sind(P2_angle)];

B0 = (1 - q).^4;
B1 = 4*(1 - q).^3.*q;
B2 = 6*(1 - q).^2.*q.^2;
B3 = 4*(1 - q).*q.^3;
B4 = q.^4;

%G H and Q2 solved using linear least squares
constant = (B0 + B1)*P1 + (B3 + B4)*P2;

matrix_x = [B1*D0(1), -B3*D4(1), B2, zeros(size(q))];
matrix_z = [B1*D0(2), -B3*D4(2), zeros(size(q)),  B2];

matrix = [matrix_x; matrix_z];
right_side = [quartic_target(:, 1) - constant(:, 1); quartic_target(:, 2) - constant(:, 2)];

solution = matrix\right_side;

G = solution(1);
H = solution(2);
Q0 = P1;
Q1 = P1 + G*D0;
Q2 = [solution(3), solution(4)];
Q3 = P2 - H*D4;
Q4 = P2;

chord_length = norm(P2 - P1);
g = G/chord_length;
h = H/chord_length;

quartic_pointedness = g^2 + h^2;
quartic_bias = g^2/quartic_pointedness;

%curvature P1
cross_1 = (Q1(1) - Q0(1))*(Q0(2) - 2*Q1(2) + Q2(2)) - (Q1(2) - Q0(2))*(Q0(1) - 2*Q1(1) + Q2(1));

curvature_P1 = 3*cross_1/(4*G^3);

%curvature P2
cross_2 = (Q4(1) - Q3(1))*(Q4(2) - 2*Q3(2) + Q2(2)) - (Q4(2) - Q3(2))*(Q4(1) - 2*Q3(1) + Q2(1));

curvature_P2 = 3*cross_2/(4*H^3);

%fitted quartic
P1_P2_curve = B0*Q0 + B1*Q1 + B2*Q2 + B3*Q3 + B4*Q4;

%P3 to P4 conic
q = linspace(0, 1, 101).';
P3_angle = beta_3 + gamma_turning_ri;
P4_angle = beta_3 - gamma_ri/2;

P3_tangent_x = cosd(P3_angle);
P3_tangent_z = sind(P3_angle);
P4_tangent_x = cosd(P4_angle);
P4_tangent_z = sind(P4_angle);
P3_to_P4 = [P4(1) - P3(1); P4(2) - P3(2)];

tangent_lengths = [P3_tangent_x, -P4_tangent_x; P3_tangent_z, -P4_tangent_z]\P3_to_P4;

S1_x = P3(1) + tangent_lengths(1)*P3_tangent_x;
S1_z = P3(2) + tangent_lengths(1)*P3_tangent_z;
S1 = [S1_x,S1_z];

%conic weight. 1 for parabola
w = 1;

%bezier eqn
P3_P4_curve = ((1 - q).^2*P3 + 2*w*q.*(1 - q)*S1 + q.^2*P4)./((1 - q).^2 + 2*w*q.*(1 - q) + q.^2);

%curvature at P3 required by the paluszny curve
conic_first_derivative = 2*w*(S1 - P3);
conic_second_derivative = 2*(P3 - 2*w*S1 +P4) - 4*(1 - w)*P3 - 4*(w - 1)*2*w*(S1 - P3);
curvature_P3 = abs(conic_first_derivative(1)*conic_second_derivative(2) - conic_first_derivative(2)*conic_second_derivative(1))/norm(conic_first_derivative)^3;

%splitting pritchards P1 to P2 to allow fitting of sudhofs P2 to P3 and P3 to P4
[~, split_i] = min(sum((suction_target - P3).^2, 2));

split_i = max(2, min(split_i, size(suction_target, 1) - 1));
target_P2_P3 = resampleCurve(suction_target(1:split_i, :), q);
target_P3_P4 = resampleCurve(suction_target(split_i:end, :), q);

%fit P2 to P3 

P2_tangent = D4;
P3_tangent = [P3_tangent_x, P3_tangent_z];

s_values = logspace(-4, 1, 101);
u_values = [-logspace(-2, 6, 161), logspace(-2, 2, 81)];

best_J = inf;
s_bar = NaN;
u_bar = NaN;
P2_P3_curve = [];

P3_P4_direction = P3_P4_curve(2, :) - P3_P4_curve(1, :);

for s = s_values
    for u = u_values

        curve = evaluatePalusznyCurve(P2, P3, P2_tangent, P3_tangent, curvature_P2, curvature_P3, s, u, q);

        if isreal(curve) && all(isfinite(curve(:)))

            P2_P3_direction = curve(end, :) - curve(end - 1, :);

            tangent_cosine = dot(P2_P3_direction, P3_P4_direction)/(norm(P2_P3_direction)*norm(P3_P4_direction));

            if tangent_cosine > 0.9999

                J = norm(curve - target_P2_P3, 'fro');

                if J < best_J
                    best_J = J;
                    s_bar = s;
                    u_bar = u;
                    P2_P3_curve = curve;
                end
            end
        end
    end
end

%fit P5 to P6
q = linspace(0, 1, size(pressure_target, 1)).';

%P5 tangent
beta_3_suction = beta_3 - gamma_ri/2;
beta_3_pressure = beta_3 + gamma_ri/2;

theta_4 = atan2d(-semi_minor_axis_ri*cosd(beta_3_suction),a_ri*sind(beta_3_suction));
theta_5 = atan2d(-semi_minor_axis_ri*cosd(beta_3_pressure),a_ri*sind(beta_3_pressure)) + 180;

delta_theta_45 = mod(theta_5 - theta_4 + 180,360) - 180;
ellipse_direction = sign(delta_theta_45);

P5_tangent = ellipse_direction*[-a_ri*sind(theta_5), semi_minor_axis_ri*cosd(theta_5)];
P5_tangent = P5_tangent/norm(P5_tangent);

%P6 tangent
theta_6 = atan2d(P6(2), P6(1));
theta_1 = atan2d(P1(2), P1(1));

delta_theta_61 = mod(theta_1 - theta_6 + 180, 360) - 180;

P6_tangent = sign(delta_theta_61)*[-sind(theta_6), cosd(theta_6)];

B0 = (1 - q).^3;
B1 = 3*q.*(1 - q).^2;
B2 = 3*q.^2.*(1 - q);
B3 = q.^3;

%lambda_5 and lambda_6 findign with linear least squares
constant = (B0 + B1)*P5 + (B2 + B3)*P6;

matrix_x = [B1*P5_tangent(1), -B2*P6_tangent(1)];
matrix_z = [B1*P5_tangent(2), -B2*P6_tangent(2)];

matrix = [matrix_x; matrix_z];
right_side = [pressure_target(:, 1) - constant(:, 1); pressure_target(:, 2) - constant(:, 2)];

solution = matrix\right_side;

lambda_5 = solution(1);
lambda_6 = solution(2);

V0 = P5;
V1 = P5 + lambda_5*P5_tangent;
V2 = P6 - lambda_6*P6_tangent;
V3 = P6;

P5_P6_curve = (1 - q).^3*V0 + 3*q.*(1 - q).^2*V1 + 3*q.^2.*(1 - q)*V2 + q.^3*V3;

%inputs for createSudhofRotor function
sudhof_inputs.d_m = d_m;
sudhof_inputs.chord = hypot(axial_chord, tangential_chord);
sudhof_inputs.stagger_angle = atan2d(tangential_chord, axial_chord);
sudhof_inputs.gamma_turning_ri = gamma_turning_ri;
sudhof_inputs.gamma_ri = gamma_ri;
sudhof_inputs.beta_3 = beta_3;
sudhof_inputs.semi_minor_axis_ri = semi_minor_axis_ri;
sudhof_inputs.eccentricity_ri = eccentricity_ri;
sudhof_inputs.gamma_turning_ro = gamma_turning_ro;
sudhof_inputs.gamma_ro = gamma_ro;
sudhof_inputs.beta_4 = beta_4;
sudhof_inputs.r_ro = r_ro;
sudhof_inputs.t_ri = t_ri;
sudhof_inputs.curvature_P1 = curvature_P1;
sudhof_inputs.curvature_P2 = curvature_P2;
sudhof_inputs.quartic_pointedness = quartic_pointedness;
sudhof_inputs.quartic_bias = quartic_bias;
sudhof_inputs.lambda_5 = lambda_5;
sudhof_inputs.lambda_6 = lambda_6;
sudhof_inputs.s_bar = s_bar;
sudhof_inputs.u_bar = u_bar;
sudhof_inputs.channel_expansion_ratio = channel_expansion_ratio;
sudhof_inputs.N = N;

fit.P1 = P1;
fit.P2 = P2;
fit.P3 = P3;
fit.P4 = P4;
fit.P5 = P5;
fit.P6 = P6;
fit.points = [P1; P2; P3; P4; P5; P6];
fit.P1_P2_curve = P1_P2_curve;
fit.P2_P3_curve = P2_P3_curve;
fit.P3_P4_curve = P3_P4_curve;
fit.P5_P6_curve = P5_P6_curve;
fit.quartic_target = quartic_target;
fit.target_P2_P3 = target_P2_P3;
fit.target_P3_P4 = target_P3_P4;
fit.pressure_target = pressure_target;
fit.J = best_J;
end

function P2_P3_curve = evaluatePalusznyCurve(P2, P3, P2_tangent, P3_tangent, curvature_P2, curvature_P3, s_bar, u_bar, q)
t_bar = 1;

%intersection of tangents at the end points
tangent_lengths = [P2_tangent(:),-P3_tangent(:)]\(P3 - P2).';
M = P2 + tangent_lengths(1)*P2_tangent;

triangle_area = abs((M(1) - P2(1))*(P3(2) - P2(2)) - (M(2) - P2(2))*(P3(1) - P2(1)))/2;

%curvature parameters
K2 = abs(curvature_P2*norm(M - P2)^3/(4*triangle_area));
K3 = abs(curvature_P3*norm(P3 - M)^3/(4*triangle_area));

X = -K3*t_bar*(K2*t_bar^2 - s_bar*u_bar)/(u_bar*(K3*t_bar^2 - s_bar*u_bar) - K3*t_bar*(K2*t_bar^2 - s_bar*u_bar));
Y = -K2*t_bar*(K3*t_bar^2 - s_bar*u_bar)/(s_bar*(K2*t_bar^2 - s_bar*u_bar) - K2*t_bar*(K3*t_bar^2 - s_bar*u_bar));

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
end

function resample_curve = resampleCurve(curve, q)

distance = [0; cumsum(sqrt(sum(diff(curve ).^2, 2)))];

sample_distance = distance(end)*q;

resampled_x = interp1(distance, curve(:, 1), sample_distance);
resampled_z = interp1(distance, curve(:, 2), sample_distance);

resample_curve = [resampled_x, resampled_z];

end