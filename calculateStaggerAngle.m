function stagger_angle = calculateStaggerAngle(beta_i, beta_o, n)

%eqn 3.6 sudhof
curve_parameter = @(turning_angle) ((turning_angle*pi/180)*(n + 1)*pi^(-n/2 - 0.5))^(1/(n + 1));

t_i = curve_parameter(90 - beta_i);
t_o = curve_parameter(90 - beta_o);

%eqn 3.5 sudhof
x_i_1 = integral(@(a) sin(a.^(n + 1))/(n + 1), 0, t_i);
x_i_2 = integral(@(a) cos(a.^(n + 1))/(n + 1), 0, t_i);
x_o_1 = integral(@(a) sin(a.^(n + 1))/(n + 1), 0, t_o);
x_o_2 = integral(@(a) cos(a.^(n + 1))/(n + 1), 0, t_o);

%mirorring inlet
x_i_1 = -x_i_1;

%chord from mirrored point
chord_1 = x_o_1 - x_i_1;
chord_2 = x_o_2 - x_i_2;

stagger_angle = atan2d(chord_2, chord_1);

end