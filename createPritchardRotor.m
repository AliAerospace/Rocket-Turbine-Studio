function pritchard = createPritchardRotor(inputs)

%11 independent inputs
radius = inputs.radius;
axial_chord = inputs.axial_chord;
tangential_chord = inputs.tangential_chord;
throat = inputs.throat;
unguided_turning = inputs.unguided_turning;
beta_in = inputs.beta_in;
epsilon_in = inputs.epsilon_in;
r_le = inputs.r_le;
beta_out = inputs.beta_out;
r_te = inputs.r_te;
N = inputs.N;

point_count = 201;

%pritchard eqn 1
pitch = 2*pi*radius/N;

%exit half wedge initial guess
epsilon_out = unguided_turning/2;

converged = false;
%exit half wedge iteration pritchard appendix F
for iteration = 1:200

    %P1 pritchard eqn 2 to 4
    beta_1 = beta_out - epsilon_out;
    x_1 = axial_chord - r_te*(1 + sind(beta_1));
    y_1 = r_te*cosd(beta_1);
    P1 = [x_1, y_1];

    %P2 pritchard eqn 5 to 7
    beta_2 = beta_out - epsilon_out + unguided_turning;
    x_2 = axial_chord - r_te + (throat + r_te)*sind(beta_2);
    y_2 = pitch - (throat + r_te)*cosd(beta_2);
    P2 = [x_2, y_2];

    %unguided turning circle appendix F
    x_0 = ((y_1 - y_2)*tand(beta_1)*tand(beta_2) +x_1*tand(beta_2) - x_2*tand(beta_1))/(tand(beta_2) - tand(beta_1));
    y_0 = -(x_0 - x_1)/tand(beta_1) + y_1;
    r_0 = sqrt((x_1 - x_0)^2 + (y_1 - y_0)^2);
    y_2_circle = y_0 + sqrt(r_0^2 - (x_2 - x_0)^2);

    %stop when P2 on circle
    if abs(y_2 - y_2_circle) < 1e-10 * axial_chord
        converged = true;
        break
    end
    %update half wedge
    epsilon_out = epsilon_out*(y_2/y_2_circle)^4;
end
assert(converged, 'Exit half wedge iteration did not converge.')

    %P3 pritchard eqn 8 to 10
    beta_3 = beta_in + epsilon_in;
    x_3 = r_le*(1 - sind(beta_3));
    y_3 = tangential_chord + r_le*cosd(beta_3);
    P3 = [x_3, y_3];

    

    %P4 pritchard eqn 11 to 13
    beta_4 = beta_in - epsilon_in;
    x_4 = r_le*(1 + sind(beta_4));
    y_4 = tangential_chord - r_le*cosd(beta_4);
    P4 = [x_4, y_4];

    %P5 prichard eqn 14 to 16
    beta_5 = beta_out + epsilon_out;
    x_5 = axial_chord - r_te*(1 - sind(beta_5));
    y_5 = -r_te*cosd(beta_5);
    P5 = [x_5, y_5];

    %suction cubic appendix C
    suction_d = (tand(beta_2) + tand(beta_3))/(x_2 - x_3)^2 - 2*(y_2 - y_3)/(x_2 - x_3)^3;

    suction_c = (y_2 - y_3)/(x_2 - x_3)^2 - tand(beta_3)/(x_2 - x_3) - suction_d*(x_2 + 2*x_3);

    suction_b = tand(beta_3) - 2*suction_c*x_3 - 3*suction_d*x_3^2;

    suction_a = y_3 - suction_b*x_3 - suction_c*x_3^2 - suction_d*x_3^3;

    %pressure cubic appendix C
    pressure_d = (tand(beta_5) + tand(beta_4))/(x_5 - x_4)^2 - 2*(y_5 - y_4)/(x_5 - x_4)^3;

    pressure_c = (y_5 - y_4)/(x_5 - x_4)^2 - tand(beta_4)/(x_5 - x_4) - pressure_d*(x_5 + 2*x_4);

    pressure_b = tand(beta_4) - 2*pressure_c*x_4 - 3*pressure_d*x_4^2;

    pressure_a = y_4 - pressure_b*x_4 - pressure_c*x_4^2 - pressure_d*x_4^3;

    %suction cubic P3 to P2
    suction_x = linspace(x_3, x_2, point_count).';
    suction_y = suction_a + suction_b*suction_x + suction_c*suction_x.^2 + suction_d*suction_x.^3;
    suction_curve = [suction_x, suction_y];

    %pressure cubic from P5 to P4
    pressure_x = linspace(x_5, x_4, point_count).';
    pressure_y = pressure_a + pressure_b*pressure_x + pressure_c*pressure_x.^2 + pressure_d*pressure_x.^3;
    pressure_curve = [pressure_x, pressure_y];

    %circle parameter
    q = linspace(0, 1, point_count).';

    %leading edge circle (drawing the circle) appendix B 
    leading_edge_centre = [r_le, tangential_chord];
    theta_3 = atan2(y_3 - leading_edge_centre(2), x_3 - leading_edge_centre(1));
    theta = pi + (theta_3 - pi)*q;
    leading_edge_upper = leading_edge_centre + r_le*[cos(theta), sin(theta)];

    theta_4 = atan2(y_4 - leading_edge_centre(2), x_4 - leading_edge_centre(1));
    theta = theta_4 + (-pi - theta_4)*q;
    leading_edge_lower = leading_edge_centre + r_le*[cos(theta), sin(theta)];

    %unguided turning circle appendix F
    unguided_centre = [x_0, y_0];
    theta_2 = atan2(y_2 - y_0, x_2 - x_0);
    theta_1 = atan2(y_1 - y_0, x_1 - x_0);
    theta = theta_2 + (theta_1 - theta_2)*q;
    unguided_curve = unguided_centre + r_0*[cos(theta), sin(theta)];

    %trailing edge circle appendix B
    trailing_edge_centre = [axial_chord - r_te, 0];
    theta_1 = atan2(y_1, x_1 - trailing_edge_centre(1));
    theta = theta_1 - theta_1*q;
    trailing_edge_upper = trailing_edge_centre + r_te*[cos(theta), sin(theta)];

    theta_5 = atan2(y_5, x_5 - trailing_edge_centre(1));
    theta = theta_5*q;
    trailing_edge_lower = trailing_edge_centre + r_te*[cos(theta), sin(theta)];

    %complete rotor profile
    profile = [leading_edge_upper; suction_curve(2:end, :); unguided_curve(2:end, :); trailing_edge_upper(2:end, :); trailing_edge_lower(2:end, :); pressure_curve(2:end, :); leading_edge_lower(2:end, :)];
    
    %outputs required later for fitting sudhofs rotor to this one
    pritchard.inputs = inputs;
    pritchard.pitch = pitch;
    pritchard.epsilon_out = epsilon_out;
    pritchard.beta_1 = beta_1;
    pritchard.beta_2 = beta_2;
    pritchard.beta_3 = beta_3;
    pritchard.beta_4 = beta_4;
    pritchard.beta_5 = beta_5;
    pritchard.P1 = P1;
    pritchard.P2 = P2;
    pritchard.P3 = P3;
    pritchard.P4 = P4;
    pritchard.P5 = P5;
    pritchard.points = [P1; P2; P3; P4; P5];
    pritchard.suction_coefficients = [suction_a, suction_b, suction_c, suction_d];
    pritchard.pressure_coefficients = [pressure_a, pressure_b, pressure_c, pressure_d];
    pritchard.suction_curve = suction_curve;
    pritchard.pressure_curve = pressure_curve;
    pritchard.unguided_curve = unguided_curve;
    pritchard.leading_edge_upper = leading_edge_upper;
    pritchard.leading_edge_lower = leading_edge_lower;
    pritchard.trailing_edge_upper = trailing_edge_upper;
    pritchard.trailing_edge_lower = trailing_edge_lower;
    pritchard.leading_edge_centre = leading_edge_centre;
    pritchard.trailing_edge_centre = trailing_edge_centre;
    pritchard.unguided_centre = unguided_centre;
    pritchard.unguided_radius = r_0;
    pritchard.profile = profile;
end

