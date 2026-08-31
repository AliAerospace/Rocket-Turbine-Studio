clear
clc
close all

%AliAerospace https://github.com/AliAerospace/rocket-turbine-studio/

%please read the instruction outlined in the github before using this
%program. As the methodologies here are novel, the github provides
%clarification.

%to generate CAD .iges file (which can be converted to .STEP in most
%programs) set to true
generate_IGES = true;
IGES_folder = fullfile(pwd, 'CAD output');

%name the .iges file outputs for rotor and stator
rotor_name = 'rotor';
stator_name = 'stator';


%meanline inputs
cycle.P = 3.076e6; %shaft power W
cycle.mdot = 2.7; %mass flow rate kg/s
cycle.d_m = 323; %mean diameter mm
cycle.omega = 400; %rotational speed rev/s
cycle.reaction = 0; %reaction at mean stage
cycle.beta_1 = 65.1; %stator inlet blade angle deg
cycle.beta_2 = 15; %stator outlet blade angle deg
cycle.delta_3 = 0; %rotor inlet flow deviation deg
cycle.gamma_ri = 8; %rotor inlet wedge angle deg
cycle.T_01 = 875; %inlet total temperature K
cycle.c_p = 8.1e3; %specific heat capacity J/kgK
cycle.gamma = 1.37; %ratio of specific heats
cycle.Mmol = 3.83; %molar mass kg/kmol
cycle.p_01 = 2.623134767e6; %inlet total pressure Pa
cycle.degree_of_admission = 1; %flow admitted to fraction of annulus
cycle.eta_m = 0.99; %mechanical efficiency
cycle.eta_n = 0.93; %nozzle efficiency
cycle.n = 1; %inlet spanwise vortex exponent
cycle.m = 1; %work distribution exponent

%rotor inputs
rotorValues.chord = 30; %rotor chord mm
rotorValues.N = 37; %number of rotor blades
rotorValues.unguided_turning = 8; %unguided turning angle deg
rotorValues.r_le_ratio = 0.02; %leading edge ratio of radius to axial chord
rotorValues.r_te_ratio = 0.01; %trailing edge ratio of radius to axial chord
rotorValues.t_ri = 9; %distance between rotors at inlet mm 
rotorValues.gamma_turning_ri = -10; %rotor inlet turning angle deg
rotorValues.n_stagger = 2; %stagger calculation exponent

%stator inputs
statorValues.r_te = 1.3712; %trailing edge radius mm
statorValues.wedge_angle_1 = 6.4329; %stator inlet wedge angle deg
statorValues.gamma_no = 0.4393; %stator outlet wedge angle deg
statorValues.semi_minor_axis_le = 6.9402; %leading edge ellipse semi minor axis mm
statorValues.eccentricity_le = 0.2347; %leading edge ellipse eccentricity
statorValues.w_45 = 0.4; %P4 to P5 rational conic weight
statorValues.Rc = 130; %radius of curvature at P2 mm
statorValues.p_quartic = 0.067355; %P1 to P2 quartic pointedness
statorValues.b_quartic = 0.733924; %P1 to P2 quartic bias
statorValues.angle_2 = 115; %tangent angle at P2 deg
statorValues.s_bar = -2.5; %paluszny curve s parameter
statorValues.t_bar = 0.8; %paluszny curve t parameter
statorValues.u_bar = -2.5; %paluszny curve u parameter
statorValues.N = 37; %number of stator blades
statorValues.retraction = 0.324746; %throat retraction fraction (P4 to P5 distance fraction at which the throat is located on the pressure side P2p)
statorValues.contraction_ratio = 0.95; %throat to reference opening ratio (t/t_0)
statorValues.chord = 41.2594; %blade chord mm
statorValues.n_stagger = 2; %stagger calculation exponent








%spanwise and meanline
[meanline, span] = solveMeanline(cycle);

%generate rotor and stator geometries
rotorGeometry = generateRotorGeometry(rotorValues, cycle, span);
statorGeometry = generateStatorGeometry(statorValues, cycle, span);


%cad generation
if generate_IGES
exportBladeRowsIGES(rotorGeometry, statorGeometry, rotorValues, statorValues, cycle, rotor_name, stator_name, IGES_folder);
end


%sudhof rotor results
disp('Root Sudhof inputs')
disp(rotorGeometry.section_sudhof_inputs{1})

disp('Mean Sudhof inputs')
disp(rotorGeometry.section_sudhof_inputs{2})

disp('Tip Sudhof inputs')
disp(rotorGeometry.section_sudhof_inputs{3})

%results
disp('Meanline outputs')
disp(meanline)

fprintf('Rotor root radius: %.3f mm\n', rotorGeometry.section_radius(1))
fprintf('Rotor mean radius: %.3f mm\n', rotorGeometry.section_radius(2))
fprintf('Rotor tip radius:  %.3f mm\n', rotorGeometry.section_radius(3))

%plot of angle variations from root to tip

figure('Color', 'w')
plot(1000*span.r, span.beta_2_span, 'LineWidth', 1.5)
hold on
plot(1000*span.r, span.beta_3_span, 'LineWidth', 1.5)
plot(1000*span.r, span.beta_4_span, 'LineWidth', 1.5)
grid on
xlabel('Radius [mm]')
ylabel('\beta [deg]')
legend('\beta_2', '\beta_3', '\beta_4')

%plot of reaction
figure('Color', 'w')
plot(1000*span.r, span.reaction, 'k', 'LineWidth', 1.5)
hold on
yline(mean(span.reaction), 'r--', 'Mean')
grid on
xlabel('Radius [mm]')
ylabel('Reaction')

%plot rotor
section_name = {'Root', 'Mean', 'Tip'};
p_colour = lines(6);

figure('Color', 'w')
r_layout = tiledlayout(1, 3);

for k = 1:3
    rotor = rotorGeometry.section_rotor{k};
    r_profile = [rotor.P1_P2_curve; rotor.P2_P3_curve(2:end, :); rotor.P3_P4_curve(2:end, :); rotor.P4_P5_curve(2:end, :); rotor.P5_P6_curve(2:end, :); rotor.P6_P1_curve(2:end, :)];
    ax = nexttile(r_layout);
    plot(ax, r_profile(:, 2), r_profile(:, 1), 'k', 'LineWidth', 1.4)
    hold on
    pritchard = rotorGeometry.section_pritchard{k};
    axial_chord = [pritchard.leading_edge_centre(1), pritchard.trailing_edge_centre(1)] - pritchard.inputs.axial_chord + pritchard.inputs.r_te;
    chord_tangential = [pritchard.leading_edge_centre(2), pritchard.trailing_edge_centre(2)];
    plot(ax, axial_chord, chord_tangential, 'r', 'LineWidth', 1.2)

    for point = 1:6
        plot(ax, rotor.points(point, 2), rotor.points(point, 1), 'o', 'MarkerFaceColor', p_colour(point, :), 'MarkerEdgeColor', p_colour(point, :))
    end

    axis equal
    grid on
    xlabel('axial z [mm]')
    ylabel('tangential x [mm]')
    title(sprintf('%s stagger angle of %.2f deg', section_name{k}, rotorGeometry.section_stagger(k)))
end

legend(r_layout.Children(end), 'Sudhof rotor', 'Chord', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'Location', 'best')

%plot stator

sp_colour = lines(5);
figure('Color', 'w')
s_layout = tiledlayout(1, 3);

for k = 1:3
    stator = statorGeometry.sections{k};
    ax = nexttile(s_layout);
    plot(stator.profile(:, 1), stator.profile(:, 2), 'k', 'LineWidth', 1.4)
    hold on
    plot([stator.P3(1), stator.te_tip(1)], [stator.P3(2), stator.te_tip(2)], 'r', 'LineWidth', 1.2)

    for point = 1:5
        plot(stator.points(point, 1), stator.points(point, 2), 'o', 'MarkerFaceColor', sp_colour(point, :), 'MarkerEdgeColor', sp_colour(point, :))
    end

    axis equal
    grid on
    xlabel('axial z [mm]')
    ylabel('tangential x [mm]')
    title(sprintf('%s stagger angle of %.2f deg', section_name{k}, statorGeometry.stagger(k)))
end
legend(s_layout.Children(end), 'Stator', 'Chord', 'P1', 'P2', 'P3', 'P4', 'P5', 'Location', 'best')

%sudhof fitted rotor versus pritchard target rotor at mean

pritchard = rotorGeometry.section_pritchard{2};
rotor = rotorGeometry.section_rotor{2};
r_profile = [rotor.P1_P2_curve; rotor.P2_P3_curve(2:end, :); rotor.P3_P4_curve(2:end, :); rotor.P4_P5_curve(2:end, :); rotor.P5_P6_curve(2:end, :); rotor.P6_P1_curve(2:end, :)];
pritchard_axial = pritchard.profile(:, 1) - pritchard.inputs.axial_chord + pritchard.inputs.r_te;
pritchard_tangential = pritchard.profile(:, 2);

figure('Color', 'w')
plot(pritchard_axial, pritchard_tangential, 'r--', 'LineWidth', 1.2)
hold on
plot(r_profile(:, 2), r_profile(:, 1), 'k', 'LineWidth', 1.5)

for point = 1:6
    plot(rotor.points(point, 2), rotor.points(point, 1), 'o', 'MarkerFaceColor', p_colour(point, :), 'MarkerEdgeColor', p_colour(point, :))
end

axis equal
grid on
xlabel('axial z [mm]')
ylabel('tangential x [mm]')
title('Sudhof rotor parameters fitted to Pritchards')
legend('Pritchard rotor', 'Sudhof rotor', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'Location', 'northeast')

%plot of cylindrical mapped rotor and stator

figure('Color', 'w')
wrap_layout = tiledlayout(1, 2);
section_color = lines(3);
r_ax = nexttile(wrap_layout);
hold on

for segment = 1:6
    surf(r_ax, rotorGeometry.mantle{segment}.X, rotorGeometry.mantle{segment}.Y, rotorGeometry.mantle{segment}.Z, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.75)
end

for k = 1:3
    radius = rotorGeometry.section_radius(k);
    rotor_outline = [rotorGeometry.section_curve{k, 1}; rotorGeometry.section_curve{k, 2}; rotorGeometry.section_curve{k, 3}; rotorGeometry.section_curve{k, 4}; rotorGeometry.section_curve{k, 5}; rotorGeometry.section_curve{k, 6}];
    theta_canvas = linspace(min(rotor_outline(:, 1))/radius - 0.03, max(rotor_outline(:, 1))/radius + 0.03, 30);
    z_canvas = linspace(min(rotor_outline(:, 2)), max(rotor_outline(:, 2)), 2);
    [theta_canvas, z_canvas] = meshgrid(theta_canvas, z_canvas);
    surf(r_ax, radius*sin(theta_canvas), radius*cos(theta_canvas), z_canvas, 'FaceColor', [0.2 0.55 0.9], 'FaceAlpha', 0.12, 'EdgeColor', 'none')

    for segment = 1:6
        curve = rotorGeometry.curve_3d{k, segment};
        plot3(r_ax, curve(:, 1), curve(:, 2), curve(:, 3), 'Color', section_color(k, :), 'LineWidth', 1.5)
    end
end

axis equal
grid on
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
title('Rotor')
view(38, 24)

s_ax = nexttile(wrap_layout);
hold(s_ax, 'on')

for segment = 1:6
    surf(s_ax, statorGeometry.mantle{segment}.X, statorGeometry.mantle{segment}.Y, statorGeometry.mantle{segment}.Z, 'FaceColor', [0.65 0.7 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.75)
end

for k = 1:3
    radius = statorGeometry.radius(k);
    stator_outline = [statorGeometry.curves_2d{k, 1}; statorGeometry.curves_2d{k, 2}; statorGeometry.curves_2d{k, 3}; statorGeometry.curves_2d{k, 4}; statorGeometry.curves_2d{k, 5}; statorGeometry.curves_2d{k, 6}];
    stator_theta = statorGeometry.theta_0(k) - stator_outline(:, 2)/radius;
    theta_canvas = linspace(min(stator_theta) - 0.03, max(stator_theta) + 0.03, 30);
    z_canvas = linspace(min(stator_outline(:, 1)), max(stator_outline(:, 1)), 2);
    [theta_canvas, z_canvas] = meshgrid(theta_canvas, z_canvas);
    surf(s_ax, radius*sin(theta_canvas), radius*cos(theta_canvas), z_canvas, 'FaceColor', [0.2 0.55 0.9], 'FaceAlpha', 0.12, 'EdgeColor', 'none')

    for segment = 1:6
        curve = statorGeometry.curves_3d{k, segment};
        plot3(s_ax, curve(:, 1), curve(:, 2), curve(:, 3), 'Color', section_color(k, :), 'LineWidth', 1.5)
    end
end

axis equal
grid on
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
title('Stator')
view(38, 24)
rotate3d on

%plot final rotor blisk (disk is just for visual) and stator

root_profile = [rotorGeometry.section_curve{1, 1}; rotorGeometry.section_curve{1, 2}; rotorGeometry.section_curve{1, 3}; rotorGeometry.section_curve{1, 4}; rotorGeometry.section_curve{1, 5}; rotorGeometry.section_curve{1, 6}];
z_height = max(root_profile(:, 2)) - min(root_profile(:, 2));
z_centre = (min(root_profile(:, 2)) + max(root_profile(:, 2)))/2;
ring_out_radius = rotorGeometry.section_radius(1);
ring_in_radius = 0.95*ring_out_radius;
hub_out_radius = 0.25*ring_out_radius;
hub_bore_radius = 0.50*hub_out_radius;
hub_z_height = z_height;
ring_half_height = 0.10*z_height;

disk_radius = [hub_bore_radius; hub_bore_radius; hub_out_radius; ring_in_radius; ring_in_radius; ring_out_radius; ring_out_radius; ring_in_radius; ring_in_radius; hub_out_radius; hub_bore_radius];
disk_z = [z_centre - hub_z_height/2; z_centre + hub_z_height/2; z_centre + hub_z_height/2; z_centre + ring_half_height; z_centre + z_height/2; z_centre + z_height/2; z_centre - z_height/2; z_centre - z_height/2; z_centre - ring_half_height; z_centre - hub_z_height/2; z_centre - hub_z_height/2];
disk_theta = linspace(0, 2*pi, 161);
[disk_theta, disk_radius] = meshgrid(disk_theta, disk_radius);
disk_z = repmat(disk_z, 1, size(disk_theta, 2));

figure('Color', 'w')
b_ax = axes('Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
hold(b_ax, 'on')
surf(b_ax, disk_radius.*cos(disk_theta), disk_radius.*sin(disk_theta), disk_z, 'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none')

for blade = 0:rotorValues.N - 1
    angle = 2*pi*blade/rotorValues.N;
    rotation = [cos(angle), -sin(angle), 0 ; sin(angle), cos(angle), 0; 0, 0, 1];
    vertices = rotorGeometry.solid.vertices*rotation.';
    patch(b_ax, 'Faces', rotorGeometry.solid.faces, 'Vertices', vertices, 'FaceColor', [0.72 0.72 0.72], 'EdgeColor', 'none')
end

stator_vertices = statorGeometry.solid.vertices;
axial_gap = 6;
stator_move = min(rotorGeometry.solid.vertices(:, 3)) - axial_gap - max(stator_vertices(:, 3));
stator_vertices(:, 3) = stator_vertices(:, 3) + stator_move;
stator_count = max(1, round(cycle.degree_of_admission*statorValues.N));
first_stator = floor((statorValues.N - stator_count)/2);

for vane = 0:stator_count - 1
    angle = 2*pi*(first_stator + vane)/statorValues.N;
    rotation = [cos(angle), -sin(angle), 0; sin(angle), cos(angle), 0; 0, 0, 1];
    vertices = stator_vertices*rotation.';
    patch(b_ax, 'Faces', statorGeometry.solid.faces, 'Vertices', vertices, 'FaceColor', [0.35 0.45 0.65], 'EdgeColor', 'none')
end

axis equal
grid on
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
title('Rotor Blisk and Stator')
camlight headlight
lighting gouraud
view(38, 24)
rotate3d on