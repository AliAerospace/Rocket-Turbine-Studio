%sudhof rotor 3d generation
function geometry = generateRotorGeometry(rotorValues, cycle, span)
section_i = [1, span.mean_radius_i, length(span.r)];
section_radius = 1000*span.r(section_i);
section = 3;
segment = 6;


section_stagger = zeros(1, section);
section_rotor = cell(1, section);
section_curve = cell(section, segment);

for k = 1:section
    i = section_i(k);
    beta_i = span.beta_3_span(i);
    beta_o = span.beta_4_span(i);
    section_stagger(k) = calculateStaggerAngle(beta_i, beta_o, rotorValues.n_stagger);
    
    inputs = struct;
    inputs.radius = section_radius(k);
    inputs.axial_chord = rotorValues.chord*cosd(section_stagger(k));
    inputs.tangential_chord = rotorValues.chord*sind(section_stagger(k));
    inputs.unguided_turning = rotorValues.unguided_turning;
    inputs.beta_in = 90 - beta_i;
    inputs.epsilon_in = cycle.gamma_ri/2;
    inputs.r_le = rotorValues.r_le_ratio*inputs.axial_chord;
    inputs.beta_out = beta_o - 90;
    inputs.r_te = rotorValues.r_te_ratio*inputs.axial_chord;
    inputs.N = rotorValues.N;
    section_pitch = 2*pi*inputs.radius/inputs.N;
    inputs.throat = section_pitch*cosd(inputs.beta_out) - 2*inputs.r_te;

    section_pritchard{k} = createPritchardRotor(inputs);
    [section_sudhof_inputs{k}, section_fit{k}] = fitPritchardToSudhof(section_pritchard{k}, rotorValues.t_ri, rotorValues.gamma_turning_ri);
    section_rotor{k} = createSudhofRotor(section_sudhof_inputs{k});
    
    section_curve{k, 1} = section_rotor{k}.P1_P2_curve;
    section_curve{k, 2} = section_rotor{k}.P2_P3_curve;
    section_curve{k, 3} = section_rotor{k}.P3_P4_curve;
    section_curve{k, 4} = section_rotor{k}.P4_P5_curve;
    section_curve{k, 5} = section_rotor{k}.P5_P6_curve;
    section_curve{k, 6} = section_rotor{k}.P6_P1_curve;
end

%cylindrical projection
theta_0 = zeros(1, section);
curve_3d = cell(section, segment);
for k = 1:section
    for j = 1:segment
        curve_2d = section_curve{k, j};
        theta = theta_0(k) + curve_2d(:, 1)/section_radius(k);
        %cylinderical to cartesian
        X = section_radius(k)*sin(theta);
        Y = section_radius(k)*cos(theta);
        Z = curve_2d(:, 2);
        curve_3d{k, j} = [X, Y, Z];
    end
end

profile_degree = 3;
point = size(curve_3d{1, 1}, 1);
u = linspace(0, 1, point).';
knots = zeros(1, point + profile_degree + 1);
knots(end - profile_degree:end) = 1;

%averaging. choosing internal knots for interpolation
for j = 1:point - profile_degree - 1
    knots(j + profile_degree + 1) = mean(u(j + 1:j + profile_degree));
end

%b spline basis function nurbs eqn 2.5 (top one)
basis = zeros(point, point + profile_degree);
for j = 1:size(basis, 2)
    basis(:, j) = u >= knots(j) & u < knots(j + 1);
end

%coxde boor nurbs eqn 2.5 (bottom one)
for degree = 1:profile_degree
    next_basis = zeros(point, point + profile_degree - degree);
    for j = 1:size(next_basis, 2)
        if knots(j + degree) - knots(j) ~=0
            next_basis(:, j) = next_basis(:, j) + ((u - knots(j))/(knots(j + degree) - knots(j))).*basis(:, j);
        end
        if knots(j + degree + 1) - knots(j + 1) ~= 0
            next_basis(:, j) = next_basis(:, j) + ((knots(j + degree + 1) - u)/(knots(j + degree + 1) - knots(j + 1))).*basis(:, j + 1);
        end
    end
    basis = next_basis;
end
basis(end, :) = 0;
basis(end, end) = 1;

%nurbs eqn 9.2
control = cell(section, segment);
for k = 1:section
    for j = 1:segment
        control{k, j} = basis\curve_3d{k, j};
    end
end

%bernstein polynomials
v = linspace(0, 1, 41).';
M0 = (1 - v).^2;
M1 = 2*v.*(1 - v);
M2 = v.^2;
mantle = cell(1, segment);

for j = 1:segment
    A0 = control{1, j};
    A2 = control{3, j};
    A1 = 2*control{2, j} - 0.5*(A0 + A2);

    mantle{j}.X = M0*(basis*A0(:, 1)).' + M1*(basis*A1(:, 1)).' + M2*(basis*A2(:, 1)).';
    mantle{j}.Y = M0*(basis*A0(:, 2)).' + M1*(basis*A1(:, 2)).' + M2*(basis*A2(:, 2)).';
    mantle{j}.Z = M0*(basis*A0(:, 3)).' + M1*(basis*A1(:, 3)).' + M2*(basis*A2(:, 3)).';
end


%mesh construction and closing the root and tip regions
ends = cell(1, 2);
end_section = [1, section];
vert = [];
face = [];

for j = 1:segment
    %creates triangle
    tri_mesh = surf2patch(mantle{j}.X, mantle{j}.Y, mantle{j}.Z, 'triangles');
    offset = size(vert, 1);
    vert = [vert; tri_mesh.vertices];
    face = [face; tri_mesh.faces + offset];
end

for end_i = 1:2
    k = end_section(end_i);
    %converting the 6 curves into 1 closed curve
    outline = [section_curve{k, 1}; 
            section_curve{k, 2}(2:end, :);
               section_curve{k, 3}(2:end, :);
               section_curve{k, 4}(2:end, :);
               section_curve{k, 5}(2:end, :);
               section_curve{k, 6}(2:end, :)];

    %ensure first and last points close
    if norm(outline(end, :) - outline(1, :)) < 1e-12
        outline(end, :) = [];
    end

    shape = polyshape(outline(:, 1), outline(:, 2), 'Simplify', false);
    tri = triangulation(shape);

    theta = theta_0(k) + tri.Points(:, 1)/section_radius(k);

    end_vert = [section_radius(k)*sin(theta), section_radius(k)*cos(theta), tri.Points(:, 2)];
    end_face = tri.ConnectivityList;
    first = end_face(1, :);
    normal = cross(end_vert(first(2), :) - end_vert(first(1), :), end_vert(first(3), :) - end_vert(first(1), :));

    centre = mean(end_vert(first, :), 1);
    radial = [centre(1), centre(2), 0];
    sense = 2*(end_i == 2) - 1;

    if sense*dot(normal, radial) < 0
        end_face(:, [2, 3]) = end_face(:, [3, 2]);
    end

    ends{end_i}.vert = end_vert;
    ends{end_i}.face = end_face;

    offset = size(vert, 1);
    vert = [vert; end_vert];
    face = [face; end_face + offset];
end

%combining vertices that are extremely close to within 10^-10m
[vert, ~, vertex_map] = uniquetol(vert, 1e-10, 'ByRows', true, 'DataScale', 1);
face = vertex_map(face);
%removing invalid triangles
nondegen_faces = face(:, 1) ~= face(:, 2) & face(:, 2) ~= face(:, 3) & face(:, 3) ~= face(:, 1);

%neither of the vertices can be the same for a respective triangle
face = face(nondegen_faces, :);

blade_solid.vertices = vert;
blade_solid.faces = face;

%ouputs
geometry.section_i = section_i;
geometry.section_radius = section_radius;
geometry.section_stagger = section_stagger;
geometry.section_rotor = section_rotor;
geometry.section_pritchard = section_pritchard;
geometry.section_sudhof_inputs = section_sudhof_inputs;
geometry.section_fit = section_fit;
geometry.section_curve = section_curve;
geometry.curve_3d = curve_3d;
geometry.control = control;
geometry.mantle = mantle;
geometry.ends = ends;
geometry.solid = blade_solid;
end
