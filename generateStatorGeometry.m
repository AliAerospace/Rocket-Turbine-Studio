%stator 3d generation
function geometry = generateStatorGeometry(statorValues, cycle, span)

stator_section_i = [1, span.mean_radius_i, length(span.r)];
stator_section_radius = 1000*span.r(stator_section_i);
stator_section = 3;
stator_segment = 6;

stator_section_stagger = zeros(1, stator_section);
stator_section_blade = cell(1, stator_section);
stator_section_curve = cell(stator_section, stator_segment);
stator_line_parameter = linspace(0, 1, 101).';

for k = 1:stator_section
    i = stator_section_i(k);

    section_inputs = statorValues;
    section_inputs.radius = stator_section_radius(k);
    section_inputs.beta_1 = cycle.beta_1;
    section_inputs.beta_2 = span.beta_2_span(i);
    stator_section_blade{k} = createSudhofStator(section_inputs);
    stator_section_stagger(k) = stator_section_blade{k}.stagger_angle;

    stator_section_curve{k, 1} = stator_section_blade{k}.P1_P2_curve;
    stator_section_curve{k, 2} = stator_section_blade{k}.P2_P3_curve;
    stator_section_curve{k, 3} = stator_section_blade{k}.P3_P4_curve;
    stator_section_curve{k, 4} = stator_section_blade{k}.P4_P5_curve;
    stator_section_curve{k, 5} = stator_section_blade{k}.P5_tip_curve;
    tip_line = stator_section_blade{k}.tip_P1_curve;
    stator_section_curve{k, 6} = (1 - stator_line_parameter)*tip_line(1, :) + stator_line_parameter*tip_line(2, :);

    for j = 1:stator_segment
        stator_section_curve{k, j} = stator_section_curve{k, j} - stator_section_blade{k}.te_tip;
    end
end

%cylindrical projection
stator_theta_0 = zeros(1, stator_section);
stator_curve_3d = cell(stator_section, stator_segment);
for k = 1:stator_section
    for j = 1:stator_segment
        curve_2d = stator_section_curve{k, j};
        theta = stator_theta_0(k) - curve_2d(:, 2)/stator_section_radius(k);
        %cylinderical to cartesian
        X = stator_section_radius(k)*sin(theta);
        Y = stator_section_radius(k)*cos(theta);
        Z = curve_2d(:, 1);

        stator_curve_3d{k, j} = [X, Y, Z];
    end
end

stator_profile_degree = 3;
stator_point = size(stator_curve_3d{1, 1}, 1);
stator_u = linspace(0, 1, stator_point).';
stator_knots = zeros(1, stator_point + stator_profile_degree + 1);
stator_knots(end - stator_profile_degree:end) = 1;

%averaging. choosing internal knots for interpolation
for j = 1:stator_point-stator_profile_degree-1
    stator_knots(j + stator_profile_degree + 1) = mean(stator_u(j + 1:j + stator_profile_degree));
end

%b spline basis function nurbs eqn 2.5 (top one)
stator_basis = zeros(stator_point, stator_point + stator_profile_degree);
for j = 1:size(stator_basis, 2)
    stator_basis(:, j) = stator_u >= stator_knots(j) & stator_u < stator_knots(j + 1);
end

%coxde boor nurbs eqn 2.5 (bottom one)
for degree = 1:stator_profile_degree
    next_basis = zeros(stator_point, stator_point + stator_profile_degree - degree);
    for j = 1:size(next_basis, 2)
        if stator_knots(j + degree) - stator_knots(j) ~=0
            next_basis(:, j) = next_basis(:, j) + ((stator_u - stator_knots(j))/(stator_knots(j + degree) - stator_knots(j))).*stator_basis(:, j);
        end
        if stator_knots(j + degree + 1) - stator_knots(j + 1) ~= 0
            next_basis(:, j) = next_basis(:, j) + ((stator_knots(j + degree + 1) - stator_u)/(stator_knots(j + degree + 1) - stator_knots(j + 1))).*stator_basis(:, j + 1);
        end
    end
    stator_basis = next_basis;
end
stator_basis(end, :) = 0;
stator_basis(end, end) = 1;

%nurbs eqn 9.2
stator_control = cell(stator_section, stator_segment);
for k = 1:stator_section
    for j = 1:stator_segment
        stator_control{k, j} = stator_basis\stator_curve_3d{k, j};
    end
end

%bernstein polynomials
stator_v = linspace(0, 1, 41).';
stator_M0 = (1 - stator_v).^2;
stator_M1 = 2*stator_v.*(1 - stator_v);
stator_M2 = stator_v.^2;
stator_mantle = cell(1, stator_segment);

for j = 1:stator_segment
    A0 = stator_control{1, j};
    A2 = stator_control{3, j};
    A1 = 2*stator_control{2, j} - 0.5*(A0 + A2);

    stator_mantle{j}.X = stator_M0*(stator_basis*A0(:, 1)).' + stator_M1*(stator_basis*A1(:, 1)).' + stator_M2*(stator_basis*A2(:, 1)).';
    stator_mantle{j}.Y = stator_M0*(stator_basis*A0(:, 2)).' + stator_M1*(stator_basis*A1(:, 2)).' + stator_M2*(stator_basis*A2(:, 2)).';
    stator_mantle{j}.Z = stator_M0*(stator_basis*A0(:, 3)).' + stator_M1*(stator_basis*A1(:, 3)).' + stator_M2*(stator_basis*A2(:, 3)).';
end

%mesh construction and closing the root and tip regions
stator_ends = cell(1, 2);
stator_end_section = [1, stator_section];
stator_vert = [];
stator_face = [];

for j = 1:stator_segment
    %creates triangle
    tri_mesh = surf2patch(stator_mantle{j}.X, stator_mantle{j}.Y, stator_mantle{j}.Z, 'triangles');
    offset = size(stator_vert, 1);
    stator_vert = [stator_vert; tri_mesh.vertices];
    stator_face = [stator_face; tri_mesh.faces + offset];
end

for end_i = 1:2
    k = stator_end_section(end_i);
    %converting all curves into 1 closed curve
    outline = [stator_section_curve{k, 1};
            stator_section_curve{k, 2}(2:end, :);
               stator_section_curve{k, 3}(2:end, :);
               stator_section_curve{k, 4}(2:end, :);
               stator_section_curve{k, 5}(2:end, :);
               stator_section_curve{k, 6}(2:end, :)];

    %ensure first and last points close
    if norm(outline(end, :) - outline(1, :)) < 1e-12
        outline(end, :) = [];
    end

    shape = polyshape(outline(:, 1), outline(:, 2), 'Simplify', false);
    tri = triangulation(shape);

    theta = stator_theta_0(k) - tri.Points(:, 2)/stator_section_radius(k);

    end_vert = [stator_section_radius(k)*sin(theta), stator_section_radius(k)*cos(theta), tri.Points(:, 1)];
    end_face = tri.ConnectivityList;
    first = end_face(1, :);
    normal = cross(end_vert(first(2), :) - end_vert(first(1), :), end_vert(first(3), :) - end_vert(first(1), :));

    centre = mean(end_vert(first, :), 1);
    radial = [centre(1), centre(2), 0];
    sense = 2*(end_i == 2) - 1;

    if sense*dot(normal,radial) < 0
        end_face(:, [2, 3]) = end_face(:, [3, 2]);
    end

    stator_ends{end_i}.vert = end_vert;
    stator_ends{end_i}.face = end_face;

    offset = size(stator_vert, 1);
    stator_vert = [stator_vert; end_vert];
    stator_face = [stator_face; end_face + offset];
end

%combining vertices that are extremely close to within 10^-10m
[stator_vert, ~, stator_vertex_map] = uniquetol(stator_vert, 1e-10, 'ByRows', true, 'DataScale', 1);
stator_face = stator_vertex_map(stator_face);
%removing invalid triangles
stator_nondegen_faces = stator_face(:, 1) ~= stator_face(:, 2) & stator_face(:, 2) ~= stator_face(:, 3) & stator_face(:, 3) ~= stator_face(:, 1);

%neither of the vertices can be the same for a respective triangle
stator_face = stator_face(stator_nondegen_faces,:);

stator_solid.vertices = stator_vert;
stator_solid.faces = stator_face;

%outputs
geometry.section_i = stator_section_i;
geometry.radius = stator_section_radius;
geometry.stagger = stator_section_stagger;
geometry.sections = stator_section_blade;
geometry.curves_2d = stator_section_curve;
geometry.curves_3d = stator_curve_3d;
geometry.theta_0 = stator_theta_0;
geometry.basis = stator_basis;
geometry.control = stator_control;
geometry.mantle = stator_mantle;
geometry.ends = stator_ends;
geometry.solid = stator_solid;
end
