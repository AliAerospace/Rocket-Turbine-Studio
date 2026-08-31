function exportBladeRowsIGES(rotorGeometry, statorGeometry, rotorValues, statorValues, cycle, rotor_name, stator_name, folder)

if ~exist(folder, 'dir')
    mkdir(folder)
end

%rotor blade row
rotor_surfs = createBladeRow(rotorGeometry.control, rotorValues.N, rotorValues.N);

%stator blade row
stator_N = round(cycle.degree_of_admission*statorValues.N);
stator_surfs = createBladeRow(statorGeometry.control, stator_N, statorValues.N);

current_folder = pwd;
cd(folder)

igesout(rotor_surfs, rotor_name);
igesout(stator_surfs, stator_name);
cd(current_folder)

fprintf('Created %s\n', fullfile(folder, rotor_name))
fprintf('Created %s\n', fullfile(folder, stator_name))

end


function surfaces = createBladeRow(control, blade_count, total_positions)
segment_count = size(control, 2);
point_count = size(control{1, 1}, 1);

%profile knot vector
degree = 3;
parameter = linspace(0, 1, point_count).';
knots = zeros(1, point_count + degree + 1);
knots(end - degree:end) = 1;

for j = 1:point_count - degree - 1
    knots(j + degree + 1) = mean(parameter(j + 1:j + degree));
end

%quadratic span knot vec
span_knots = [0 0 0 1 1 1];
surfaces = cell(1, blade_count*segment_count);
surface_number = 0;

for blade = 0:blade_count - 1
    theta = 2*pi*blade/total_positions;
    rotation = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    for j = 1:segment_count

        %control points root to tip
        A0 = control{1, j};
        Am = control{2, j};
        A2 = control{3, j};

        %quadratic span control point
        A1 = 2*Am - 0.5*(A0 + A2);

        %rotate blade around axis
        A0 = A0*rotation.';
        A1 = A1*rotation.';
        A2 = A2*rotation.';

        %nurbs control points
        coefficients = zeros(4, point_count, 3);
        coefficients(1:3, :, 1) = A0.';
        coefficients(1:3, :, 2) = A1.';
        coefficients(1:3, :, 3) = A2.';
        coefficients(4, :, :) = 1;

        surface_number = surface_number + 1;
        surfaces{surface_number} = nrbmak(coefficients, {knots, span_knots});
    end
end

end