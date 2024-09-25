function [model_name] = generate_helmholtz_tree_model(N, K, tree_type, Fs)
% GENERATE_HELMHOLTZ_TREE_MODEL Generate an NxK Helmholtz Tree Simulink
% Model.
%   Detailed explanation goes here

% Specify the name of the model to create
model_name = sprintf('helmholtz_resonator_tree_%s_%d_%d', tree_type, N, K);

% Check if the file already exists and delete it if it does
if exist(model_name, 'file') == 4
    return;
end

% Create the system
new_system(model_name);

% Blocks
add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference', ...
    [gcs, '/Ref']);
add_block('polimi/Discrete Impulse', [gcs, '/In']);
add_block('nesl_utility/Simulink-PS Converter', [gcs, '/InConv']);
add_block('nesl_utility/PS-Simulink Converter', [gcs, '/pOutConv_0']);
add_block('simulink/Sinks/To Workspace', [gcs, '/pOut_0'], ...
    'VariableName', 'pressure_0');
add_block('fl_lib/Electrical/Electrical Sources/Controlled Voltage Source', ...
    [gcs, '/CVS']);
add_block('fl_lib/Electrical/Electrical Sensors/Voltage Sensor', [gcs, '/VS']);
add_block('nesl_utility/Solver Configuration', [gcs, '/Solver'], ...
    'ResidualTolerance', '1e-09');

% Connections
add_line(gcs, 'In/1', 'InConv/1');
add_line(gcs, 'InConv/RConn1', 'CVS/RConn1');
add_line(gcs, 'VS/LConn1', 'CVS/LConn1');
add_line(gcs, 'VS/RConn2', 'Ref/LConn1');
add_line(gcs, 'CVS/RConn2', 'Ref/LConn1');
add_line(gcs, 'VS/RConn1', 'pOutConv_0/LConn1');
add_line(gcs, 'pOutConv_0/1', 'pOut_0/1');
add_line(gcs, 'Solver/RConn1', 'Ref/LConn1');

% Helmholtz resonators tree
% Root node

for n = 0:N
    for k = 0:K-1

        % IDs
        suffix = sprintf('_%d_%d', [n, k]);
        resonator_id = strcat('Resonator', suffix);
        p_conv_id = strcat('pOutConv', suffix);
        p_out_id = strcat('pOut', suffix);
        U_conv_id = strcat('UOutConv', suffix);
        U_out_id = strcat('UOut', suffix);
        
        % Blocks
        add_block('polimi/Helmholtz Resonator', [gcs, strcat('/', resonator_id)]);
        add_block('nesl_utility/PS-Simulink Converter', [gcs, strcat('/', p_conv_id)]);
        add_block('simulink/Sinks/To Workspace', [gcs, strcat('/', p_out_id)], ...
            'VariableName', strcat('pressure', suffix));
        add_block('nesl_utility/PS-Simulink Converter', [gcs, strcat('/', U_conv_id)]);
        add_block('simulink/Sinks/To Workspace', [gcs, strcat('/', U_out_id)], ...
            'VariableName', strcat('flow', suffix));
        
        % Connections
        add_line(gcs, strcat(U_conv_id, '/LConn1'), strcat(resonator_id, '/LConn3'));
        add_line(gcs, strcat(p_conv_id, '/LConn1'), strcat(resonator_id, '/LConn4'));
        add_line(gcs, strcat(U_conv_id, '/1'), strcat(U_out_id, '/1'));
        add_line(gcs, strcat(p_conv_id, '/1'), strcat(p_out_id, '/1'));
        if n == 0
            add_line(gcs, 'CVS/LConn1', strcat(resonator_id, '/RConn1'));
            add_line(gcs, 'Ref/LConn1', strcat(resonator_id, '/RConn2'));
            break;
        else
            if strcmp(tree_type, 'balanced') && n > 1
                parent_k = k;
            else
                parent_k = 0;
            end
            parent_resonator_id = sprintf('Resonator_%d_%d', [n-1, parent_k]); 
            add_line(gcs, strcat(parent_resonator_id, '/LConn1'), strcat(resonator_id, '/RConn1'));
            add_line(gcs, strcat(parent_resonator_id, '/LConn2'), strcat(resonator_id, '/RConn2'));
        end
    end

end

% Model parameters
set_param(gcs,...
    'Solver', 'FixedStepAuto', ...
    'FixedStep', sprintf('1/%d', Fs));

% Auto arrange
Simulink.BlockDiagram.arrangeSystem(model_name);

% Save the model
save_system(model_name);

end

