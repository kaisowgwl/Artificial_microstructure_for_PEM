function generate_additive_and_save_tif(microstructure_initial, background_id, additive_id, target_volume_fraction, minimum_dist, randomize_level, erosiondilatation_depth, output_filename)
    % generate_additive_and_save_tif - Adds an additive phase to a microstructure and saves as a .tif stack.
    %
    % Syntax: generate_additive_and_save_tif(microstructure_initial, background_id, additive_id, target_volume_fraction, minimum_dist, randomize_level, erosiondilatation_depth, output_filename)
    %
    % Inputs:
    %   microstructure_initial - 3D matrix representing the initial microstructure
    %   background_id - ID of the background phase (void)
    %   additive_id - ID of the additive phase
    %   target_volume_fraction - Fraction of the volume to be modified
    %   minimum_dist - Minimum distance between selected particles
    %   randomize_level - Level of randomness in selecting particles
    %   erosiondilatation_depth - Depth for erosion and dilation processes
    %   output_filename - The filename for the output .tif stack
    %
    % Outputs:
    %   None (saves output to a .tif stack file)

    % Reset the microstructure to its initial state
    microstructure = microstructure_initial;
    
    % Initialization
    binary_phase = zeros(size(microstructure));
    binary_phase(microstructure == background_id) = 1;
    [distance_bnd2bnd] = Function_particle_size_CPSD_Algorithm(binary_phase); % Calculate distance from boundary to boundary
    
    % Remove small features from the potential additive location
    distance_bnd2bnd(distance_bnd2bnd == 0) = 9e9; % Remove 0 (background) from the potential additive location
    distance_bnd2bnd(distance_bnd2bnd <= minimum_dist * sqrt(3)) = 9e9;
    
    % Initialize variables
    ntot = numel(microstructure);
    current_additive_volumefraction = 0; % Initialize
    increment = 1;
    iter = 0;
    d = 0; % Initialize
    progress = []; % Initialize progress tracking
    
    % Main loop to generate additive phase
    while current_additive_volumefraction < target_volume_fraction
        d = d + increment; % Increment distance
        iter = iter + 1;
        
        if d > minimum_dist * sqrt(3)
            potential_additives = distance_bnd2bnd <= d; % Select all voxels that could be assigned to additive
            n = sum(potential_additives(:)); % How many of them?
            
            if n > 0 % At least one
                if randomize_level == 0
                    if n / ntot <= target_volume_fraction - current_additive_volumefraction % All potential voxel are not enough to match target volume fraction
                        microstructure(potential_additives) = additive_id;
                    else % Randomize among possible location is mandatory not to overshoot the target volume fraction
                        [L, ncluster] = bwlabeln(potential_additives, 18); % Identify cluster among these potential locations
                        cluster_possibility = 1:ncluster;
                        
                        while current_additive_volumefraction <= target_volume_fraction
                            cluster_choice = randi(length(cluster_possibility)); % Random selection of cluster
                            microstructure(L == cluster_possibility(cluster_choice)) = additive_id;
                            current_additive_volumefraction = sum(microstructure(:) == additive_id) / ntot; % Update current volume fraction
                            cluster_possibility(cluster_choice) = []; % Avoid re-assignment
                        end
                    end
                else % Use a portion of them, per cluster randomly chosen
                    [L, ncluster] = bwlabeln(potential_additives, 18); % Identify cluster among these potential locations
                    cluster_possibility = 1:ncluster;
                    volume_increment = 0; % Initialize
                    
                    while volume_increment <= (1 - randomize_level) * n
                        cluster_choice = randi(length(cluster_possibility)); % Random selection of cluster
                        loc = L == cluster_possibility(cluster_choice);
                        microstructure(loc) = additive_id;
                        volume_increment = volume_increment + sum(loc(:));
                        cluster_possibility(cluster_choice) = []; % Avoid re-assignment
                    end
                end
            end
        end
        
        distance_bnd2bnd(microstructure == additive_id) = 9e9; % Remove already assigned additive from the potential additive new location
        current_additive_volumefraction = sum(microstructure(:) == additive_id) / ntot; % Update current volume fraction
        
        % Update progress
        fprintf('Progress: %.2f%% - Target Volume Fraction: %.2f%%\n', 100 * current_additive_volumefraction, 100 * target_volume_fraction);
        
        % Track progress
        progress(iter, 1) = current_additive_volumefraction;
        
        % Break if progress is stuck
        if iter > 10 && length(unique(progress(max(1, end-5):end, 1))) == 1
            break;
        end
    end
    
    % Display final volume fractions
    final_volume_fraction = sum(microstructure(:) == additive_id) / ntot;
    fprintf('Final Additive Volume Fraction: %.2f%%\n', 100 * final_volume_fraction);
    fprintf('Target Volume Fraction: %.2f%%\n', 100 * target_volume_fraction);

    % Save the modified microstructure as a .tif stack
    save_as_tif_stack(microstructure, output_filename);
end

function [distance_bnd2bnd] = Function_particle_size_CPSD_Algorithm(binary_phase)
    % Function_particle_size_CPSD_Algorithm - Placeholder for the actual function to calculate distance from boundary to boundary
    %
    % Inputs:
    %   binary_phase - Binary matrix representing the phase
    %
    % Outputs:
    %   distance_bnd2bnd - Distance from boundary to boundary matrix

    % Replace this with the actual implementation
    % For now, just returning a dummy matrix of the same size as binary_phase
    distance_bnd2bnd = bwdist(~binary_phase); % Example using Euclidean distance transform
end

function save_as_tif_stack(microstructure, filename)
    % save_as_tif_stack - Save the 3D microstructure as a .tif stack
    %
    % Inputs:
    %   microstructure - 3D matrix representing the microstructure
    %   filename - The filename for the output .tif stack

    for z = 1:size(microstructure, 3)
        if z == 1
            imwrite(microstructure(:, :, z), filename, 'WriteMode', 'overwrite', 'Compression', 'none');
        else
            imwrite(microstructure(:, :, z), filename, 'WriteMode', 'append', 'Compression', 'none');
        end
    end
end
