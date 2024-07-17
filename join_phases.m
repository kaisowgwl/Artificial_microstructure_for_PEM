function microstructure_jointed = join_phases(microstructure)
    % Convert the image to 8-bit (with scaling)
    img = im2uint8(microstructure);

    % Binarize the image: Convert all non-zero values to 1
    img = img > 0;

    % Generate a temporary file path to save the intermediate result
    outputFilePath = tempname + ".tif";

    % Step 4: Save the processed binary image as a TIFF file using Tiff class
    % Create a new Tiff object for writing
    t = Tiff(outputFilePath, 'w');

    % Set up the TIFF tags
    tagstruct.ImageLength = size(img, 1);
    tagstruct.ImageWidth = size(img, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 1; % Update for binary image
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';

    % Write each slice to the TIFF file
    for k = 1:size(img, 3)
        setTag(t, tagstruct);
        write(t, img(:, :, k));
        if k < size(img, 3)
            writeDirectory(t); % Create a new directory for the next slice
        end
    end

    % Close the TIFF file
    close(t);

    % Load the processed image back
    info = imfinfo(outputFilePath);
    num_images = numel(info);
    microstructure_jointed = zeros(size(img), 'uint8');
    for k = 1:num_images
        microstructure_jointed(:, :, k) = imread(outputFilePath, k, 'Info', info);
    end

    % Convert back to double format
    microstructure_jointed = double(microstructure_jointed);
end
