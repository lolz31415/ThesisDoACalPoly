%% Write hdf5 covariance matrix to Images 
%Just use imageDatastore, arrayDatastore has been headache with 3D data to
%get loaded into trainnet & for validation data with labels 
% Try First saving into .mat file extension, if this doesn't work then go
% for subfolders and imageDatastore
% note: if error occurs, ensure move back to ThesisDatasets folder for
% h5reading with the HDF2Struct() 

close all; clc;
% clear all; 

%%% User Settings: %%%
% 1 - 3 = training for 25k, 2.5k, 100k samples 
% 4 - 6 = same sotring for validation (alphabetically) 
chosen_train_HDF5_File = 1; 

% 1 - 6: Real+Imag+Phase, Real+Imag+Abs, Real+Imag, Abs+Phase, Real,
% RealUpperTriangular
folder_name_chosen = 5; 
% 3,4 and 5,6 require different ways to make images
% 3,4 use an avg of channel 1+2 for 3rd image (absolute or normal average?)
% Upper Triangular Matrix is Conjugate Symmetric - worth using? 
% 5,6 use single grayscale (note that 

%%% End Settings %%% 


% Entire Absolute Directory Name: Note the "HDF5_File_Location" For
% Organization Purposes
cd("C:\Users\mathe\Documents\MATLAB\Thesis Codes and Figures\ThesisDatasetsAndNeuralNetwork\HDF5_File_Location");
% Obtain file Names
% Normal with 25,000 are named normally, 'ORIGINAL' = 2,500 samples,
% -90_90deg contains 100k train having 2,000 samples for 50 average array
% perturbation multipliers
h5_files = dir(fullfile(pwd,'*.h5')); 
h5_file_names = strings(1,length(h5_files));
for i = 1:length(h5_files)
    h5_file_names(i) = h5_files(i).name;
end

% read HDF5 File Names 
train_filename = h5_file_names(chosen_train_HDF5_File);
% add 3 to get correct validation file name
validate_filename = h5_file_names(chosen_train_HDF5_File+3);

% Just use other datasets for 25,000 case - not 100,000 all angular locations 
% For 2-channeled approach w/ images try 2 methods: 
% Array Datastore and taking and Average of the 2 channels as the 
% 3rd channel -> Not sure how this would be different / compare convolution
% filters? 
foldner_names = ["Training_And_Validation_DataReal_Imag_Phase", ...
    "Training_And_Validation_DataReal_Imag_Abs", ...
    "Training_And_Validation_DataReal_Imag", ...
    "Training_And_Validation_DataAbs_Phase", ...
    "Training_And_Validation_DataReal", ...
    "Training_And_Validation_DataReal_Upper_Triangular"];
chosen_folder_name = foldner_names(folder_name_chosen);

% retry this NOT USING imagesc, but imwrite
if contains(train_filename,"ORIGINAL")==1
    train_folder = 'TrainingData5x5Image_INITIAL2500';
    validate_folder = 'TrainingData5x5Image_INITIAL2500';
elseif contains(train_filename,"-90_90_deg")==1
    train_folder = 'TrainingData-90_90degImage';
    validate_folder = 'ValidationData-90_90degImage';
else % case of original 25,00 samples
    train_folder = 'TrainingData5x5Image';
    validate_folder = 'ValidationData5x5Image';
end

train_data = HDF2Struct(train_filename,train_folder); 
validate_data = HDF2Struct(validate_filename);

% can use later like such: train_data.(train_fieldnames{1}) -> ensure () 
train_fieldnames = fieldnames(train_data);
validate_fieldnames = fieldnames(validate_data);

% check if folders exist (overwrite files, but not folders) 
% create a folder for each perturbation multiplier image 
cd .. % back out to main directory (not HDF5 Directory) 

if ~exist(chosen_folder_name,'dir') % or 'file'
    mkdir(chosen_folder_name)
end
% move directories 
cd(chosen_folder_name)
% create training and validation directories
if ~exist(train_folder,'dir') % or 'file'
    mkdir(train_folder)
end
if ~exist(validate_folder,'dir') % or 'file'
    mkdir(validate_folder)
end

% assume always in [groups# samples#]
train_num_groups_samples =  is_desired_struct(train_data); 
validate_num_groups_samples =  is_desired_struct(validate_data); 

train_group_number = train_data.(train_fieldnames{train_num_groups_samples})(1);
train_sample_number = train_data.(train_fieldnames{train_num_groups_samples})(2);
train_num_samples = train_group_number * train_sample_number;

validate_group_number = validate_data.(validate_fieldnames{validate_num_groups_samples})(1);
validate_sample_number = validate_data.(validate_fieldnames{validate_num_groups_samples})(2);
validate_num_samples = validate_group_number * validate_sample_number;

% Normalize Data before stacking it into images 
% function norm_output = normalize_data(input_data)
% with real,imag,phase = 3 channels
% channels = 3;    
train_data.dataset_cov_real = normalize_data(train_data.dataset_cov_real); 
train_data.dataset_cov_imag = normalize_data(train_data.dataset_cov_imag); 
train_data.dataset_cov_abs = normalize_data(train_data.dataset_cov_abs); 
train_data.dataset_cov_angle = normalize_data(train_data.dataset_cov_angle); 
    
validate_data.dataset_cov_real = normalize_data(validate_data.dataset_cov_real); 
validate_data.dataset_cov_imag = normalize_data(validate_data.dataset_cov_imag); 
validate_data.dataset_cov_abs = normalize_data(validate_data.dataset_cov_abs); 
validate_data.dataset_cov_angle = normalize_data(validate_data.dataset_cov_angle); 

%% Try just writing all files to TrainingData & ValidationData Folders 

% 3 channels, pre-allocate array
N = size(train_data.dataset_arr_perturb_elem,2);

if ~exist(train_folder,'dir')
    mkdir(train_folder);
end

% Can create upscaled images from this, first just try making 5x5 dataset 
%upscaled_img = upsample_image(current_img, 5);
%imshow(upscaled_img)

cd(train_folder);

% 3-Channeled Img, 2-Channeled Img with 3rd Channel Avg. or 1 Channel Img
switch folder_name_chosen 
    case 1 % real, imag, phase 
        train_channel1 = train_data.dataset_cov_real;
        train_channel2 = train_data.dataset_cov_imag;
        train_channel3 = train_data.dataset_cov_angle;
        validate_channel1 = validate_data.dataset_cov_real;
        validate_channel2 = validate_data.dataset_cov_imag;
        validate_channel3 = validate_data.dataset_cov_angle;
        num_channels = 3; 
        
train_data.dataset_cov_abs = normalize_data(train_data.dataset_cov_abs); 

validate_data.dataset_cov_abs = normalize_data(validate_data.dataset_cov_abs); 

    case 2 % real, imag, abs
        train_channel1 = train_data.dataset_cov_real;
        train_channel2 = train_data.dataset_cov_imag;
        train_channel3 = train_data.dataset_cov_abs;
        validate_channel1 = validate_data.dataset_cov_real;
        validate_channel2 = validate_data.dataset_cov_imag;
        validate_channel3 = validate_data.dataset_cov_abs;
        num_channels = 3; 

    % For case 3,4 just call it still 3 channel image created by avg of
    % other 2 channels that contain "real data" using arrayDatastore for
    % training. Compare the 2 results (loss) 
    % This could also be used to average the real/imag into 1D grayscale
    % img. Would we lose features in this scenario? 
    % Weighted avg as in: https://mmuratarat.github.io/2020-05-13/rgb_to_grayscale_formulas 

    case 3 % real, imag
        train_channel1 = train_data.dataset_cov_real;
        train_channel2 = train_data.dataset_cov_imag;
        train_channel3 = (train_channel1+train_channel2)/2; % avg
        validate_channel1 = validate_data.dataset_cov_real;
        validate_channel2 = validate_data.dataset_cov_imag;
        validate_channel3 = (validate_channel1+validate_channel2)/2; % avg
        num_channels = 3; 

    case 4 % abs, phase
        train_channel1 = train_data.dataset_cov_abs;
        train_channel2 = train_data.dataset_cov_angle;
        train_channel3 = (train_channel1+train_channel2)/2; % avg
        validate_channel1 = validate_data.dataset_cov_abs;
        validate_channel2 = validate_data.dataset_cov_angle;
        validate_channel3 = (validate_channel1+validate_channel2)/2; % avg
        num_channels = 3; 

    case 5  % real 
        train_channel1 = train_data.dataset_cov_real;
        validate_channel1 = validate_data.dataset_cov_real;
        num_channels = 1; 

    % case 6 % ignore for now
  
end

if num_channels == 3
    % we want [5,5,3,#samples] for our image - iterate with for loop
    train_start_data_write_time = tic;
    for i = 1:train_num_samples
        
        % FOR STORING ALL IMAGES IN WORKSPACE WHILE COMPUTING 
        % real/imag/phase
        %img_train(:,:,:,i) = cat(3 , train_data.dataset_cov_real(:,:,i) , ...
        %train_data.dataset_cov_imag(:,:,i) , train_data.dataset_cov_angle(:,:,i));
        %current_img = img_train(:,:,:,i); 
        current_img = cat(3 , train_channel1(:,:,i) , ...
        train_channel2(:,:,i) , train_channel3(:,:,i));
        
        % Current Method: only 5x5 images 
        % create file name, then use imwrite
        element_perturbations = join(string(train_data.dataset_arr_perturb_elem(:,:,i)),' ');
        file_name = strcat('real_imag_phase_',string(train_data.dataset_arr_perturb_mult(i)), ...
            '_[',element_perturbations,'].png');
        imwrite(current_img,file_name,'PNG');
    
    end
    % For Creating Dataset: 43.39min = 2.6e3 sec - not image conversion
    train_elapsed_time = toc(train_start_data_write_time); 
    % with 5x5 image size: 8.0171s 2,500 samples
    disp('Training Data Finished Writing');
    disp("Time to write Training Data: "+string(train_num_samples)+"Samples: "+string(train_elapsed_time)+" seconds");

% GrayScale Case: 
elseif num_channels == 1
    % we want [5,5,1,#samples] for our grayscale image 
    train_start_data_write_time = tic;
    for i = 1:train_num_samples
        current_img = train_channel1(:,:,i);
        % create file name, then use imwrite
        element_perturbations = join(string(train_data.dataset_arr_perturb_elem(:,:,i)),' ');
        file_name = strcat('real_imag_phase_',string(train_data.dataset_arr_perturb_mult(i)), ...
            '_[',element_perturbations,'].png');
        imwrite(current_img,file_name,'PNG'); 
    end
    train_elapsed_time = toc(train_start_data_write_time);
    % with 5x5 image size: 8.0171s 
    disp('Training Data Finished Writing');
    disp("Number of Channels: "+string(num_channels));
    disp("Time to write Training Data: "+string(train_num_samples)+"Samples: "+string(train_elapsed_time)+" seconds");
end

cd ..; % go back up one folder

%% Write the Validation Data Images

% uncomment this to run create validation images

if ~exist(validate_folder,'dir')
    mkdir(validate_folder);
end
cd(validate_folder);

if num_channels == 3
    
    validate_start_data_write_time = tic; 
    for i = 1:validate_num_samples
    
        % current use for 5x5 image creation instead of IMAGESC
        current_img = cat(3 , validate_channel1(:,:,i) , ...
        validate_channel2(:,:,i) , validate_channel3(:,:,i));
        
        % Current Method: only 5x5 images 
        % create file name, then use imwrite
        element_perturbations = join(string(validate_data.dataset_arr_perturb_elem(:,:,i)),' ');
        file_name = strcat('real_imag_phase_',string(validate_data.dataset_arr_perturb_mult(i)), ...
            '_[',element_perturbations,'].png');
        imwrite(current_img,file_name,'PNG');
    
    end
    validate_elapsed_time = toc(validate_start_data_write_time);

% real dataset - grayscale
elseif num_channels == 1
    validate_start_data_write_time = tic; 
    for i = 1:validate_num_samples
    
        % current use for 5x5 image creation instead of IMAGESC
        current_img = validate_channel1(:,:,i);
        
        % Current Method: only 5x5 images 
        % create file name, then use imwrite
        element_perturbations = join(string(validate_data.dataset_arr_perturb_elem(:,:,i)),' ');
        file_name = strcat('real_imag_phase_',string(validate_data.dataset_arr_perturb_mult(i)), ...
            '_[',element_perturbations,'].png');
        imwrite(current_img,file_name,'PNG');
    
    end
    validate_elapsed_time = toc(validate_start_data_write_time);
end

disp("Validation Data Finished Writing");
disp("Number of Channels: "+string(num_channels));
disp("Time to write Validation Data: "+string(validate_num_samples)+"Samples: "+string(validate_elapsed_time)+" seconds");

% Tell the user where their produced images are located: 
disp("Training Dataset Images Located in: "+string(fullfile(pwd,chosen_folder_name,train_folder)));
disp("Validation Dataset Images Located in: "+string(fullfile(pwd,chosen_folder_name,validate_folder)));

% go up to main folder (-/ThesisDatasets)
cd ../..

%% Check for string in a desired struct - for finding number of groups and datasets if the value changes
function desired_struct_location = is_desired_struct(data)
    % finds the location where "groups" and "num_samples" is contained in
    % case this string value changes for the name of the field in the
    % struct containing the samples
    % easier than using regular expressions
    field_names = fieldnames(data); 
    desired_struct_location = find(contains(field_names,"groups")&contains(field_names,"num_samples"));
end
%% Function to Normalize Input Data - use for all potential training samples (train + validation) 
function norm_output = normalize_data(input_data)
    % normalizes input data to a range of [0,1] 
    min_val = min(input_data,[],'all');
    max_val = max(input_data,[],'all');
    norm_output = (input_data - min_val) / ( max_val - min_val );
end

%% Viewing x Images on tiles output by random -> Training Imgs
num_imgs_view = 2; 
train_imgs_loc = fullfile(pwd,chosen_folder_name,train_folder); 
img_names = dir(fullfile(train_imgs_loc,'*.png'));
randomIndex = randi(length(img_names), 1, num_imgs_view); % Get random number.
for i = 1:num_imgs_view
    fullFileName = fullfile(train_imgs_loc, img_names(randomIndex(i)).name);
    curr_fig = figure(i*10); 
    img = imread(fullFileName);
    imshow(img)
    % view 100x100 representation of image instead of 5x5 grid
    truesize(curr_fig,[100,100])
    title("Image # "+string(randomIndex(i)));
    impixelinfo
end 

% Previous Upscaling method created a blurred image compared to using the
% "truesize()" function -> Probably because resizing axis handles as well

% Don't need to upscale: can use: 
% Can view a figure in larger size using: "truesize(fig,[mrows,ncols])" -> 
% mrows/ncols is the size of the display of the current figure 
% (100,100 instead of 5x5) 

