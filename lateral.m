% Script to simulate a Fourier light field psf for an emitter at
% position (x,y,z) in object space. The emitter can be isotropic (using an
% aplanatic collimation apodization instead of flat field at bfp) or a
% dipole with orientation (phi,theta).
% 
% The emitter can be in an index-matched medium, i.e. light doesn't
% encounter any refractive index interfaces as it travels from the emitter
% to the objective lens.
% 
% In case of the dipole, it can also lie on a planar refractive index interface
% (e.g. glass-water, glass-air). It should also be able to lie some
% distance Delta away from the interface, but I have to double-check some
% signs.
% 
% The bfp will shrink when there is a refractive index mismatch.
% 
% The MLA can be hexagonal or square. The number of lenses in the pupil can
% be changed by adjusting the focal length of the fourier lens as would be
% done experimentally. The MLA can be positioned with a microlens centered
% on the optical axis, or a corner between microlenses centered on the
% axis. From there, the array can be rotated around the optical axis and
% translated some amount in x and y.
% 
% Different camera models can be used:
%   - nonoise: no noise, no camera model
%   - ideal: a perfect camera, i.e. only poisson noise
%   - sCMOS: simplified sCMOS noise model that ignores pixel-dependence
%   - EMCCD: model for emccd camera (after the SMLM fightclub methods section)

clear all; close all; clc
addpath('lib') % add the utility functions

% System parameters (all in si units unless otherwise specified in comments)
array_type    = 'square'; % 'square' or 'hexagonal'
array_pitch   = 1*1e-3; % along x-direction (parallel with optical table)
array_origin  = 'center'; % 'center' or 'corner' of a microlens in the origin
array_rot     = 0; % rotation of the MLA in degrees (not radians)
array_dx      = 0; % displacement in x
array_dy      = 0; % displacement in y
array_size    = 10e-3; % physical dimensions of the optic (assumes it's a square)

NA            = 1.42; % numerical aperture of the objective
f_obj         = 3.3e-3; % focal length of the objective lens (f_obj = M/f_tube)
f_tube        = 200e-3; % focal length of the tube lens
f_4f          = 75e-3; % focal length of the fourier lens
f_u           = 36.7e-3; % focal length of the microlenses in the microlens array
cam_pixsize   = 11*1e-6;

emitterType   = 'isotropic'; % 'isotropic' or 'dipole'
wavelength    = 630e-9; % wavelength of emitted light
n_medium      = 1.33; % refractive index of the medium in which the molecule has been embedded
n_coverglass  = 1.518; % refractive index of the cover glass and immersion medium
Delta         = 0; % distance between dipole emitter and refractive index interface
                   
photons       = 10000; % mean of poisson dist for number of photons/bead/frame (or exact if cameraType is 'nonoise')

% camera parameters
cameraType = 'nonoise'; % 'nonoise', 'ideal' (only poisson), 'sCMOS' or 'EMCCD'
cameraNoiseParams.QE = 1;
cameraNoiseParams.offset = 0;
cameraNoiseParams.e_adu = 1;
cameraNoiseParams.sigmaReadNoise = 0;
cameraNoiseParams.bit = 16;
% cameraNoiseParams.c = 250; % for EMCCD, electron-multiplying gain
% cameraNoiseParams.emgain = 0.002; % for EMCCD, spurious charge (clock-induced charge only, dark counts negligible)


outputdir = fullfile(pwd,'results');
if ~exist(outputdir,'dir'); mkdir(outputdir); end

% set filepath equal to 'false' if you don't want to write the results away
% as a tif stack and return it as a variable 'stack' instead
% filepath = fullfile(outputdir,'psf_stack.tif');
filepath = false;

%% Script

tic

% generate coordinates
% zRange = (-3000:200:3000)*1e-9;

% 2D case
numRepeat = 50;
numBeads = 2; % number of emitters
separations = linspace(1, 10, 10) * 1e-6;
results = zeros(numRepeat, length(separations));

for i = 1 : numRepeat
    for j = 1 : length(separations) % separation between emitters (density)
        t = zeros(numBeads, 1);
        x = zeros(numBeads, 1);
        y = zeros(numBeads, 1);
        z = zeros(numBeads, 1);
        phi = zeros(numBeads, 1);
        theta = (pi/6) * ones(numBeads, 1);
        delta = Delta * ones(numBeads, 1);

        x(1) = -separations(j)/2;
        x(2) = separations(j)/2;

        coordinates = [t, x, y, z, phi, theta, delta];

        % create objects
        MLA = MicroLensArray(array_type,f_u,array_pitch,array_origin,array_rot,[array_dx array_dy],array_size);
        LFM = FourierLFM(NA,f_obj,f_tube,f_4f,MLA,cam_pixsize,n_coverglass,n_medium);
        ACQ = Acquisition(emitterType,coordinates,cameraType,cameraNoiseParams,LFM,MLA,wavelength,photons,Delta);

        % simulate images
        if filepath
            ACQ.simulateDataset(filepath);
            stack = File.readTifStack(filepath);
        else
            stack = ACQ.simulateDataset(filepath);
        end

        % figure;
        % imshow(stack, []);
        % saveas(gcf, 'image.png');

        %% Gaussian fitting
        third = round(length(stack)/3);
        sublens = stack(third:2*third, third:2*third);
        % mu = mean(sublens, 'all');
        % sigma = std(sublens,1,'all');
        % sublens = sublens - mu;
        % [x_candidate, y_candidate] = find(sublens > (mu + 8*sigma));

        surrounding = 4; % surrounding pixels
        x_predict = [];
        y_predict = [];

        maximum = max(sublens, [], 'all');
        [x_max, y_max] = find(sublens == maximum);
        % first point
        x_max_1 = x_max(1);
        y_max_1 = y_max(1);

        % second point
        if length(x_max) > 1 % if more than 1 maximum (2)
            x_max_2 = x_max(2);
            y_max_2 = y_max(2);
        else
            temp = sublens(sublens < maximum);
            maximum = max(temp, [], 'all');
            [x_max, y_max] = find(sublens == maximum);
            x_max_2 = x_max(1);
            y_max_2 = y_max(1);
        end

        while true
            if sqrt((x_max_1-x_max_2)^2 + (y_max_1-y_max_2)^2) < 2 % within 1 pixels
                temp = temp(temp < maximum);
                maximum = max(temp, [], 'all');
                [x_max_2, y_max_2] = find(sublens == maximum);
            else
                break;
            end
        end

        subsection_1 = sublens(x_max_1-surrounding:x_max_1+surrounding, y_max_1-surrounding:y_max_1+surrounding);
        xx_1 = 1 : length(subsection_1);
        yy_1 = 1 : length(subsection_1);
        [fitresult, ~, ~, ~, ~, ~] = fmgaussfit(xx_1, yy_1, subsection_1);
        x_predict(1) = fitresult(5) + x_max_1 - surrounding - 1;
        y_predict(1) = fitresult(6) + y_max_1 - surrounding - 1;

        subsection_2 = sublens(x_max_2-surrounding:x_max_2+surrounding, y_max_2-surrounding:y_max_2+surrounding);
        xx_2 = 1 : length(subsection_2);
        yy_2 = 1 : length(subsection_2);
        [fitresult, ~, ~, ~, ~, ~] = fmgaussfit(xx_2, yy_2, subsection_2);
        x_predict(2) = fitresult(5) + x_max_2 - surrounding - 1;
        y_predict(2) = fitresult(6) + y_max_2 - surrounding - 1;

        temp = x_predict;
        x_predict = y_predict;
        y_predict = temp;

        % pixel -> 3D position
        M = f_tube * f_u / (f_obj * f_4f);

        x_predict = x_predict + third - 1; % sublens -> stack
        y_predict = y_predict + third - 1;

        mid = (length(stack) + 1) / 2;
        x_predict = (x_predict - mid) * cam_pixsize / M; % camera -> sample
        y_predict = (y_predict - mid) * cam_pixsize / M;
        
        separation_predict = sqrt((x_predict(1) - x_predict(2))^2 + (y_predict(1) - y_predict(2))^2);
        error = sqrt((separations(j) - separation_predict)^2);

        results(i, j) = error;

    end
end

mu = mean(results, 1);
sigma = std(results, 1);

%% plot
figure(1);
plot(separations/1e-6, mu/1e-9, 'x-');
xlabel('separation (\mum)');
ylabel('accuracy (nm)');
set(gca, 'fontname', 'consolas', 'fontsize', 11);
grid on;

saveas(gcf, 'accuracy_2D_lateral.png');

figure(2);
plot(separations/1e-6, sigma/1e-9, 'x-');
xlabel('separation (\mum)');
ylabel('precision (nm)');
set(gca, 'fontname', 'consolas', 'fontsize', 11);
grid on;

saveas(gcf, 'precision_2D_lateral.png');

toc