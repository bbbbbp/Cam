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

numRepeat = 50;
numBeads = 2; % number of emitters
separations = linspace(1, 10, 10) * 1e-6;
results = zeros(numRepeat, length(separations));

for i = 1 : numRepeat
    for j = 1 : length(separations) % separation between emitters (density)
        % generate coordinates
        t = zeros(numBeads, 1);
        x = zeros(numBeads, 1);
        y = zeros(numBeads, 1);
        z = zeros(numBeads, 1);
        phi = zeros(numBeads, 1);
        theta = (pi/6) * ones(numBeads, 1);
        delta = Delta * ones(numBeads, 1);

        z(1) = -separations(j)/2;
        z(2) = separations(j)/2;

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
        
        figure;
        imshow(stack, []);
        saveas(gcf, 'image.png');

        %% Gaussian fitting
        loc_2d_1 = [];
        loc_2d_2 = [];

        no_lens = 1;
        for u = -1 : 1
            for v = -1 : 1
                [x_predict, y_predict] = localise(stack, u, v);
                loc_2d_1(no_lens, :) = [x_predict(1), y_predict(1), u, v];
                loc_2d_2(no_lens, :) = [x_predict(2), y_predict(2), u, v];
                no_lens = no_lens + 1;
            end
        end

        mid = ceil(length(stack) / 2);
        M = f_tube * f_u / (f_obj * f_4f);
        bfp = 2 * NA * f_obj * f_4f / f_tube; % back focal plane diameter in nm
        u_scaling = array_pitch / bfp * 2;

        [loc_3d_1, stdx_1, mse_1] = lfm_fit(loc_2d_1, mid, M, cam_pixsize, array_pitch, u_scaling, NA, n_coverglass);
        [loc_3d_2, stdx_2, mse_2] = lfm_fit(loc_2d_2, mid, M, cam_pixsize, array_pitch, u_scaling, NA, n_coverglass);
        loc_3d = [loc_3d_1'; loc_3d_2'];
        % stdx = [stdx_1'; stdx_2'];
        % mse = [mse_1; mse_2];

        separation_predict = sqrt((loc_3d(1, 1) - loc_3d(2, 1))^2 + (loc_3d(1, 2) - loc_3d(2, 2))^2 + (loc_3d(1, 3) - loc_3d(2, 3))^2);
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

saveas(gcf, 'accuracy_2D_axial.png');

figure(2);
plot(separations/1e-6, sigma/1e-9, 'x-');
xlabel('separation (\mum)');
ylabel('precision (nm)');
set(gca, 'fontname', 'consolas', 'fontsize', 11);
grid on;

saveas(gcf, 'precision_2D_axial.png');

toc