function [loc_3d, stdx, mse] = lfm_fit(loc_2d, mid, M, cam_pixsize, array_pitch, u_scaling, NA, n)
    x = loc_2d(:, 1);
    y = loc_2d(:, 2);
    u = loc_2d(:, 3);
    v = loc_2d(:, 4);
    
    % pixel -> 3D position
    x = (x - mid) * cam_pixsize / M; % camera -> sample
    y = (y - mid) * cam_pixsize / M;

    x = x - array_pitch / M * u;
    y = y - array_pitch / M * v;
    u = u * u_scaling;
    v = v * u_scaling;
    
    b = [x; y];
    
    zero = zeros(size(loc_2d, 1), 1);
    one = ones(size(loc_2d, 1), 1);
    
    alpha_u = u;
    alpha_v = v;
    for i = 1 : numel(alpha_u)
        u_min = alpha_u(i) - u_scaling / 2;
        u_max = alpha_u(i) + u_scaling / 2;
        v_min = alpha_v(i) - u_scaling / 2;
        v_max = alpha_v(i) + u_scaling / 2;
        [um, vm] = meshgrid(u_min : u_scaling/21 : u_max, v_min : u_scaling/21 : v_max);
        rho = sqrt(um.^2 + vm.^2);
        
        phi_u = -(NA/n) .* um ./ sqrt(1 - (NA/n)^2.*rho.^2);
        phi_u(imag(phi_u) ~= 0) = NaN;
        alpha_u(i) = nanmean(phi_u(:));
        phi_v = -(NA/n) .* vm ./ sqrt(1 - (NA/n)^2.*rho.^2);
        phi_v(imag(phi_v) ~= 0) = NaN;
        alpha_v(i) = nanmean(phi_v(:));
    end
    
    A = [[one, zero, alpha_u]; [zero, one, alpha_v]];
    
    [loc_3d, stdx, mse] = lscov(A, b); % Ax = b
end

