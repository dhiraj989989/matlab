% Copyright 2010 George F. Burkhard, Eric T. Hoke, Stanford University
% Licensed under GNU GPL v3 – this means you can use, modify, and share this code freely,
% provided you give proper credit and share derivatives under the same license.

% This program calculates the short-circuit current (Jsc) for a solar cell device
% where the thickness of one layer is varied, assuming 100% internal quantum efficiency (IQE),
% using the Transfer Matrix Method (TMM).

% Based on methods from:
% J. Appl. Phys Vol 86, No. 1 (1999) p.487 and JAP 93 No. 7 p. 3693

function TransferMatrix_VaryThickness

%-------------------- USER DEFINED PARAMETERS ----------------------

lambda = 350:1:1000; % Wavelength range (nm) for simulation

stepsize = 1; % Spatial resolution (nm) for E-field and absorption calculations

% Layer names in the device from incident light side to back contact
layers = {'Glass', 'ITO', 'PEDOT', 'FASnI3', 'ICBA', 'Ag'}; 

% Thickness of each corresponding layer in nm
thicknesses = [0, 110, 40, 120, 40, 100];

activeLayer = 4; % Layer index where absorption leads to photocurrent
LayerVaried = 5; % Layer index whose thickness is to be varied

% Range of thickness values for the layer being varied (in nm)
VaryThickness = 10:30:200;

%---------------------------------------------------------------

% Load AM1.5 solar spectrum (in mW/cm²·nm)
AM15_data = xlsread('AM15.xls'); 
AM15 = interp1(AM15_data(:,1), AM15_data(:,2), lambda, 'linear', 'extrap');

% Load complex refractive index (n + ik) for each material in each wavelength
n = zeros(length(layers), length(lambda));
for index = 1:length(layers)
    n(index,:) = LoadRefrIndex(layers{index}, lambda);
end

% Physical constants
h = 6.62606957e-34;  % Planck's constant (J·s)
c = 2.99792458e8;    % Speed of light (m/s)
q = 1.60217657e-19;  % Electron charge (C)

% Fresnel coefficients at air/glass interface
T_glass = abs(4 * 1 .* n(1,:) ./ (1 + n(1,:)).^2); 
R_glass = abs((1 - n(1,:)) ./ (1 + n(1,:))).^2;

% Initialize layer thickness and output current vector
t = thicknesses;
t(1) = 0; % First layer's thickness (glass) set to 0 for modeling purposes
Jsc = 0 * VaryThickness; % Initialize Jsc output array

%------------------ SWEEP OVER VARIED THICKNESS ------------------
for ThickInd = 1:length(VaryThickness)
    t(LayerVaried) = VaryThickness(ThickInd) % Set current thickness value

    % Cumulative thickness to locate layer boundaries
    t_cumsum = cumsum(t);
    
    % Position grid within the device
    x_pos = (stepsize/2):stepsize:sum(t); 

    % Determine which layer each position belongs to
    x_mat = sum(repmat(x_pos, length(t), 1) > repmat(t_cumsum', 1, length(x_pos)), 1) + 1; 

    % Initialize reflection (R) and field matrix (E)
    R = lambda * 0;
    E = zeros(length(x_pos), length(lambda));

    %-------- LOOP OVER WAVELENGTHS -------------
    for l = 1:length(lambda)
        % Build the global transfer matrix S across all layers
        S = I_mat(n(1,l), n(2,l));
        for matindex = 2:(length(t)-1)
            S = S * L_mat(n(matindex,l), t(matindex), lambda(l)) * ...
                    I_mat(n(matindex,l), n(matindex+1,l));
        end

        % Calculate reflectance and transmission through glass
        R(l) = abs(S(2,1)/S(1,1))^2;
        T(l) = abs(2 / (1 + n(1,l))) / sqrt(1 - R_glass(l) * R(l));

        % Compute normalized electric field inside each layer
        for material = 2:length(t)
            xi = 2 * pi * n(material,l) / lambda(l);
            dj = t(material);
            x_indices = find(x_mat == material); 
            x = x_pos(x_indices) - t_cumsum(material - 1); 

            S_prime = I_mat(n(1,l), n(2,l));
            for matindex = 3:material
                S_prime = S_prime * ...
                    L_mat(n(matindex-1,l), t(matindex-1), lambda(l)) * ...
                    I_mat(n(matindex-1,l), n(matindex,l));
            end

            S_doubleprime = eye(2);
            for matindex = material:(length(t)-1)
                S_doubleprime = S_doubleprime * ...
                    I_mat(n(matindex,l), n(matindex+1,l)) * ...
                    L_mat(n(matindex+1,l), t(matindex+1), lambda(l));
            end

            % Compute the electric field at each position (normalized)
            E(x_indices,l) = T(l) * ...
                (S_doubleprime(1,1)*exp(-1i*xi*(dj-x)) + ...
                 S_doubleprime(2,1)*exp(1i*xi*(dj-x))) ./ ...
                (S_prime(1,1)*S_doubleprime(1,1)*exp(-1i*xi*dj) + ...
                 S_prime(1,2)*S_doubleprime(2,1)*exp(1i*xi*dj));
        end
    end

    % Calculate absorption coefficient (alpha, in cm^-1) for each layer
    a = zeros(length(t), length(lambda));
    for matindex = 2:length(t)
        a(matindex,:) = 4 * pi * imag(n(matindex,:)) ./ (lambda * 1e-7);
    end

    % Power absorbed in the active layer (Q in mW/cm^3/nm)
    ActivePos = find(x_mat == activeLayer);
    Q = repmat(a(activeLayer,:) .* real(n(activeLayer,:)) .* AM15, length(ActivePos), 1) ...
        .* abs(E(ActivePos,:)).^2;

    % Exciton generation rate Gxl (per s·cm^3·nm)
    Gxl = (Q * 1e-3) .* repmat(lambda * 1e-9, length(ActivePos), 1) / (h * c);

    % Integrate Gxl over wavelength to get Gx (per s·cm^3)
    if length(lambda) == 1
        lambdastep = 1;
    else
        lambdastep = (max(lambda) - min(lambda)) / (length(lambda) - 1);
    end
    Gx = sum(Gxl, 2) * lambdastep;

    % Compute short-circuit current Jsc for this thickness
    Jsc(ThickInd) = sum(Gx) * stepsize * 1e-7 * q * 1e3; % mA/cm^2
end

%------------------ PLOT RESULTS -------------------
figure(1)
plot(VaryThickness, Jsc, 'LineWidth', 2)
title('Current Density vs. Layer Thickness')
xlabel('Thickness of Varied Layer (nm)')
ylabel('J_{SC} (mA/cm^2)')

% Export data to base workspace for further analysis
assignin('base', 'Thickness', VaryThickness);
assignin('base', 'Jsc', Jsc);

%------------------ HELPER FUNCTIONS -------------------

function I = I_mat(n1,n2)
% Returns interface matrix for two layers with complex indices n1 and n2
r = (n1 - n2) / (n1 + n2); 
t = 2 * n1 / (n1 + n2); 
I = [1 r; r 1] / t;

function L = L_mat(n, d, lambda)
% Returns propagation matrix for layer with index n, thickness d, and wavelength lambda
xi = 2 * pi * n / lambda;
L = [exp(-1i * xi * d), 0; 0, exp(1i * xi * d)];

function ntotal = LoadRefrIndex(name, wavelengths)
% Reads complex index (n + ik) of a material from the Excel library
[IndRefr, IndRefr_names] = xlsread('Index_of_Refraction_library.xls');
file_wavelengths = IndRefr(:, strmatch('Wavelength', IndRefr_names));
n = IndRefr(:, strmatch(strcat(name,'_n'), IndRefr_names));
k = IndRefr(:, strmatch(strcat(name,'_k'), IndRefr_names));
n_interp = interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp = interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');
ntotal = n_interp + 1i * k_interp;

