% Copyright information and license statement
% This code is open-source and distributed under GNU GPL v3 license, which means it's free to use, modify, and distribute with proper attribution.

% This program calculates optical electric field profiles, exciton generation,
% and the expected short-circuit current (Jsc) in multilayer thin-film devices (e.g., solar cells)
% using the Transfer Matrix Method (TMM).

% It assumes incident light comes from air (n=1) and the first layer is a thick superstrate.
% If the first layer is not thick, 'Air' should be input as the first layer for accurate modeling.

% The method is adapted from prior literature (JAP articles), which are also cited for proper academic attribution.

function TransferMatrix % Main function definition

%------------------------ USER INPUT SECTION --------------------------
lambda1 = 350:1:1000; % Define wavelength range in nm for which simulation is run
stepsize = 1; % Distance between spatial sampling points in nm for field evaluation

plotWavelengths = [400 600 800 900]; % Wavelengths to visualize E-field intensity plots

% Define material stack and corresponding thickness (nm)
layers = {'Glass'  'ITO' 'PEDOT' 'FASnI3' 'ICBA'   'Ag'}; 
thicknesses = [0 10 10 500 10 100]; 

% Flag to enable/disable generation rate and Jsc calculation
plotGeneration = true; 
activeLayer = 4; % Index of the layer where photo-generation occurs (e.g., perovskite)

%-------------------- END USER INPUT SECTION --------------------------

% Load refractive index data for each layer from external Excel library
n = zeros(size(layers,2),size(lambda1,2)); % Initialize index matrix
for index = 1:size(layers,2)
    n(index,:) = LoadRefrIndex(layers{index},lambda1); % Load n + ik for each layer
end
t = thicknesses; % Assign layer thickness array to variable 't'

% Define fundamental constants
h = 6.62606957e-34; % Planck’s constant in Js
c = 2.99792458e8; % Speed of light in m/s
q = 1.60217657e-19; % Elementary charge in C

% Calculate transmission/reflection at the air/glass interface
T_glass = abs(4*1*n(1,:)./(1+n(1,:)).^2); 
R_glass = abs((1-n(4,:))./(1+n(4,:))).^2; 

% Define spatial grid across the thickness
t(1)=0;
t_cumsum=cumsum(t); % Cumulative thicknesses
x_pos=(stepsize/2):stepsize:sum(t); % Spatial coordinates where E-field is computed
x_mat = sum(repmat(x_pos,length(t),1)>repmat(t_cumsum',1,length(x_pos)),1)+1; % Identify which layer each spatial point belongs to

% Initialize field and reflection variables
R=lambda1*0; 
E=zeros(length(x_pos),length(lambda1)); 

% Loop over each wavelength to compute transfer matrix and E-field
for l = 1:length(lambda1)
    S = I_mat(n(1,l),n(2,l)); % First interface
    for matindex=2:(length(t)-1)
        S = S * L_mat(n(matindex,l),t(matindex),lambda1(l)) * I_mat(n(matindex,l),n(matindex+1,l));
    end
    R(l) = abs(S(2,1)/S(1,1))^2; % Power reflection (JAP Eq. 9)
    T(l) = abs(2/(1+n(1,l)))/sqrt(1 - R_glass(l)*R(l)); % Corrected transmission through substrate

    % Compute field in each layer
    for material = 2:length(t)
        xi = 2*pi*n(material,l)/lambda1(l);
        dj = t(material);
        x_indices = find(x_mat == material); % Points in current material layer
        x = x_pos(x_indices) - t_cumsum(material-1); % Local depth from interface
        S_prime = I_mat(n(1,l),n(2,l));
        for matindex=3:material
            S_prime = S_prime * L_mat(n(matindex-1,l),t(matindex-1),lambda1(l)) * I_mat(n(matindex-1,l),n(matindex,l));
        end
        S_doubleprime = eye(2);
        for matindex=material:(length(t)-1)
            S_doubleprime = S_doubleprime * I_mat(n(matindex,l),n(matindex+1,l)) * L_mat(n(matindex+1,l),t(matindex+1),lambda1(l));
        end
        E(x_indices,l) = T(l)*(S_doubleprime(1,1)*exp(-1i*xi*(dj-x)) + S_doubleprime(2,1)*exp(1i*xi*(dj-x))) ./ ...
                         (S_prime(1,1)*S_doubleprime(1,1)*exp(-1i*xi*dj) + S_prime(1,2)*S_doubleprime(2,1)*exp(1i*xi*dj));
    end 
end

% Total reflection including incoherent reflection at air/glass interface
Reflection = R_glass + T_glass.^2 .* R ./ (1 - R_glass .* R);

% Plot normalized electric field intensity |E|^2
close all
figure(1)
plotString = '';
legendString = cell(1,length(plotWavelengths));
for index = 1:length(plotWavelengths)
    plotString = strcat(plotString, ['x_pos,abs(E(:,', num2str(find(lambda1 == plotWavelengths(index))), ').^2),']);
    legendString{index} = [num2str(plotWavelengths(index)), ' nm'];
end
eval(['plot(', plotString, '''LineWidth'',2)'])
axislimit1 = axis;
for matindex=2:length(t)
    line([sum(t(1:matindex)) sum(t(1:matindex))],[0 axislimit1(4)]);
    text((t_cumsum(matindex)+t_cumsum(matindex-1))/2,0,layers{matindex},'HorizontalAlignment','center','VerticalAlignment','bottom')
end
title('E-field intensity in device');
xlabel('Position (nm)');
ylabel('|E|^2');
legend(legendString);

% Absorption coefficient for each material
a = zeros(length(t),length(lambda1));
for matindex = 2:length(t)
    a(matindex,:) = 4*pi*imag(n(matindex,:)) ./ (lambda1*1e-7);
end

% Plot absorption spectrum for each layer and total reflection
figure(2)
Absorption = zeros(length(t),length(lambda1));
plotString = '';
for matindex=2:length(t)
    Pos = find(x_mat == matindex);
    AbsRate = repmat(a(matindex,:) .* real(n(matindex,:)),length(Pos),1) .* (abs(E(Pos,:)).^2);
    Absorption(matindex,:) = sum(AbsRate,1) * stepsize * 1e-7;
    plotString = strcat(plotString, ['lambda1,Absorption(', num2str(matindex), ',:),']);
end
eval(['plot(', plotString, 'lambda1,Reflection,''LineWidth'',2)'])
title('Absorbed and reflected light fractions');
xlabel('Wavelength (nm)');
ylabel('Intensity Fraction');
legend(layers{2:end}, 'Reflectance');

% Calculate generation rate and Jsc if plotGeneration is true
if plotGeneration == true
    AM15_data = xlsread('AM15.xls'); % Load AM1.5G solar spectrum
    AM15 = interp1(AM15_data(:,1), AM15_data(:,2), lambda1, 'linear', 'extrap');

    figure(3)
    ActivePos = find(x_mat == activeLayer);
    Q = repmat(a(activeLayer,:) .* real(n(activeLayer,:)) .* AM15, length(ActivePos), 1) .* (abs(E(ActivePos,:)).^2);

    Gxl = (Q*1e-3).*repmat(lambda1*1e-9,length(ActivePos),1)/(h*c); % Exciton generation rate

    if length(lambda1) == 1
        lambda1step = 1;
    else
        lambda1step = (max(lambda1)-min(lambda1))/(length(lambda1)-1);
    end
    Gx = sum(Gxl,2)*lambda1step; % Position-dependent generation profile
    plot(x_pos(ActivePos), Gx, 'LineWidth', 2)
    axislimit3 = axis;
    axis([axislimit1(1:2) axislimit3(3:4)])
    for matindex=2:length(t)
        line([sum(t(1:matindex)) sum(t(1:matindex))],[0 axislimit3(4)]);
        text((t_cumsum(matindex)+t_cumsum(matindex-1))/2,0,layers{matindex},'HorizontalAlignment','center','VerticalAlignment','bottom')
    end
    title('Generation Rate in Device') 
    xlabel('Position (nm)');
    ylabel('Generation rate /(s·cm³)');

    Jsc = sum(Gx) * stepsize * 1e-7 * q * 1e3 % Output: Short-circuit current in mA/cm²

    % Calculate parasitic absorption outside active layer
    parasitic_abs = (1 - Reflection - Absorption(activeLayer,:))';

    % Export data to MATLAB base workspace
    assignin('base','absorption',Absorption');
    assignin('base','reflection',Reflection');
    assignin('base','parasitic_abs',parasitic_abs);
    assignin('base','lambda1',lambda1');
end

%------------------ Helper Functions -----------------------

function I = I_mat(n1,n2)
% Interface matrix between n1 and n2
r = (n1 - n2) / (n1 + n2); % Fresnel reflection
t = 2*n1 / (n1 + n2);      % Fresnel transmission
I = [1 r; r 1] / t;        % Normalized interface matrix

function L = L_mat(n,d,lambda1)
% Propagation matrix through layer with thickness d and index n
xi = 2*pi*n/lambda1;
L = [exp(-1i*xi*d) 0; 0 exp(1i*xi*d)];

function ntotal = LoadRefrIndex(name,wavelengths)
% Load complex index n+ik from Excel and interpolate for input wavelengths
[IndRefr,IndRefr_names]=xlsread('Index_of_Refraction_library.xls');
file_wavelengths = IndRefr(:,strmatch('Wavelength',IndRefr_names));
n = IndRefr(:,strmatch(strcat(name,'_n'),IndRefr_names));
k = IndRefr(:,strmatch(strcat(name,'_k'),IndRefr_names));
n_interp = interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp = interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');
ntotal = n_interp + 1i*k_interp; % Complex refractive index

