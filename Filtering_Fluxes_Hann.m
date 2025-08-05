%% Heat Flux Filtering Using Temporal and Spatial Hann Windows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------- %
% ----- Extraction of Heat Flux Anomalies Associated with Eddies ------  %
% ---------- Temporal and Spatial Filtering Using Hann Window ---------- %
% ---------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Rafaela Rizzi (rafaela.rizzi@gmail.com)
% Created: 2023
%
% Description:
%   This script filters daily heat flux data (LHF or SHF) using both 
%   temporal and spatial Hann windows to isolate heat fluxes associated 
%   with submesoscale eddies.
%
% Inputs:
%   - HeatFlux: 3D matrix (x, y, time) - LHF_2019_2021.mat or 
%   SHF_2019_2021.mat (pre-calculated heat fluxes)
%   - Eddy polygon data in Excel format (planilha_ods.xlsx)
%
% Outputs:
%   - 3D matrices of filtered fluxes (temporally, spatially, or both)
%   - Time vector corresponding to valid days with eddy detection
%
% Requirements:
%   - Pre-processed daily LHF or SHF data stored as HeatFlux
%   - Eddy contour data structured as alternating rows of lon/lat points
%   - Coordinates (xs, ys) must match the heat flux grid
%
% Notes:
%   - Filtering uses a 1D temporal Hann window (temporal filtering) and a 
%   2D spatial Hann window (spatial filtering)
%   - Spatial filtering targets submesoscale (e.g. < 23 km)
%   - Temporal filtering isolates shorter-term anomalies (~5-day window)
%   - Ensure `TIME_REC`, `xs`, and `ys` are available in the workspace
%
% Usage:
%   - Set `flux_type = 'LAT'` for Latent Heat Flux
%   - Set `flux_type = 'SENS'` for Sensible Heat Flux
%   - Adjust `window_days` and `submesoscale_radius_km` as needed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
clc

% --------------------------- USER SETTINGS --------------------------- %

%  Options: 'SENS' for Sensible Heat Flux, 'LAT' for Latent Heat Flux   %
WHICH_FLUX = 'LAT';  
% WHICH_FLUX = 'SENS';  


% --------------------------- Load Data ------------------------------- %

if strcmp(WHICH_FLUX, 'LAT')
    load LHF_2019_2021.mat
elseif strcmp(WHICH_FLUX, 'SENS')
    load SHF_2019_2021.mat
end

%% ---------------------- Temporal Filtering --------------------------- %%

window_days = 5;
temporal_window = hanning(window_days);
temporal_window = temporal_window / sum(tempora_window);

filtered_temp = zeros(size(heat_flux));

for lon = 1:size(heat_flux, 1)
    for lat = 1:size(heat_flux, 2)
        time_series = squeeze(heat_flux(lon, lat, :));
        filtered_series = conv(time_series, temporal_window, 'same');
        filtered_temp(lon, lat, :) = time_series - filtered_series; % Residuals
    end
end

clear lon lat time_series filtered_series temporal_window window_days


%% ---------------- Load Eddy Polygon Data from Excel ------------------ %%

[~, text_data] = xlsread('planilha_ods.xlsx', 2);

lon_poly = str2double(text_data(1:2:end, 4:end));
lat_poly = str2double(text_data(2:2:end, 4:end));
eddy_dates = datetime(text_data(1:2:end, 1), 'InputFormat', 'dd.MM.yyyy');

% Build polygons
for idx = 1:numel(eddy_dates)
    eddy_shapes(idx) = polyshape(lon_poly(idx, :), lat_poly(idx, :));
end

% Remove data within eddy polygons
masked_flux = zeros(size(filtered_temp));
valid_days_mask = zeros(size(TIME_REC));

for t = 1:numel(TIME_REC)
    matched_idx = find(eddy_dates == TIME_REC(t));
    
    if isempty(matched_idx)
        continue;
    end

    day_flux = squeeze(filtered_temp(:, :, t));

    for j = matched_idx'
        mask = inpolygon(xs, ys, lon_poly(j, :), lat_poly(j, :));
        day_flux(mask) = 0;
    end

    masked_flux(:, :, t) = day_flux;
    valid_days_mask(t) = 1;
end

valid_time_idx = find(valid_days_mask);
time_anomalies = TIME_REC(valid_time_idx);
filtered_temp = filtered_temp(:, :, valid_time_idx);
masked_flux = masked_flux(:, :, valid_time_idx);

clear t j matched_idx day_flux mask eddy_shapes lon_poly lat_poly text_data

%% ------------------------ Spatial Filtering -------------------------- %%

submesoscale_radius_km = 23;
window_size = round([submesoscale_radius_km, submesoscale_radius_km]);
spatial_window = hann(window_size(1)) * hann(window_size(2))';
spatial_window = spatial_window / sum(spatial_window(:));

filtered_spatial = zeros(size(masked_flux));

for k = 1:size(filtered_spatial, 3)
    flux_snapshot = squeeze(masked_flux(:, :, k));
    flux_snapshot(isnan(flux_snapshot)) = 0;
    flux_snapshot(isnan(filtered_temp(:, :, 1))) = NaN;

    filtered = conv2(flux_snapshot, spatial_window, 'same');
    filtered(isnan(filtered_temp(:, :, 1))) = NaN;
    filtered_spatial(:, :, k) = filtered;
end

clear k flux_snapshot filtered spatial_window window_size


%% ----------------- Submesoscale Heat Flux Extraction ------------------ %%

HeatFlux_submeso = filtered_temp - filtered_spatial;

%% ---------------------------- Save Output ----------------------------- %%

output_filename = ['HeatFlux_submeso_', WHICH_FLUX, '.mat'];
save(output_filename, 'HeatFlux_submeso', '-v7.3');
