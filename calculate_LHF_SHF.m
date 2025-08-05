%% LHF and SHF calculation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------- %
% ----------- Calculation of LHF and SHF using daily data of ---------- %
% -------- air temperature, dew point temperature, air density, ------- %
% ------- air pressure, wind speed, and sea surface temperature ------- %
% --------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Rafaela Rizzi (rafaela.rizzi@gmail.com)
% Created: 2022
%
% Description:
%   This script calculates either the latent (LHF) or sensible (SHF) 
%   heat flux from daily atmospheric and oceanographic data. The user 
%   can define which flux to compute by setting the WHICH_FLUX parameter.
%   
%
% Inputs:
%   - ERA5 Product (Air temperature at 2 m; Dew point temperature at 2 m; 
%   Superficial air pressure);
%   - MUR Product (Sea surface temperature)
%   - CMEMS Product (Wind speed at 10 m and Air density)
%
% Outputs:
%   - HeatFlux: 3D matrix (x, y, time) containing either LHF or SHF
%
% Requirements:
%   - Input data must be previously downloaded and organized by date
%   - Variables must follow standard units:
%       * Temperature in K (converted to °C where needed)
%       * Wind speed in m/s
%       * Surface air pressure in Pa
%       * Air density in kg/m³
%    - The resulting heat flux is given in W/m²
%
% Note:
%   - The geographic domain and time period are defined by the user
%   - Ensure consistent grid resolution and coverage across variables
%   - Coordinate system assumed to be regular lat/lon grid
%
% Usage:
%   Set WHICH_FLUX = 'SENS' for sensible heat flux
%   Set WHICH_FLUX = 'LAT' for latent heat flux

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all

% ----------------------------- EDIT HERE------------------------------ %

%  Options: 'SENS' for Sensible Heat Flux, 'LAT' for Latent Heat Flux   %
WHICH_FLUX = 'LAT';  
WHICH_FLUX = 'SENS';  


% Define the analysis period 
Time_Interval = (datetime(2019,01,01):datetime(2021,12,31))';
Time_Interval.Format = 'dd.MM.yyyy';


% Define the geographic region of interest
region_box = [-50 -38; -29 -22]; % [lon_min lon_max; lat_min lat_max]


% Define start and end indices for looping through Time_Interval 
BEGIN = 1; % starting index; 
F = length(Time_Interval); % final index


% Preallocate spatial dimensions (x = longitude, y = latitude)
% Modify these according to the spatial resolution of the dataset used
xl = 701; % Number of grid points in x-direction (longitude)
yl = 1201; % Number of grid points in y-direction (latitude)

% --------------------------------------------------------------------- %

HeatFlux = ones(xl,yl,F); 

for f = BEGIN:F % Latent and Sensible Heat Flux Calculation for F day
    
% --------------------------------------------------------------------- %
% --------------------------- Input Files: ---------------------------- %
% -------------------------- SST Data (MUR) --------------------------- %
% --------------------------------------------------------------------- %

DataFolder = "..";  % Path to the folder containing SST data
FolderContents = dir(DataFolder);
FileNames = cell(0); % Initialize cell array to store filenames

% Store filenames
for i = 1:numel(FolderContents)
    fileName = FolderContents(i).name;
    
% Ignores '.' and '..' folders

    if strcmp(fileName, '.') || strcmp(fileName, '..')
        continue;
    end
    
    FileNames{end+1} = fileName; % end+1 is used to access the next empty 
    % index in the cell array, ensuring that each filename is stored in a 
    % unique position
end
clear fileName
FileNames = FileNames';

% Match file by date from filename structure
AnalysedDate = Time_Interval(f); 
StartFileName = 'sst';
YearFile = sprintf('%04d', year(AnalysedDate));
MonthFile =  sprintf('%02d', month(AnalysedDate));
DayFile = sprintf('%02d', day(AnalysedDate));
EndFileName = '090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc';

AnalysedFile_Ix = [YearFile,MonthFile,DayFile,EndFileName];

% Find the index of the file that matches the date pattern
Ix = find(strcmp(FileNames,AnalysedFile_Ix));
AnalisedFile = cell2mat(FileNames(Ix));

if isempty(AnalisedFile)
    HeatFlux(:,:,f) = NaN(xl,yl);
    [YearFile, MonthFile, DayFile]
    continue;
end

cd(DataFolder)
lon_s = ncread(AnalisedFile, 'lon');
lat_s = ncread(AnalisedFile, 'lat');
sst = ncread(AnalisedFile, 'analysed_sst'); 
time = ncread(AnalisedFile,'time');

% Select data subset limited to specified region coordinates
lon_mask = (lon_s >= region_box(1,1)) & (lon_s <= region_box(1,2));
lat_mask = (lat_s >= region_box(2,1)) & (lat_s <= region_box(2,2));

sst = sst(lon_mask, lat_mask);
lon_s = lon_s(lon_mask, lat_mask);
lat_s = lat_s(lon_mask, lat_mask);

[xs, ys] = meshgrid(lon_s,lat_s);

% --------------------------------------------------------------------- %
% --------------------------- Input Files: ---------------------------- %
% ----------------- Wind and Air Density Data (CMEMS) ----------------- %
% --------------------------------------------------------------------- %

DataFolder = "..";  % Path to the folder containing Wind data
FolderContents = dir(DataFolder);
FileNames = cell(0); % Initialize cell array to store filenames

% Store filenames
for i = 1:numel(FolderContents)
    fileName = FolderContents(i).name;
        if strcmp(fileName, '.') || strcmp(fileName, '..')
        continue;
        end
    FileNames{end+1} = fileName;
end
FileNames = FileNames';

% Match file by date from filename structure
AnalysedDate = Time_Interval(f); 
previousDate = AnalysedDate-1;
StartFileName = 'cmems_obs-wind_glo_phy_my_l4_0.125deg_PT1H_';
YearFile = sprintf('%04d', year(AnalysedDate));
MonthFile =  sprintf('%02d', month(AnalysedDate));
DayFile = sprintf('%02d', day(AnalysedDate));
HourFile = '08_R';
YearFilePrev = sprintf('%04d', year(previousDate));
MonthFilePrev = sprintf('%02d', month(previousDate));
DayFilePrev = sprintf('%02d', day(previousDate));
EndFileName = 'T18_14.nc';

AnalysedFile_Ix = [StartFileName,YearFile,MonthFile,DayFile,HourFile,...
    YearFilePrev,MonthFilePrev,DayFilePrev,EndFileName];

% Find the index of the file that matches the date pattern
Ix = find(strcmp(FileNames,AnalysedFile_Ix));
AnalisedFile = cell2mat(FileNames(Ix));

if isempty(AnalisedFile)
    HeatFlux(:,:,f) = NaN(xl,yl);
     [StartFileName, DayFile, MonthFile, YearFile]
    continue;
end

% Change to data folder and read variables from NetCDF file
cd(DataFolder)
Dair = ncread(AnalisedFile,'air_density'); 
v = ncread(AnalisedFile,'northward_wind');
u = ncread(AnalisedFile,'eastward_wind');
lon_wind = ncread(AnalisedFile,'lon');
lat_wind = ncread(AnalisedFile,'lat');

% Select data subset limited to specified region coordinates
lon_mask = (lon_wind >= region_box(1,1)) & ...
(lon_wind <= region_box(1,2));
lat_mask = (lat_wind >= region_box(2,1)) & ...
(lat_wind <= region_box(2,2));

v = v(lon_mask,lat_mask);
u = u(lon_mask,lat_mask);
V = sqrt((u.^2+v.^2)); 

lon_wind = lon_wind(lon_mask);
lat_wind = lat_wind(lat_mask);
[xv,yv] = meshgrid(lon_wind,lat_wind);

clear u v

% --------------------------------------------------------------------- %
% --------------------------- Input Files: ---------------------------- %
% ------------- Dew Point Temperature, Air Temperature and ------------ %
% ----------------------- Air Pressure (ERA5) ------------------------- %
% --------------------------------------------------------------------- %

DataFolder = ""; % Path to the folder containing the files with air 
% density, dew point temperature and air temperature data
FolderContents = dir(DataFolder); 
FileNames = cell(0); 

for i = 1:numel(FolderContents)
    fileName = FolderContents(i).name;
        if strcmp(fileName, '.') || strcmp(fileName, '..')
        continue;
    end
    
    FileNames{end+1} = fileName; 
end
FileNames = FileNames';

% Match file by date from filename structure
AnalysedDate = Time_Interval(f); 
StartFileName = 'era';
YearFile = sprintf('%04d', year(AnalysedDate));
MonthFile =  sprintf('m%02d', month(AnalysedDate));
DayFile = sprintf('%02d', day(AnalysedDate));
FormatFile = 'nc';

AnalysedFile_Ix = strjoin({StartFileName, YearFile, MonthFile,DayFile, ...
    FormatFile}, '.');

% Find the index of the file that matches the date pattern
Ix = find(strcmp(FileNames,AnalysedFile_Ix));
AnalisedFile = cell2mat(FileNames(Ix));

lon_e = ncread(AnalisedFile,'longitude'); 
lat_e = ncread(AnalisedFile,'latitude'); 
time_e = ncread(AnalisedFile,'time'); 
d2m = ncread(AnalisedFile, 'd2m'); 
Tar = ncread(AnalisedFile, 't2m'); 
surfacep=ncread(AnalisedFile,'sp'); 

% Select data subset limited to specified region coordinates
lon_mask = (lon_e >= region_box(1,1)) & (lon_e <= region_box(1,2));
lat_mask = (lat_e >= region_box(2,1)) & (lat_e <= region_box(2,2));

d2m = d2m(lon_mask, lat_mask);
Tar = Tar(lon_mask, lat_mask);
surfacep = surfacep(lon_mask, lat_mask);

[xe,ye] = meshgrid(lon_e,lat_e);

% --------------------------------------------------------------------- %
% ----------------------- Data standardization ------------------------ %
% --------------------------------------------------------------------- %
surfacep = rot90(surfacep);
d2m = rot90(d2m);
Tar = rot90(Tar);
V = flipud(rot90(V));
sst = flipud(rot90(sst));

% --------------------------------------------------------------------- %
% --------------------------- Interpolation --------------------------- %
% --------------------------------------------------------------------- %
surfacepq = interp2(xe,flipud(ye),surfacep,xs,ys,'linear');   
Tarq = interp2(xe,flipud(ye),Tar,xs,ys,'linear');
vq = interp2(xv,yv,V,xs,ys,'linear');
d2mq = interp2(xe,flipud(ye),d2m,xs,ys,'linear');

% --------------------------------------------------------------------- %
% ------------------------ Humidity Calculation ----------------------- %
% --------------------------------------------------------------------- %

qs = zeros(length(lat_s), length(lon_s));
qa = zeros(length(lat_s), length(lon_s));

for i = 1:length(lat_s)
    for j = 1:length(lon_s)
        sst_celsius = sst(i,j)-273.15;
        surface_pressure = surfacepq(i,j);
        dewpoint_celsius = d2mq(i,j) - 273.15; % Celsius
        
        % Specific humidity at sea surface (qs)
        A =  0.7859 + (0.03477 .* sst_celsius);
        B = 1+(0.00412 .* sst_celsius);
        saturation_vapor_pressure = 10.^(A./B)*100;  % Pa (converted by 
        % *100)
        qs_numerator = (saturation_vapor_pressure .* 0.622);
        qs_denominator = (surface_pressure - ...
            (0.378 .* saturation_vapor_pressure));
        qs_Ts = qs_numerator/qs_denominator;
        qs (i,j) = 0.98 .* qs_Ts .*1000; % g/kg 
    
        % Specific humidity of the air (qa)
        e1 = 17.67 .* dewpoint_celsius;
        e2 = dewpoint_celsius+243.5;
        vapor_pressure = (6.112 .* exp(e1./e2)).*100; % Pa
        qa_numerator = 0.622 .* vapor_pressure;
        qa_denominator = surface_pressure - (0.378 .* vapor_pressure);
        qa (i,j) = (qa_numerator./qa_denominator).*1000; % g/kg 
               
    end
end   

clear E Es saturation_vapor_pressure i j k air_temp sst_celsius ...
dewpoint_celsius SEALEVELP E cE Es RH qa_denominator qa_numerator e1 e2 ...
vapor_pressure qs_numerator qs_denominator qs_Ts A B

% --------------------------------------------------------------------- %
% ---------------------- Coefficients Calculation --------------------- %
% --------------------------------------------------------------------- %

% Coefficient calculation based on Isemer and Hasse (1987)
Ce_table = [0.06 0.26 0.63 1.15 1.44 1.78 2.19;...
    0.19 0.58 0.97 1.17 1.26 1.46 1.75; ...
    0.6 1.02 1.18 1.25 1.27 1.37 1.56; ...
    0.92 1.18 1.29 1.33 1.37 1.44 1.56; ...
    1.21 1.37 1.40 1.43 1.46 1.51 1.60; ...
    1.38 1.46 1.52 1.57 1.58 1.62 1.69; ...
    1.51 1.56 1.59 1.62 1.62 1.62 1.68; ...
    1.57 1.60 1.61 1.62 1.63 1.64 1.68; ...
    1.62 1.62 1.62 1.62 1.62 1.62 1.62];

% Required parameters
wind_speed = vq; 
sea_surface_temp = (sst)- 273.15; % Celsius
air_temp = Tarq - 273.15; % Celsius
delta_temp = air_temp-sea_surface_temp;

latent_heat_vaporization = ((2.501 - (0.00237 .* sea_surface_temp)) ...
    * 1e6); % J/kg

mean_air_density = nanmean(Dair); 
air_density_field = ones(size(sst)) .* mean_air_density;

% Wind speed categories (rows for Ce_table)
wind_category = nan(size(wind_speed));
wind_category(wind_speed <= 3) = 1;
wind_category(wind_speed > 3 & wind_speed <= 6) = 2;
wind_category(wind_speed > 6 & wind_speed <= 9) = 3;
wind_category(wind_speed > 9 & wind_speed <= 12) = 4;
wind_category(wind_speed > 12 & wind_speed <= 15) = 5;
wind_category(wind_speed > 15 & wind_speed <= 20) = 6;
wind_category(wind_speed > 20 & wind_speed <= 25) = 7;
wind_category(wind_speed > 25 & wind_speed <= 30) = 8;

% Temperature difference categories (columns for Ce_table)
temp_diff_category = nan(size(delta_temp));
temp_diff_category(delta_temp >= 5) = 1;
temp_diff_category(delta_temp >= 1 & delta_temp < 5) = 2;
temp_diff_category(delta_temp < 1 & delta_temp >= 0.2) = 3;
temp_diff_category(delta_temp < 0.2 & delta_temp >= -0.2) = 4;
temp_diff_category(delta_temp < -0.2 & delta_temp >= -1) = 5;
temp_diff_category(delta_temp < -1 & delta_temp > -5) = 6;
temp_diff_category(delta_temp <= -5) = 7;

% Initialize coefficient matrix
coeff_Ce = nan(size(sst));

% Assign coefficients based on wind and temperature categories
for i = 1:size(temp_diff_category,1)
    for j = 1:size(temp_diff_category,2)
        if isnan(wind_category(i,j)) || isnan(temp_diff_category(i,j))
            coeff_Ce(i,j) = NaN;   
        else
            coeff_Ce(i,j) = Ce_table(wind_category(i,j), ...
                temp_diff_category(i,j));
        end 
    end
end 

clearvars -except sst vq Tarq qs qa Darqq latent_heat_vaporization coeff_Ce


% --------------------------------------------------------------------- %
% ----------------------- Heat Flux Calculation ----------------------- %
% --------------------------------------------------------------------- %

% Parameters for flux calculation
coeff_Ce = coeff_Ce * 1e-3;
coeff_Ch = 0.94 * coeff_Ce;
Cp = 1004; % Specific heat capacity of air (J/kg/K)
Cp = ones(size(sst)) * Cp;

qs_kgkg = qs / 1000;
qa_kgkg = qa / 1000;

if strcmp(WHICH_FLUX, 'LAT')
    latent_heat_flux = air_density_field .* latent_heat_vaporization .* ...
        coeff_Ce .* wind_speed .* (qs_kgkg - qa_kgkg);
    HeatFlux(:,:,f) = latent_heat_flux;
    save('LHF_2019_2021.mat','-v7.3')

elseif strcmp(WHICH_FLUX, 'SENS')
    sensible_heat_flux = air_density_field .* Cp .* coeff_Ch .* ...
        wind_speed .* (sea_surface_temp - air_temp);
    HeatFlux(:,:,f) = sensible_heat_flux;
    save('SHF_2019_2021.mat','-v7.3')

else
    error('Invalid WHICH_FLUX parameter. Use ''LAT'' or ''SENS''.');
end

end



%%


