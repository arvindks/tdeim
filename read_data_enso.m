function [ Lat, Lon, time, mask, X ] = read_data_enso( filename, maskname )
%   read_data Read in ENSO data
%   Return meshgrid of lat,long,times and sea-surface temp

X = ncread(filename,'sst');
time = ncread(filename,'time'); %days

lat = ncread(filename,'lat');
lon = ncread(filename,'lon');

mask = ncread(maskname,'mask');
mask = mask(:,:,1);
[Lat,Lon] = meshgrid(lat,lon);

end

