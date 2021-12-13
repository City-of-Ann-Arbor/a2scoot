clc;clear;close all;
%inputs
m= 70; %mass of rider (kg)
PointA= 'NC48 Shortcut, Ann Arbor, MI 48109';   %trip start address
PointB= '2101 Bonisteel Boulevard, Ann Arbor, MI 48109'; %trip end address

% Parameters
L = 0.97;       % Length from front to back wheel (m)
g = 9.81;       % Gravity (m/s^2)
mu = 0.004;     % Friction coeff for bike tire on asphalt
rho = 1.225;    % Density of air (kg/m^3)
M = 15+m;       % Total mass of scooter with rider (kg)  **if no m assume default value
Cd = 1.1; %1.5  % Aerodynamic Drag coefficient
A = 0.6;        % Cross-sectional Area of a person (m^2)
v_c = 11;       % Assumed constant maximum velocity from data (m/s)

% Wheel parameters
r = .1;         % Wheel radius (m)
Mw = 4;         % Wheel mass (kg)
Iw = Mw*L^2/2;  % Wheel Inertia (both) (kg m^2)

%% Get Distance and elevation (from elevations_and_dist)
        %output distTotal and elevationsTotal

        [polyline, durationStr] = getPolyline(PointA, PointB);
        [elevations, resolution, latitudes, longitudes] = getElevationsPolyline(polyline, 'key', 'YOUR_API_KEY');

        SAMPLES = 10;
        API_KEY = 'YOUR_API_KEY';

        cumulativeDist = 0;
        cumElevation = 0;
        distTotal = zeros(1, 40*(length(elevations)-1)); % preallocate to append distances to array to get total dist
        elevationsTotal = zeros(1, 40*(length(elevations)-1));
        trip_time = durationStr; %min
      for segmentNum = 1:(length(elevations) - 1)
            % loop all latlong coordinates with getElevationsTwo to get coarser data
            % points for elevation 
            %  Coord1 = [latitudes(segmentNum), longitudes(segmentNum)];
            %  Coord2 = [latitudes(segmentNum + 1), longitudes(segmentNum + 1)];
            lat = linspace(latitudes(segmentNum), latitudes(segmentNum+1), 40); % lat of point a and b
            lon = linspace(longitudes(segmentNum), longitudes(segmentNum+1), 40); % long of point a and b
            try
            % Calculate distance between two coordinates 
            radius=6371; %radius of earth in km 
            lat1=lat(1)*pi/180;
            lat2=lat(40)*pi/180;
            lon1=lon(1)*pi/180;
            lon2=lon(40)*pi/180;
            deltaLat=lat2-lat1;
            deltaLon=lon2-lon1;
            a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2;
            c=2*atan2(sqrt(a),sqrt(1-a));
            d1km=radius*c;    %Haversine distance
            d1m = d1km * 1000; % convert to meters 
            coordToCoordDist = linspace(0, d1m, 40); % Distance between coord points in
            %m w/ step size of 40
            % Get elevation between two coordinates 
            [elevation, resolution] = getElevationsTwo(lat, lon, 'key', 'YOUR_API_KEY' );

            % For loop to get total distance and total elevation
            for elementInDistArray = [1+40*(segmentNum-1)]:[segmentNum*40]
                distTotal(elementInDistArray)= coordToCoordDist((elementInDistArray - (40*[segmentNum-1])));
                elevationsTotal(elementInDistArray)= elevation(elementInDistArray - (40*[segmentNum-1]));
                if segmentNum > 1
                    distTotal(elementInDistArray) = distTotal(elementInDistArray) + distTotal(40*(segmentNum-1));
                end
            end

            % Calculate differences in elevation (meters)
            format long g
            elevation_m(segmentNum) = (elevations(segmentNum+1) - elevations(segmentNum)) ;
            % Calculate grade 
            format long g
            grade_decimal(segmentNum) = elevation_m(segmentNum) / distance_m(segmentNum);
            format long g
            grade_percent(segmentNum) = grade_decimal(segmentNum) * 100;
            end

        end
 
%% Calculate speed and accel (state space)

%Get Time, assuming constant velocity = 30km/hr about 8.5 m/s
[d,i] = unique(distTotal);       %distance vector
h = elevationsTotal(i); %elevation vector

n=2000;
accel = zeros(n,1) ;
velocity = zeros(n,1) ;
pos = zeros(n,1);
timestep=1;

Fmotor=0;
Fad = 0;
theta=atan(diff(h)./diff(d));
%index = find(isnan(theta));
%theta([index]) = 0;
velocity(1)= v_c;
kp=60;
ki = 13;

error_accum = 0;
for i = 1:n-1
    error = velocity(i) - v_c;
    error_accum = error_accum + error;
    Fmotor(i) = -kp*error - ki*error_accum;
    Pr(i) = Fmotor(i)*velocity(i);
    if Pr(i) > 500
        %set torque limit
        Fmotor(i) = 500./velocity(i);
        Pr(i) = 500;
        
        %anti wind up scheme for error accum
        error_accum = error_accum - error;
        %error = velocity(i) - v_c;
        
    end
    if Fmotor(i) > 160
            Fmotor(i) = 160;
            Pr(i) = Fmotor(i)*velocity(i);
    end
    if Fmotor(i) < -160
            Fmotor(i) = -160;
            Pr(i) = Fmotor(i)*velocity(i);
    end
    %large drag
    Fad(i) = (0.5*rho*A*Cd).*velocity(i)^2;  %unit [N]
    theta1(i) = interp1(d(1:end-1),theta,pos(i));
    Frr(i) = mu*M*g*cos(theta1(i));  %unit [N]
    Wx(i) = M*g*sin(theta1(i)); %unit [N]
    
    accel(i) = (Fmotor(i)-Frr(i)-Fad(i)-Wx(i))./M;
    velocity(i+1)=timestep*accel(i)+velocity(i);
    pos(i+1)=timestep*velocity(i)+pos(i);
end

index = find(isnan(Fmotor));
Fmotor([index]) = [];
index = find(isnan(velocity));
velocity([index]) = [];

%
for i = 1:numel(Fmotor)
    eff(i) = gradient(Fmotor(i), velocity(i));
end

for z = 1:numel(Pr)
    if Pr(z) > 0
        Pe(z) = Pr(z)./eff(z);
    elseif Pr(z) < 0 
        Pe(z) = Pr(z).*eff(z);
    elseif Pr(z) == 0
        Pe(z) = Pr(z);
    end
end
t = ones(1,numel(Pe));

E = 0;
for i = 1:numel(Pe)
    E = E + Pe(i)*t(i);
end

%% Calculate Energy and Emissions
E = E/3600; %Wh
trip_dist = (distTotal(end)./1609.34); %mi
E_per_mi = E./trip_dist;  %Wh/mi

%convert to emissions
% conversion source: https://www.epa.gov/energy/greenhouse-gases-equivalencies-calculator-calculations-and-references
E_co2 = E*0.8491032;  %gCO2
E_co2_mi = E_per_mi*0.8491032; %gCO2/pmt (passenger-mile travelled)


%% Print outputs ---------------------------------------------------------
%total trip distance = trip_dist
%total trip time = trip_time  %trip_time = durationStr
%total trip energy per mile = E_per_mi
%total trip energy = E
%total trip emissions per mile = E_co2_mi
%total trip emissions = E_co2

fprintf('The Time required to complete this %4.2f mi long trip is %s.', trip_dist, trip_time)
fprintf('\nThe Energy required to complete this trip is approximately %4.2f Watt-hours per mile, or %4.2f Watt-hours total.', E_per_mi, E)
fprintf('\nThe Emissions generated to complete this trip is approximately %4.2f gCO2/pmt, or %4.2f gCO2 total.\n', E_co2_mi, E_co2)

%% function to find efficiency from ansys map
%force should be in N and vel should be in m/s
function eff = gradient(force, vel)
force = abs(force);
vel = abs(vel);
mapped_force = interp1([0 160], [10 310], force);
mapped_vel = interp1([0 12], [200 7000], vel);

rpm_range1 = [200 400];
torque_range1 = [10 30 100 310];

if mapped_vel >= rpm_range1(1) && mapped_vel < rpm_range1(2)
    if mapped_force >= torque_range1(1) && mapped_force < torque_range1(2)
        eff = .92;
    elseif mapped_force >= torque_range1(2) && mapped_force < torque_range1(3)
        eff = .88;
    elseif mapped_force >= torque_range1(3) && mapped_force <= torque_range1(4)
        eff = .8;
    else
        eff = .8;
    end
end

rpm_range2 = [400 600];
torque_range2 = [10 20 50 80 200 310];

if mapped_vel >= rpm_range2(1) && mapped_vel < rpm_range2(2)
    if mapped_force >= torque_range2(1) && mapped_force < torque_range2(2)
        eff = .92;
    elseif mapped_force >= torque_range2(2) && mapped_force < torque_range2(3)
        eff = .95;
    elseif mapped_force >= torque_range2(3) && mapped_force < torque_range2(4)
        eff = .92;
    elseif mapped_force >= torque_range2(4) && mapped_force < torque_range2(5)
        eff = .88;
    elseif mapped_force >= torque_range2(5) && mapped_force <= torque_range2(6)
        eff = .8;
    else
        eff = .8;
    end
end

rpm_range3 = [600 800];
torque_range3 = [10 20 60 100 280 310];

if mapped_vel >= rpm_range3(1) && mapped_vel < rpm_range3(2)
    if mapped_force >= torque_range3(1) && mapped_force < torque_range3(2)
        eff = .92;
    elseif mapped_force >= torque_range3(2) && mapped_force < torque_range3(3)
        eff = .95;
    elseif mapped_force >= torque_range3(3) && mapped_force < torque_range3(4)
        eff = .92;
    elseif mapped_force >= torque_range3(4) && mapped_force < torque_range3(5)
        eff = .88;
    elseif mapped_force >= torque_range3(5)
        eff = .8;
    end
end

rpm_range4 = [800 1000];
torque_range4 = [10 20 130 180 310];

if mapped_vel >= rpm_range4(1) && mapped_vel < rpm_range4(2)
    if mapped_force >= torque_range4(1) && mapped_force < torque_range4(2)
        eff = .88;
    elseif mapped_force >= torque_range4(2) && mapped_force < torque_range4(3)
        eff = .95;
    elseif mapped_force >= torque_range4(3) && mapped_force < torque_range4(4)
        eff = .92;
    elseif mapped_force >= torque_range4(4)
        eff = .88;
    end
end

rpm_range5 = [1000 2000];
torque_range5 = [10 30 200 250 310];

if mapped_vel >= rpm_range5(1) && mapped_vel < rpm_range5(2)
    if mapped_force >= torque_range5(1) && mapped_force < torque_range5(2)
        eff = .88;
    elseif mapped_force >= torque_range5(2) && mapped_force < torque_range5(3)
        eff = .95;
    elseif mapped_force >= torque_range5(3) && mapped_force < torque_range5(4)
        eff = .92;
    elseif mapped_force >= torque_range5(4)
        eff = .88;
    end
end

rpm_range6 = [2000 3000];
torque_range6 = [10 20 40];

if mapped_vel >= rpm_range6(1) && mapped_vel < rpm_range6(2)
    if mapped_force >= torque_range6(1) && mapped_force < torque_range6(2)
        eff = .8;
    elseif mapped_force >= torque_range6(2) && mapped_force < torque_range6(3)
        eff = .88;
    elseif mapped_force >= torque_range6(3)
        eff = .95;
    end
end

rpm_range7 = [3000 4000];
torque_range7 = [10 20 60];

if mapped_vel >= rpm_range7(1) && mapped_vel < rpm_range7(2)
    if mapped_force >= torque_range7(1) && mapped_force < torque_range7(2)
        eff = .8;
    elseif mapped_force >= torque_range7(2) && mapped_force < torque_range7(3)
        eff = .88;
    elseif mapped_force >= torque_range7(3)
        eff = .92;
    end
end

rpm_range8 = [4000 5000];
torque_range8 = [10 20];

if mapped_vel >= rpm_range8(1) && mapped_vel < rpm_range8(2)
    if mapped_force >= torque_range8(1) && mapped_force < torque_range8(2)
        eff = .8;
    elseif mapped_force >= torque_range8(2)
        eff = .88;
    end
end

rpm_range9 = [5000 7000];
torque_range9 = [10];

if mapped_vel >= rpm_range9(1)
    if mapped_force >= torque_range9(1)
        eff = .8;
    end
end
end

%% function get Polyline
function [polyline, durationStr] = getPolyline(initialInput, finalInput)
% takes initial and final addresses and outputs their place IDs via google
% place search api 
% https://developers.google.com/maps/documentation/places/web-service/search#FindPlaceRequests

% Check parameters 
%prompt_1 = 'Enter address 1'; 
%initialInput = input(prompt_1, 's'); % makes sure input is a string
initialInput = regexprep(initialInput, ' ', '%'); % replaces space with & 
inputType = '&inputtype=textquery';
key = '&key=YOUR_API_KEY';
locationBias = '&locationbias=circle:500@42.2808,83.7430' ;
%centered at ann arbor's coordinates with radius of 500 m; doesn't work if
%user's ip address is not in AA, so you need city and state if so

website = 'https://maps.googleapis.com/maps/api/place/findplacefromtext/json?input=';
url = [website, initialInput, inputType, locationBias, key];
str = urlread(url);
str;
data = jsondecode(str);

% Parse results
status = data.status;
switch status
case 'OK'
  initialPlaceID = [data.candidates.place_id];
case 'INVALID_REQUEST'
  error('Google Maps API request was malformed or Invalid Resolution');
case 'OVER_QUERY_LIMIT'
  error('Google Maps API requestor has exceeded quota');
case 'REQUEST_DENIED'
  error('Google Maps API did not complete the request (invalid sensor parameter?)');
case 'UNKNOWN_ERROR'
  error('Google Maps API: an unknown error.');
end
% Repeat for address 2 
%prompt_2 = 'Enter address 2'; 
%finalInput = input(prompt_2, 's'); % makes sure input is a string
finalInput = regexprep(finalInput, ' ', '%'); % replaces space with & 

inputType = '&inputtype=textquery';

locationBias = '&locationbias=circle:500@42.2808,83.7430' ;

website = 'https://maps.googleapis.com/maps/api/place/findplacefromtext/json?input=';
url = [website, finalInput, inputType, locationBias, key];
str2 = urlread(url);
data2 = jsondecode(str2);

% Parse results
status = data2.status;
switch status
case 'OK'
  finalPlaceID = [data2.candidates.place_id];
case 'INVALID_REQUEST'
  error('Google Maps API request was malformed or Invalid Resolution');
case 'OVER_QUERY_LIMIT'
  error('Google Maps API requestor has exceeded quota');
case 'REQUEST_DENIED'
  error('Google Maps API did not complete the request (invalid sensor parameter?)');
case 'UNKNOWN_ERROR'
  error('Google Maps API: an unknown error.');
end

% Get directions 
% uses intial and final place IDs to get the directions in an encoded
% polyline
directionSite = 'https://maps.googleapis.com/maps/api/directions/json?';
origin = 'origin=place_id:';
destination = '&destination=place_id:';
travelMode = '&mode=bicycling';
url = [directionSite, origin, initialPlaceID, destination, finalPlaceID, travelMode, key];
str3 = urlread(url);
data3 = jsondecode(str3);
polyline = data3.routes.overview_polyline.points;

%get trip duration
polyline = data3.routes.overview_polyline.points;
polylineOutput = strcat("Polyline: ", polyline);
            app.PolylineLabel.Text=polylineOutput; %outputs polyline using from and to points
            durationStr = erase(data3.routes.legs.duration.text, " mins"); % duration of trip
            app.DurationLabel.Text= strcat("Duration: ", durationStr, " minutes");
distanceTransit = 0;
end


%% function getElevationsTwo
function [elevation, resolution] = getElevationsTwo(latitude, longitude, varargin)
%GETELEVATIONS queries Google Maps API webservice for ground elevations.
%
%   elevation = GETELEVATIONS(latitude, longitude) returns ground 
%   elevation for latitude and longitude arrays.
%
%   [elevation, resolution] = getElevations(latitude, longitude, ...
%     'key', 'YOUR_API_KEY' );
%   is an example of a call passing additional attributes to Google Maps 
%   API webservice, and capturing also the resolution of the data. 
%   See https://developers.google.com/maps/documentation/elevation/
%   for details.
%
% Author: Jarek Tuszynski (jaroslaw.w.tuszynski@leidos.com)
% License: BSD (Berkeley Software Distribution)
% Documentation: https://developers.google.com/maps/documentation/elevation/

% process varargin
keyStr = '';
if nargin>2
  p = inputParser;
  p.addParameter('key', '', @ischar)
  p.FunctionName = 'getElevations';
  p.parse(varargin{:})
  results = p.Results;
  if ~isempty(results.key)
    keyStr = sprintf('&key=%s', results.key);
  end
end

radius=6371;
lat1=latitude(1)*pi/180;
lat2=latitude(40)*pi/180;
lon1=longitude(1)*pi/180;
lon2=longitude(40)*pi/180;
deltaLat=lat2-lat1;
deltaLon=lon2-lon1;
a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2;
c=2*atan2(sqrt(a),sqrt(1-a));
d1km=radius*c;    %Haversine distance
d1m = d1km * 1000; 
%x=deltaLon*cos((lat1+lat2)/2);
%y=deltaLat;
%d2km=radius*sqrt(x*x + y*y); %Pythagoran distance
d1mi = d1km * 0.621371; % 1 km = 0.621371 mi, distance in miles (more accurate for local distances)
%d2mi = d2km * 0.621371; % 1 km = 0.621371 mi, distance in miles 
dist= linspace(0, d1m, 40);
% Check inputs
nPos = numel(latitude);
assert(nPos>0, 'Latitude and longitude inputs can not be empty')
assert(nPos==numel(longitude), 'Latitude and longitude inputs are not of the same length')
assert(min(latitude(:)) >= -90 && max(latitude(:)) <= 90, 'Latitudes has to be between -90 and 90')
assert(min(longitude(:))>=-180 && max(longitude(:))<=180, 'Longitude has to be between -180 and 180')

% Query Google
elevation  = zeros(size(latitude))*nan;
resolution = zeros(size(latitude))*nan;
batch = [1:50:nPos nPos+1]; % group in batches of 50
for iBatch=2:length(batch)
  idx = batch(iBatch-1):batch(iBatch)-1;
  coord = '';
  for k = 1:length(idx)
    coord = sprintf('%s%9.6f,%9.6f|',coord,latitude(idx(k)),longitude(idx(k)));
  end
  
  % create query string and run a query
  website = 'https://maps.googleapis.com/maps/api/elevation/xml?locations=';
  url = [website, coord(1:end-1), keyStr];
  str = urlread(url);
  
  % Parse results
  status = regexp(str, '<status>([^<]*)<\/status>', 'tokens');
  switch status{1}{1}
    case 'OK'
      res = regexp(str, '<elevation>([^<]*)<\/elevation>', 'tokens');
      elevation(idx) = cellfun(@str2double,res);
      if nargout>1
       res = regexp(str, '<resolution>([^<]*)<\/resolution>', 'tokens');
       resolution(idx) = cellfun(@str2double,res);
      end
    case 'INVALID_REQUEST'
      error('Google Maps API request was malformed');
    case 'OVER_QUERY_LIMIT'
      error('Google Maps API requestor has exceeded quota');
    case 'REQUEST_DENIED'
      error('Google Maps API did not complete the request (invalid sensor parameter?)');
    case 'UNKNOWN_ERROR'
      error('Google Maps API: an unknown error.');
  end
end
end
