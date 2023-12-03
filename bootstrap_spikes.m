% Bootstrap analysis for significance testing of direction tuning curves
% Wei Wei 2014-09-11

% load spike count data from ans.spikeCounts (the output of loadClampexData
d1 = SortedPeakData.responses{1,4}(1:8,:);
d2=reshape(d1,size(d1,1)*size(d1,2),1);

% specify data directions
data_directions = [0:45:315];
rad_dirs = degToRad(data_directions);

% resample n times
n = 5000;
[bootstat,bootsam] = bootstrp(n,[],d2);

bootsam
% for example the first resampled data is: d2(bootsam(:,1))
dsi_values = [];
theta = [];
radius = [];
for i = 1:n
    data_resam = d2(bootsam(:,i)); % 36x1 
    data_resam_matrix = reshape(data_resam,size(d1,1), size(d1,2));
    mean_measurements_list=mean(data_resam_matrix')';
    norm_resps = (mean_measurements_list)/sum(mean_measurements_list);
    [x_coords,y_coords] = pol2cart(rad_dirs,norm_resps');
vector_sum = [sum(x_coords) sum(y_coords)];
[th,r]=cart2pol(vector_sum(1), vector_sum(2));
theta = [theta;th];
radius = [radius;r];
% find dsi value
if th >= 0
deg = th/pi*180;
else
    deg_raw = th/pi*180;
    deg = 360 + deg_raw;
end
%find the prefered direction, use 30 if there're 12 directions; 45 for 8
%dirs
n1 = floor(deg/45);
n2= n1+1;
a= deg - n1*45;
b = n2*45 - deg;
if a>=b
    pref_dir = n2*45;
else 
    pref_dir = n1*45;
end
%note: need to account for 360 vs 330 situation:
if pref_dir == 360
    pref_dir = 0;
end
% calculate DSI
null_dir = mod(pref_dir+180,360); % make sure the null dir is less than 360
null_ind = find(abs((data_directions-null_dir)) < eps('single')); %hack to find match
pref_ind = find(abs((data_directions-pref_dir)) < eps('single')); %hack to find match
null_response = mean_measurements_list(null_ind);
pref_response = mean_measurements_list(pref_ind);
DSI = (pref_response - null_response)/(pref_response + null_response);
dsi_values = [dsi_values;DSI];
end
% Find the range of confidence intervals 
ci_dsi = prctile(dsi_values, [2.5 50 97.5]);
ci_theta = prctile(theta, [2.5 50 97.5]);
ci_radius = prctile(radius, [2.5 50 97.5]);

% Calculate the dsi theta, and vs of the raw data
[th,r]=cart2pol(cell.meanVectorSum(1),cell.meanVectorSum(2))

% save('C:\Users\wei\Documents\Work\data\light response exp in 2PT rig\091201 recordings\09d01014.mat')

%find the direction of the vector sum
if th >= 0
deg = th/pi*180;
else
    deg_raw = th/pi*180;
    deg = 360 + deg_raw;
end
%find the prefered direction, use 30 if there're 12 directions; 45 for 8
%dirs
n1 = floor(deg/30);
n2= n1+1;
a= deg - n1*30;
b = n2*30 - deg;
if a>=b
    pref_dir = n2*30;
else 
    pref_dir = n1*30;
end
%note: need to account for 360 vs 330 situation:
if pref_dir == 360
    pref_dir = 0;
end
% calculate DSI
null_dir = mod(pref_dir+180,360); % make sure the null dir is less than 360
null_ind = find(abs((cell.directions-null_dir)) < eps('single')); %hack to find match
pref_ind = find(abs((cell.directions-pref_dir)) < eps('single')); %hack to find match
null_response = cell.meanSpikeCounts(null_ind);
pref_response = cell.meanSpikeCounts(pref_ind);
DSI = (pref_response - null_response)/(pref_response + null_response);

%plot histogram of boostraped dsi values
figure(1)
hist(dsi_values,100)
%find confidence interval 1-alpha
% alpha = 0.05;
% sorted_dsi=sort(dsi_values);
% dsi_low = sorted_dsi(alpha/2*n);
% dsi_high = sorted_dsi(n-alpha/2*n);


% Compare raw value with CI from bootstrap
if DSI > ci_dsi(3)
    DSI_sig = 1 %this cell is DS
else
    DSI_sig = 0 %this cell is not DS
end

%do the same for vector sum length r
 figure(2)
hist(radius,100)   
    if r > ci_radius(3)
    r_sig = 1 %this cell is DS
else
    r_sig = 0 %this cell is not DS
end

