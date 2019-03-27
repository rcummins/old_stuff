% A script to plot error ellipses showing the correlated errors in two
% variables. In this case, the two variables are temperature (temp) and
% d18O_water. This needs to be in the same folder as your data and the
% matlab script error_ellipse.m
% 
% By Renata Cummins 2013

clear all

test_definitions;

% load processed CIDS data 
CIDS_data = load('fig_data_ellipses','Samples');
S = CIDS_data.Samples;
num_of_samps = size(S,1);

% prepare the figure
clf;
hold;

% initalize the counters of the different sample types
cem_counter = 1;
mic_counter = 1;
b_counter = 1;
r_counter = 1;
bry_counter = 1;
cri_counter = 1;
gas_counter = 1;
ost_counter = 1;

% calculate and plot the error ellipse for each data point
for i = 1:num_of_samps
    
    cov_mat = cov(S{i}{10}(:,temp_aq),S{i}{10}(:,d18O_water_aq));
    cov_mat(1,1) = cov_mat(1,1)/S{i}{12}(aqs);
    cov_mat(2,2) = cov_mat(2,2)/S{i}{12}(aqs);
    cov_mat(1,2) = cov_mat(1,2)/S{i}{12}(aqs);
    cov_mat(2,1) = cov_mat(2,1)/S{i}{12}(aqs);

    position = [S{i}{12}(temp) S{i}{12}(d18O_water)];
    
    if S{i}{12}(sample_type) == 7;
        cem_d18Ow(cem_counter) = position(2);
        cem_temp(cem_counter)  = position(1);
        cem_counter = cem_counter +1;
    elseif S{i}{12}(sample_type) == 5;
        mic_d18Ow(mic_counter) = position(2);
        mic_temp(mic_counter)  = position(1);
        mic_counter = mic_counter +1;
    elseif S{i}{12}(sample_type) == 19;
        b_d18Ow(b_counter) = position(2);
        b_temp(b_counter)  = position(1);
        b_counter = b_counter +1;
    elseif S{i}{12}(sample_type) == 23;
        r_d18Ow(r_counter) = position(2);
        r_temp(r_counter)  = position(1);
        r_counter = r_counter +1;
    elseif S{i}{12}(sample_type) == 18;
        bry_d18Ow(bry_counter) = position(2);
        bry_temp(bry_counter)  = position(1);
        bry_counter = bry_counter +1;
    elseif S{i}{12}(sample_type) == 2;
        cri_d18Ow(cri_counter) = position(2);
        cri_temp(cri_counter)  = position(1);
        cri_counter = cri_counter +1;
    elseif S{i}{12}(sample_type) == 15;
        gas_d18Ow(gas_counter) = position(2);
        gas_temp(gas_counter)  = position(1);
        gas_counter = gas_counter +1;
    elseif S{i}{12}(sample_type) == 16;
        ost_d18Ow(ost_counter) = position(2);
        ost_temp(ost_counter)  = position(1);
        ost_counter = ost_counter +1;
    end

    error_ellipse(cov_mat,position)

end

% plot the points in the middle of the ellipses, and label the axes
plot(cem_temp, cem_d18Ow, 'd', 'MarkerEdgeColor', [1 0.502 0])
plot(mic_temp, mic_d18Ow, 's', 'MarkerEdgeColor', [1 0.502 0])
plot(b_temp, b_d18Ow, 'ok', 'MarkerFaceColor','k') 
plot(r_temp, r_d18Ow, 'vk', 'MarkerFaceColor','k')
plot(bry_temp, bry_d18Ow, 'dk', 'MarkerFaceColor','k') 
plot(cri_temp, cri_d18Ow, '*k', 'MarkerFaceColor','k') 
plot(gas_temp, gas_d18Ow, 'xk', 'MarkerFaceColor','k') 
plot(ost_temp, ost_d18Ow, 'sk', 'MarkerFaceColor','k') 
%xlabel('Temperature (°C)')
%ylabel('d18O water (VSMOW)')

