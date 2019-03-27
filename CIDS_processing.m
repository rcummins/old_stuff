% Script to import CIDS data into Matlab
%
% First I opened the CIDS in excel and saved them as a csv. Then I deleted
% all the columns to the right of the last column of voltages you copy from
% isodat. Then I selected all the columns containing the voltages you copy
% from isodat, clicked on Format in the top menu, and then clicked on
% Cells. I changed the format to 0.000 which is an option in the custom
% menu. This makes sure the full decimal accuracy is preserved in the csv
% when you import it into matlab. Then I saved the CIDS again. I saved each
% CIDS with a different number in the name, and then added it to the if 
% else-if statement below % import data. 
% 
% If the CIDS is from the 01 autoline, sometimes it has the string 
%'Rx @90C, 15 min GC, -20C' in the cell next to 'Notes'. This messes up the
% data import because of the commas, so I used find and replace to delete
% all these strings. Don't use find and replace to delete all the commas,
% because this script uses the commas in the Background report to separate
% the background voltages. 
%
% Some older CIDS also have a space before the word 'Background'. This 
% space will mess up the import, so I used find and replace to delete this 
% space. 
%
% I also deleted the rows of samples from the CIDS files that ended up
% being way off the D48 line, or anything else funky. 
% 
% I also made a file called HG_residuals.csv that has one line of headers
% and then the heated gas slope and intercept and the average standard
% residual by date and machine (old or new). 
% 
% Made by Renata Cummins in 2013

clear all

CIDS_definitions;

samp_counter = 1;
number_of_CIDS = 6; % change this if you have more or fewer than six CIDS to import
total_samples = 307; % change if you know how many samples you ran, otherwise don't worry
S = cell(total_samples,1); 

% load the heated gas line and standard residuals with corresponding dates
fileID = fopen('HG_residuals.csv');
HG_residuals = textscan(fileID, '%s %s %s %s %s %s %s %s %s %s %s %s','delimiter',',');
fclose(fileID);
fd = str2double(HG_residuals{1}(2:end));
fm = str2double(HG_residuals{2}(2:end));
fy = str2double(HG_residuals{3}(2:end));
ld = str2double(HG_residuals{4}(2:end));
lm = str2double(HG_residuals{5}(2:end));
ly = str2double(HG_residuals{6}(2:end));
mach      = HG_residuals{7}(2:end);
slope_HG     = str2double(HG_residuals{8}(2:end));
intercept_HG = str2double(HG_residuals{9}(2:end));
residual  = str2double(HG_residuals{10}(2:end));
slope_STF = str2double(HG_residuals{11}(2:end));
intercept_STF = str2double(HG_residuals{12}(2:end));

for CIDS = 1:number_of_CIDS

% import data
% if you have more or fewer than 6 CIDS, modify accordingly
  if CIDS == 1                
    fileID = fopen('CIDS_1.csv');
  elseif CIDS == 2
    fileID = fopen('CIDS_2.csv');
  elseif CIDS == 3
    fileID = fopen('CIDS_3.csv');
  elseif CIDS == 4
    fileID = fopen('CIDS_4.csv');
  elseif CIDS == 5
    fileID = fopen('CIDS_5.csv');
  elseif CIDS == 6
    fileID = fopen('CIDS_6.csv');
  end
data_raw = textscan(fileID, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','delimiter',',');
fclose(fileID);

% figure out how many samples there are
[data_raw_rows data_raw_columns] = size(data_raw{1});
number_of_samples = 0;
  for i = 1:data_raw_rows
    if (isempty(regexp(data_raw{7}{i},'d45','ONCE')) == 0)
        number_of_samples = number_of_samples + 1;
    end
  end

% pre-allocate the cell array of aquisitions, sample names, and sample summary
  for i = samp_counter : samp_counter + number_of_samples - 1
    S{i} = cell(12,1);
    for j = 1:9
        S{i}{j} = zeros(7,38);
    end
    S{i}{10} = zeros(9,110);
    S{i}{11} = cell(1,4);
    S{i}{12} = zeros(1,35);
  end

% grab the data and save it into the preallocated cell arrays
aq_counter = 1;
numbers = zeros(8,12);
date = zeros(1,6);
  for i = 1:data_raw_rows
    if (isempty(regexp(data_raw{4}{i},'pre','ONCE')) == 0)
       
        % change voltages to numbers
       for j = 1:12                                           
           numbers(1:8,j) = str2double(data_raw{j+4}(i:i+7));
       end
       
       % store voltages in array
       for j = 1:7                                            
           S{samp_counter}{aq_counter}(j,1:6) = numbers(j,7:12);    
           S{samp_counter}{aq_counter}(j,7:18) = numbers(j+1,1:12);
       end
       S{samp_counter}{11}{mach_id} = data_raw{2}{i+7};
       S{samp_counter}{11}{analyst} = data_raw{2}{i+8};
       date = datevec(datenum(data_raw{2}{i+3}));
       S{samp_counter}{10}(aq_counter,year) = date(1);
       S{samp_counter}{10}(aq_counter,month) = date(2);
       S{samp_counter}{10}(aq_counter,day) = date(3);
       S{samp_counter}{11}{type} = data_raw{2}{i+5};
       S{samp_counter}{11}{ID} = data_raw{2}{i};
       S{samp_counter}{10}(aq_counter,aq) = str2double(data_raw{3}{i});
       S{samp_counter}{10}(aq_counter,d13C) = str2double(data_raw{2}{i+9});
       S{samp_counter}{10}(aq_counter,d18O_SMOW) = str2double(data_raw{2}{i+10});
       S{samp_counter}{10}(aq_counter,ref_d13C) = str2double(data_raw{2}{i+11});
       S{samp_counter}{10}(aq_counter,ref_d18O_SMOW) = str2double(data_raw{2}{i+12});
       S{samp_counter}{10}(aq_counter,peak_center) = sscanf(data_raw{1}{i+13},'Peak Center found at [%u');
       S{samp_counter}{10}(aq_counter,background_44mV) = sscanf(data_raw{1}{i+14},'"Background: %f');
       S{samp_counter}{10}(aq_counter,background_45mV) = sscanf(data_raw{2}{i+14},'%f');
       S{samp_counter}{10}(aq_counter,background_46mV) = sscanf(data_raw{3}{i+14},'%f');
       S{samp_counter}{10}(aq_counter,background_47mV) = sscanf(data_raw{4}{i+14},'%f');
       S{samp_counter}{10}(aq_counter,background_48mV) = sscanf(data_raw{5}{i+14},'%f');
       S{samp_counter}{10}(aq_counter,background_49mV) = sscanf(data_raw{6}{i+14},'%f');
       aq_counter = aq_counter +1;
    end
    
    % check if a sample is done
    if (isempty(regexp(data_raw{7}{i},'d45','ONCE')) == 0)     
       samp_counter = samp_counter +1;
       aq_counter = 1;
    end
  end
end

% calculate d45, d46, d47, d48, d49 for each acquisition
for i = 1:total_samples
  for j = 1:9
    
    % check that the acquisition is not empty
    if S{i}{j}(1,1) ~= 0
            
    % first calculate the values for each cycle of the acquisition
    
    % calculate R values for ref gas measured before sample
    S{i}{j}(:,R45_pre) = S{i}{j}(:,pre45mV)./S{i}{j}(:,pre44mV);
    S{i}{j}(:,R46_pre) = S{i}{j}(:,pre46mV)./S{i}{j}(:,pre44mV);
    S{i}{j}(:,R47_pre) = S{i}{j}(:,pre47mV)./S{i}{j}(:,pre44mV);
    S{i}{j}(:,R48_pre) = S{i}{j}(:,pre48mV)./S{i}{j}(:,pre44mV);
    S{i}{j}(:,R49_pre) = S{i}{j}(:,pre49mV)./S{i}{j}(:,pre44mV);
    
    % calculate R values for sample
    S{i}{j}(:,R45_samp) = S{i}{j}(:,samp45mV)./S{i}{j}(:,samp44mV);
    S{i}{j}(:,R46_samp) = S{i}{j}(:,samp46mV)./S{i}{j}(:,samp44mV);
    S{i}{j}(:,R47_samp) = S{i}{j}(:,samp47mV)./S{i}{j}(:,samp44mV);
    S{i}{j}(:,R48_samp) = S{i}{j}(:,samp48mV)./S{i}{j}(:,samp44mV);
    S{i}{j}(:,R49_samp) = S{i}{j}(:,samp49mV)./S{i}{j}(:,samp44mV);
    
    % calculate R values for ref gas measured after sample
    S{i}{j}(:,R45_post) = S{i}{j}(:,post45mV)./S{i}{j}(:,post44mV);
    S{i}{j}(:,R46_post) = S{i}{j}(:,post46mV)./S{i}{j}(:,post44mV);
    S{i}{j}(:,R47_post) = S{i}{j}(:,post47mV)./S{i}{j}(:,post44mV);
    S{i}{j}(:,R48_post) = S{i}{j}(:,post48mV)./S{i}{j}(:,post44mV);
    S{i}{j}(:,R49_post) = S{i}{j}(:,post49mV)./S{i}{j}(:,post44mV);
    
    % calculate delta values for sample
    S{i}{j}(:,d45_cycle) = 1000*((S{i}{j}(:,R45_samp)./((S{i}{j}(:,R45_pre)+S{i}{j}(:,R45_post))/2)-1));
    S{i}{j}(:,d46_cycle) = 1000*((S{i}{j}(:,R46_samp)./((S{i}{j}(:,R46_pre)+S{i}{j}(:,R46_post))/2)-1));
    S{i}{j}(:,d47_cycle) = 1000*((S{i}{j}(:,R47_samp)./((S{i}{j}(:,R47_pre)+S{i}{j}(:,R47_post))/2)-1));
    S{i}{j}(:,d48_cycle) = 1000*((S{i}{j}(:,R48_samp)./((S{i}{j}(:,R48_pre)+S{i}{j}(:,R48_post))/2)-1));
    S{i}{j}(:,d49_cycle) = 1000*((S{i}{j}(:,R49_samp)./((S{i}{j}(:,R49_pre)+S{i}{j}(:,R49_post))/2)-1));
    
    % now take the mean of all cycles in each acquisition
    S{i}{10}(j,d45) = mean(S{i}{j}(:,d45_cycle));
    S{i}{10}(j,d46) = mean(S{i}{j}(:,d46_cycle));
    S{i}{10}(j,d47) = mean(S{i}{j}(:,d47_cycle));
    S{i}{10}(j,d48) = mean(S{i}{j}(:,d48_cycle));
    S{i}{10}(j,d49) = mean(S{i}{j}(:,d49_cycle));
    
    % else if the acquisition is empty, delete that row from the aq array
    else
    index = true(1,size(S{i}{10},1));
    if j > size(S{i}{10},1)
        index(size(S{i}{10},1)) = false;
    else
        index(j) = false;
    end
    S{i}{10} = S{i}{10}(index,:);
    end
    
  end
end

% pre-allocate the cell array of only samples
num_of_samps = 87;
Samples = cell(num_of_samps,1);
  for i = 1:length(Samples)
    Samples{i} = cell(12,1);
    for j = 1:9
        Samples{i}{j} = zeros(1,1);
    end
    Samples{i}{10} = zeros(9,106);
    Samples{i}{11} = cell(1,4);
    Samples{i}{12} = zeros(1,35);
  end
  sample_counter = 1;
  
% calculate d13C, d18O, d47, D47, d48, D48, d18O_water for each acquisition
for i = 1:total_samples
    
    % calculate the R13, R18, R17 of sample and ref gas
    S{i}{10}(:,R13_r) = ((S{i}{10}(:,ref_d13C)./1000)+1).*VPDB_R13;
    S{i}{10}(:,R13_s) = ((S{i}{10}(:,d13C)./1000)+1).*VPDB_R13;
    S{i}{10}(:,R18_r) = ((S{i}{10}(:,ref_d18O_SMOW)./1000)+1).*VSMOW_R18;
    S{i}{10}(:,R18_s) = ((S{i}{10}(:,d18O_SMOW)./1000)+1).*VSMOW_R18;
    S{i}{10}(:,R17_r) = ((S{i}{10}(:,R18_r)./VSMOW_R18).^beta).*VSMOW_R17;
    S{i}{10}(:,R17_s) = ((S{i}{10}(:,R18_s)./VSMOW_R18).^beta).*VSMOW_R17;
    
    % calculate abundances of 12C, 13C, 16O, 17O, 18O in sample and ref gas
    S{i}{10}(:,a12C_r) = (1+S{i}{10}(:,R13_r)).^-1;
    S{i}{10}(:,a12C_s) = (1+S{i}{10}(:,R13_s)).^-1;
    S{i}{10}(:,a13C_r) = -1.*S{i}{10}(:,a12C_r)+1;
    S{i}{10}(:,a13C_s) = -1.*S{i}{10}(:,a12C_s)+1;
    S{i}{10}(:,a16O_r) = (1+S{i}{10}(:,R18_r)+S{i}{10}(:,R17_r)).^-1;
    S{i}{10}(:,a16O_s) = (1+S{i}{10}(:,R18_s)+S{i}{10}(:,R17_s)).^-1;
    S{i}{10}(:,a17O_r) = S{i}{10}(:,R17_r).*S{i}{10}(:,a16O_r);
    S{i}{10}(:,a17O_s) = S{i}{10}(:,R17_s).*S{i}{10}(:,a16O_s);
    S{i}{10}(:,a18O_r) = S{i}{10}(:,R18_r).*S{i}{10}(:,a16O_r);
    S{i}{10}(:,a18O_s) = S{i}{10}(:,R18_s).*S{i}{10}(:,a16O_s);
    
    % calculate abundances of CO2 isotopologues in sample and ref gas
    S{i}{10}(:,a266_r) = S{i}{10}(:,a12C_r).*S{i}{10}(:,a16O_r).*S{i}{10}(:,a16O_r);
    S{i}{10}(:,a266_s) = S{i}{10}(:,a12C_s).*S{i}{10}(:,a16O_s).*S{i}{10}(:,a16O_s);
    S{i}{10}(:,a267_r) = 2.*S{i}{10}(:,a12C_r).*S{i}{10}(:,a16O_r).*S{i}{10}(:,a17O_r);
    S{i}{10}(:,a267_s) = 2.*S{i}{10}(:,a12C_s).*S{i}{10}(:,a16O_s).*S{i}{10}(:,a17O_s);
    S{i}{10}(:,a366_r) = S{i}{10}(:,a13C_r).*S{i}{10}(:,a16O_r).*S{i}{10}(:,a16O_r);
    S{i}{10}(:,a366_s) = S{i}{10}(:,a13C_s).*S{i}{10}(:,a16O_s).*S{i}{10}(:,a16O_s);
    S{i}{10}(:,a268_r) = 2.*S{i}{10}(:,a12C_r).*S{i}{10}(:,a16O_r).*S{i}{10}(:,a18O_r);
    S{i}{10}(:,a268_s) = 2.*S{i}{10}(:,a12C_s).*S{i}{10}(:,a16O_s).*S{i}{10}(:,a18O_s);
    S{i}{10}(:,a277_r) = S{i}{10}(:,a12C_r).*S{i}{10}(:,a17O_r).*S{i}{10}(:,a17O_r);
    S{i}{10}(:,a277_s) = S{i}{10}(:,a12C_s).*S{i}{10}(:,a17O_s).*S{i}{10}(:,a17O_s);
    S{i}{10}(:,a376_r) = 2.*S{i}{10}(:,a13C_r).*S{i}{10}(:,a17O_r).*S{i}{10}(:,a16O_r);
    S{i}{10}(:,a376_s) = 2.*S{i}{10}(:,a13C_s).*S{i}{10}(:,a17O_s).*S{i}{10}(:,a16O_s);
    S{i}{10}(:,a278_r) = 2.*S{i}{10}(:,a12C_r).*S{i}{10}(:,a17O_r).*S{i}{10}(:,a18O_r);
    S{i}{10}(:,a278_s) = 2.*S{i}{10}(:,a12C_s).*S{i}{10}(:,a17O_s).*S{i}{10}(:,a18O_s);
    S{i}{10}(:,a368_r) = 2.*S{i}{10}(:,a13C_r).*S{i}{10}(:,a16O_r).*S{i}{10}(:,a18O_r);
    S{i}{10}(:,a368_s) = 2.*S{i}{10}(:,a13C_s).*S{i}{10}(:,a16O_s).*S{i}{10}(:,a18O_s);
    S{i}{10}(:,a377_r) = S{i}{10}(:,a13C_r).*S{i}{10}(:,a17O_r).*S{i}{10}(:,a17O_r);
    S{i}{10}(:,a377_s) = S{i}{10}(:,a13C_s).*S{i}{10}(:,a17O_s).*S{i}{10}(:,a17O_s);
    S{i}{10}(:,a288_r) = S{i}{10}(:,a12C_r).*S{i}{10}(:,a18O_r).*S{i}{10}(:,a18O_r);
    S{i}{10}(:,a288_s) = S{i}{10}(:,a12C_s).*S{i}{10}(:,a18O_s).*S{i}{10}(:,a18O_s);
    S{i}{10}(:,a378_r) = 2.*S{i}{10}(:,a13C_r).*S{i}{10}(:,a17O_r).*S{i}{10}(:,a18O_r);
    S{i}{10}(:,a378_s) = 2.*S{i}{10}(:,a13C_s).*S{i}{10}(:,a17O_s).*S{i}{10}(:,a18O_s);
    S{i}{10}(:,a388_r) = S{i}{10}(:,a13C_r).*S{i}{10}(:,a18O_r).*S{i}{10}(:,a18O_r);
    S{i}{10}(:,a388_s) = S{i}{10}(:,a13C_s).*S{i}{10}(:,a18O_s).*S{i}{10}(:,a18O_s);
    
    % sum abundances of CO2 isotopologues w/ same cardinal mass in sample and ref gas
    S{i}{10}(:,a44_r) = S{i}{10}(:,a266_r);
    S{i}{10}(:,a44_s) = S{i}{10}(:,a266_s);
    S{i}{10}(:,a45_r) = S{i}{10}(:,a267_r)+S{i}{10}(:,a366_r);
    S{i}{10}(:,a45_s) = S{i}{10}(:,a267_s)+S{i}{10}(:,a366_s);
    S{i}{10}(:,a46_r) = S{i}{10}(:,a268_r)+S{i}{10}(:,a277_r)+S{i}{10}(:,a376_r);
    S{i}{10}(:,a46_s) = S{i}{10}(:,a268_s)+S{i}{10}(:,a277_s)+S{i}{10}(:,a376_s);
    S{i}{10}(:,a47_r) = S{i}{10}(:,a278_r)+S{i}{10}(:,a368_r)+S{i}{10}(:,a377_r);
    S{i}{10}(:,a47_s) = S{i}{10}(:,a278_s)+S{i}{10}(:,a368_s)+S{i}{10}(:,a377_s);
    S{i}{10}(:,a48_r) = S{i}{10}(:,a288_r)+S{i}{10}(:,a378_r);
    S{i}{10}(:,a48_s) = S{i}{10}(:,a288_s)+S{i}{10}(:,a378_s);
    S{i}{10}(:,a49_r) = S{i}{10}(:,a388_r);
    S{i}{10}(:,a49_s) = S{i}{10}(:,a388_s);
    
    % calculate the stochastic distribution in sample and ref gas
    S{i}{10}(:,R45_r) = S{i}{10}(:,a45_r)./S{i}{10}(:,a44_r);
    S{i}{10}(:,R45_s) = S{i}{10}(:,a45_s)./S{i}{10}(:,a44_s);
    S{i}{10}(:,R46_r) = S{i}{10}(:,a46_r)./S{i}{10}(:,a44_r);
    S{i}{10}(:,R46_s) = S{i}{10}(:,a46_s)./S{i}{10}(:,a44_s);
    S{i}{10}(:,R47_r) = S{i}{10}(:,a47_r)./S{i}{10}(:,a44_r);
    S{i}{10}(:,R47_s) = S{i}{10}(:,a47_s)./S{i}{10}(:,a44_s);
    S{i}{10}(:,R48_r) = S{i}{10}(:,a48_r)./S{i}{10}(:,a44_r);
    S{i}{10}(:,R48_s) = S{i}{10}(:,a48_s)./S{i}{10}(:,a44_s);
    S{i}{10}(:,R49_r) = S{i}{10}(:,a49_r)./S{i}{10}(:,a44_r);
    S{i}{10}(:,R49_s) = S{i}{10}(:,a49_s)./S{i}{10}(:,a44_s);
    
    % calculate R45, R46, R47, R48, R49 of sample relative to stochastic ref gas
    S{i}{10}(:,R45_star) = ((S{i}{10}(:,d45)./1000)+1).*S{i}{10}(:,R45_r);
    S{i}{10}(:,R46_star) = ((S{i}{10}(:,d46)./1000)+1).*S{i}{10}(:,R46_r);
    S{i}{10}(:,R47_star) = ((S{i}{10}(:,d47)./1000)+1).*S{i}{10}(:,R47_r);
    S{i}{10}(:,R48_star) = ((S{i}{10}(:,d48)./1000)+1).*S{i}{10}(:,R48_r);
    S{i}{10}(:,R49_star) = ((S{i}{10}(:,d49)./1000)+1).*S{i}{10}(:,R49_r);
    
    % calculate the pieces of D47, D48, D49
    S{i}{10}(:,D_45) = ((S{i}{10}(:,R45_star)./S{i}{10}(:,R45_s))-1).*1000;
    S{i}{10}(:,D_46) = ((S{i}{10}(:,R46_star)./S{i}{10}(:,R46_s))-1).*1000;
    S{i}{10}(:,D_47) = ((S{i}{10}(:,R47_star)./S{i}{10}(:,R47_s))-1).*1000;
    S{i}{10}(:,D_48) = ((S{i}{10}(:,R48_star)./S{i}{10}(:,R48_s))-1).*1000;
    S{i}{10}(:,D_49) = ((S{i}{10}(:,R49_star)./S{i}{10}(:,R49_s))-1).*1000;
    
    % calculate the D47, D48, D49
    S{i}{10}(:,D47_raw) = S{i}{10}(:,D_47) - S{i}{10}(:,D_46) - S{i}{10}(:,D_45);
    S{i}{10}(:,D48_raw) = S{i}{10}(:,D_48) - 2.*S{i}{10}(:,D_46);
    S{i}{10}(:,D49_raw) = S{i}{10}(:,D_49) - 2.*S{i}{10}(:,D_46) - S{i}{10}(:,D_45);
    
    % look up HG line and secondary transfer function based on the date
    sd = S{i}{10}(1,day);
    sm = S{i}{10}(1,month);
    sy = S{i}{10}(1,year);
    m  = S{i}{11}(mach_id);
    for j = 1:size(fd,1)
        if (fd(j)<=sd&&sd<=ld(j)&&fm(j)<=sm&&sm<=lm(j)&&fy(j)<=sy&&sy<=ly(j)&&strcmp(m,mach(j)))
            S{i}{10}(:,HG_slope) = slope_HG(j);
            S{i}{10}(:,HG_int)   = intercept_HG(j);
            %S{i}{10}(:,resid)    = residual(j);
            S{i}{10}(:,STF_slope) = slope_STF(j);
            S{i}{10}(:,STF_int)   = intercept_STF(j);
        end
    end
    
    % correct D47 for HG line, and stretching
    S{i}{10}(:,D47_HG) = S{i}{10}(:,D47_raw)-(S{i}{10}(:,d47).*S{i}{10}(:,HG_slope)+S{i}{10}(:,HG_int));
    S{i}{10}(:,D47_stretch) = S{i}{10}(:,D47_HG).*-0.8453./S{i}{10}(:,HG_int);
    %S{i}{10}(:,D47_acid) = S{i}{10}(:,D47_stretch) + 0.08;
    %S{i}{10}(:,D47_aq) = S{i}{10}(:,D47_acid) - S{i}{10}(:,resid);
    
    % convert D47 to absolute reference frame, then correct for acid
    S{i}{10}(:,D47_RF) = S{i}{10}(:,D47_stretch).*S{i}{10}(:,STF_slope)+S{i}{10}(:,STF_int);
    S{i}{10}(:,D47_RF_AC) = S{i}{10}(:,D47_RF) + 0.092;
    
    % calculate temperature using equation 9 from Dennis et al. (2011)
    %S{i}{10}(:,temp_aq) = (59200.*(S{i}{10}(:,D47_aq)+0.02).^-1).^0.5-273.15;
    S{i}{10}(:,temp_aq) = (63600./(S{i}{10}(:,D47_RF_AC)+0.0047)).^0.5-273.15;
    
    % calculate the d18O_mineral_SMOW from the d18O_gas_SMOW
    S{i}{10}(:,d18O_min_PDB_aq) = ((((S{i}{10}(:,d18O_SMOW)-30.86)./1.03086)+1000)./1.00821)-1000;
    S{i}{10}(:,d18O_min_SMOW_aq) = S{i}{10}(:,d18O_min_PDB_aq).*1.03086+30.86;
    
    % calculate d18O_water
    S{i}{10}(:,d18O_water_aq) = ((1./exp((18030./(S{i}{10}(:,temp_aq)...
        +273.15)-32.42)./1000)).*(S{i}{10}(:,d18O_min_SMOW_aq)+1000))-1000;
    
    % grab only my samples (excluding standards & HGs & equil gases)
    if strcmp(S{i}{11}(type),'sample') && isempty(regexp(S{i}{11}{ID},'G\d* ','ONCE'))==0
        for j = 10:12
            Samples{sample_counter}{j} = S{i}{j};
        end
        sample_counter = sample_counter + 1;
    end
    
end

% if a sample has 9 acquisitions, delete the first acquisition
for i = 1:num_of_samps
    if size(Samples{i}{10},1)==9
        Samples{i}{10} = Samples{i}{10}(2:end,:);
    end
end

% find samples with the same IDs and concatenate their acquisitions
counter = 0;
while_stop = num_of_samps;
i = 1;
while i <= while_stop
    
    match = Samples{i}{11}{ID};
    j = 1;
    
    while j <= while_stop
        if strcmp(match,Samples{j}{11}{ID}) && i~=j
            
            % concatenate replicate samples
            Samples{i}{10} = [Samples{i}{10}; Samples{j}{10}];
            Samples{i}{11} = [Samples{i}{11}; Samples{j}{11}];
            
            % delete the other replicates
            index = true(1,size(Samples,1));
            index(j) = false;
            Samples = Samples(index);
            while_stop = while_stop - 1;
            
            counter = counter + 1;
        end
        % advance the j while loop
        j = j+1;
    end
    % calculate how many replicates and acquisitions there are
    Samples{i}{12}(reps) = size(Samples{i}{11},1);
    Samples{i}{12}(aqs) = size(Samples{i}{10},1);
    
    % advance the i while loop
    i = i+1;
end

% load trace metal data 
Metals_data = load('trace_metals_data.mat','Metals');
Metals = Metals_data.Metals;

% fill out the flat list array
for i = 1:size(Samples,1)
    
    % calculate the mean d13C, d18O, d47, D47, d48, D48, and sigmas
    Samples{i}{12}(d13C_PDB)   = mean(Samples{i}{10}(:,d13C));
    Samples{i}{12}(sigma_d13C) =  std(Samples{i}{10}(:,d13C));
    Samples{i}{12}(sterr_d13C) =  std(Samples{i}{10}(:,d13C))/sqrt(Samples{i}{12}(aqs));
    Samples{i}{12}(d18O_gas_SMOW)  = mean(Samples{i}{10}(:,d18O_SMOW));
    Samples{i}{12}(d18O_min_PDB)   = mean(Samples{i}{10}(:,d18O_min_PDB_aq));
    Samples{i}{12}(sigma_d18O_min) =  std(Samples{i}{10}(:,d18O_min_PDB_aq));
    Samples{i}{12}(sterr_d18O_min) =  std(Samples{i}{10}(:,d18O_min_PDB_aq))/sqrt(Samples{i}{12}(aqs));
    Samples{i}{12}(d18O_min_SMOW)       = mean(Samples{i}{10}(:,d18O_min_SMOW_aq));
    Samples{i}{12}(sigma_d18O_min_SMOW) =  std(Samples{i}{10}(:,d18O_min_SMOW_aq));
    Samples{i}{12}(d47_flat)  = mean(Samples{i}{10}(:,d47));
    Samples{i}{12}(sigma_d47) =  std(Samples{i}{10}(:,d47));
    Samples{i}{12}(sterr_D47) =  std(Samples{i}{10}(:,d47))/sqrt(Samples{i}{12}(aqs));
    Samples{i}{12}(D47)       = mean(Samples{i}{10}(:,D47_RF_AC));
    Samples{i}{12}(sigma_D47) =  std(Samples{i}{10}(:,D47_RF_AC));
    Samples{i}{12}(d48_flat)  = mean(Samples{i}{10}(:,d48));
    Samples{i}{12}(D48_flat)  = mean(Samples{i}{10}(:,D48_raw));
    
    % calculate the temp and d18O_water and sigmas from the mean values
    Samples{i}{12}(temp) = (63600*(Samples{i}{12}(D47)+0.0047)^-1)^0.5-273.15;
    Samples{i}{12}(sigma_temp) = ((63600/(Samples{i}{12}(D47)+0.0047))^0.5)...
        *0.5/(Samples{i}{12}(D47)+0.0047)*Samples{i}{12}(sigma_D47);
    Samples{i}{12}(sterr_temp) = Samples{i}{12}(sigma_temp)/sqrt(Samples{i}{12}(aqs));
    Samples{i}{12}(d18O_water) = ((1/exp((18030/(Samples{i}{12}(temp)...
        +273.15)-32.42)/1000))*(Samples{i}{12}(d18O_min_SMOW)+1000))-1000;
    Samples{i}{12}(sigma_d18O_water) = ((1/exp((18030/(Samples{i}{12}(temp)...
        +273.15)-32.42)/1000))*(Samples{i}{12}(d18O_min_SMOW)+1000))...
        *((18.03*Samples{i}{12}(sigma_temp)/(Samples{i}{12}(temp)+273.15)...
        ^2)^2+(Samples{i}{12}(sigma_d18O_min_SMOW)/(Samples{i}{12}(d18O_min_SMOW)+1000))^2)^0.5;
    Samples{i}{12}(sterr_d18Ow) = Samples{i}{12}(sigma_d18O_water)/sqrt(Samples{i}{12}(aqs));

    % translate the sample name into section height & sample type
    [match_loc nomatch_loc] = regexp(Samples{i}{11}{1,ID},'G\d* ','match','split');
    location_string = regexp(match_loc,'\d*','match');
    location = str2double(location_string{1});
    if regexp(nomatch_loc{2},'LVF')
        loc_height = -2;
        [match_h nomatch_h] = regexp(nomatch_loc{2},'LVF','match','split');
        type_string = nomatch_h{2};  
    elseif regexp(nomatch_loc{2},'UF')
        loc_height = 32;
        [match_h nomatch_h] = regexp(nomatch_loc{2},'UF','match','split');
        type_string = nomatch_h{2};
    elseif isempty(regexp(nomatch_loc{2},' to ','ONCE'))==0
        [match_h nomatch_h]= regexp(nomatch_loc{2},'-?\d*\.\d*','match','split');
        loc_height_base = str2double(match_h{1});
        loc_height_top  = str2double(match_h{2});
        loc_height = mean([loc_height_base loc_height_top]);
        type_string = nomatch_h{3};
    else
        [match_h nomatch_h]= regexp(nomatch_loc{2},'-?\d*\.\d*','match','split');
        if isempty(match_h)
            no_height = 1;
            type_string = nomatch_h{1};
        else 
            no_height = 0;
            loc_height = str2double(match_h{1});
            type_string = nomatch_h{2};
        end
    end
    
    if location==1 || location==2 || location==3
        Samples{i}{12}(height) = loc_height;
    elseif location==4
        Samples{i}{12}(height) = loc_height + 6;
    elseif location==5
        Samples{i}{12}(height) = loc_height + 46.1;
    elseif location==6
        if no_height==1
        Samples{i}{12}(height) = 13;
        else 
        Samples{i}{12}(height) = loc_height + 11.4;
        end
    elseif location==7
        Samples{i}{12}(height) = loc_height + 17;
    elseif location==8
        Samples{i}{12}(height) = loc_height + 100;
    elseif location==9
        Samples{i}{12}(height) = loc_height + 313;
    elseif location==10
        Samples{i}{12}(height) = loc_height + 319;
    elseif location==11
        Samples{i}{12}(height) = loc_height + 1;
    elseif location==12
        if no_height==1
        Samples{i}{12}(height) = -5.5;
        else
        Samples{i}{12}(height) = loc_height - 9.5;
        end
    elseif location==13
        Samples{i}{12}(height) = loc_height + 75;
    end
    
    if isempty(regexp(type_string,'cement','ONCE'))==0
        Samples{i}{12}(sample_type) = 7;
    elseif isempty(regexp(type_string,'micrite','ONCE'))==0
        Samples{i}{12}(sample_type) = 5;
    elseif isempty(regexp(type_string,'M\d*','ONCE'))==0
        Samples{i}{12}(sample_type) = 5;
    elseif isempty(regexp(type_string,'B\d*','ONCE'))==0
        Samples{i}{12}(sample_type) = 19;
    elseif isempty(regexp(type_string,'R\d*','ONCE'))==0
        Samples{i}{12}(sample_type) = 23;
    elseif isempty(regexp(type_string,'bryozoan','ONCE'))==0
        Samples{i}{12}(sample_type) = 18;
    elseif isempty(regexp(type_string,'crinoid','ONCE'))==0
        Samples{i}{12}(sample_type) = 2;
    elseif isempty(regexp(type_string,'gastropod','ONCE'))==0
        Samples{i}{12}(sample_type) = 15;
    elseif isempty(regexp(type_string,'ostracod','ONCE'))==0
        Samples{i}{12}(sample_type) = 16;
    end
    
    if isempty(regexp(type_string,'leftover','ONCE'))==0
        Samples{i}{12}(bad) = 1;
    elseif isempty(regexp(type_string,'bad','ONCE'))==0
        Samples{i}{12}(bad) = 1;
    elseif isempty(regexp(type_string,'top','ONCE'))==0
        Samples{i}{12}(bad) = 1;
    elseif isempty(regexp(type_string,'micrite','ONCE'))==0
        Samples{i}{12}(bad) = 1;
    elseif isempty(regexp(type_string,'M\d*','ONCE'))==0
        Samples{i}{12}(bad) = 1;
    elseif isempty(regexp(type_string,'cement','ONCE'))==0
        Samples{i}{12}(bad) = 1;
    else
        Samples{i}{12}(bad) = 0;
    end

    % look up bulk trace metal concentrations for each sample
    for j = 1:size(Metals,1)
        if strcmp(Samples{i}{11}{1,ID},Metals{j}{5})
            Samples{i}{12}(Fe)       = Metals{j}{6}(1);
            Samples{i}{12}(Fe_sterr) = Metals{j}{6}(2);
            Samples{i}{12}(Mg)       = Metals{j}{6}(4);
            Samples{i}{12}(Mg_sterr) = Metals{j}{6}(5);
            Samples{i}{12}(Mn)       = Metals{j}{6}(7);
            Samples{i}{12}(Mn_sterr) = Metals{j}{6}(8);
            Samples{i}{12}(Sr)       = Metals{j}{6}(10);
            Samples{i}{12}(Sr_sterr) = Metals{j}{6}(11);
        end
    end

end

% preallocate export arrays
section_csv      = zeros(size(Samples,1),10);
comparisons_txt  =  cell(8,9);
trace_metals_txt =  cell(size(Samples,1),17);
Brand_csv        = zeros(size(Samples,1),6);
summary_txt      =  cell(size(Samples,1),1);
summary_csv      = zeros(size(Samples,1),size(Samples{1}{12},2));
comparisons_row = 1;
metals_row = 1;

% loop over Samples to fill export arrays row-by-row
for i = 1:size(Samples,1)
    
    % fill export array for section figure
    section_csv(i,1)  = Samples{i}{12}(height);
    section_csv(i,2)  = Samples{i}{12}(d13C_PDB);
    section_csv(i,3)  = Samples{i}{12}(sterr_d13C);
    section_csv(i,4)  = Samples{i}{12}(d18O_min_PDB);
    section_csv(i,5)  = Samples{i}{12}(sterr_d18O_min);
    section_csv(i,6)  = Samples{i}{12}(temp);
    section_csv(i,7)  = Samples{i}{12}(sterr_temp);
    section_csv(i,8)  = Samples{i}{12}(sample_type);
    section_csv(i,9)  = Samples{i}{12}(bad);
    section_csv(i,10) = Samples{i}{12}(reps);
    
    if i==10||i==16||i==19||i==22||i==37||i==40||i==43||i==47
    % fill export array for comparisons figure
    comparisons_txt{comparisons_row,1} = Samples{i}{11}{1,ID};
    comparisons_txt{comparisons_row,2} = Samples{i}{12}(sample_type);
    comparisons_txt{comparisons_row,3} = Samples{i}{12}(d13C_PDB);
    comparisons_txt{comparisons_row,4} = Samples{i}{12}(sterr_d13C);
    comparisons_txt{comparisons_row,5} = Samples{i}{12}(d18O_min_PDB);
    comparisons_txt{comparisons_row,6} = Samples{i}{12}(sterr_d18O_min);
    comparisons_txt{comparisons_row,7} = Samples{i}{12}(temp);
    comparisons_txt{comparisons_row,8} = Samples{i}{12}(sterr_temp);
    comparisons_txt{comparisons_row,9} = Samples{i}{12}(reps);
    if Samples{i}{12}(sample_type) == 19
        comparisons_txt{comparisons_row,10} = 'brachiopod';
    elseif Samples{i}{12}(sample_type) == 5
        comparisons_txt{comparisons_row,10} = 'micrite';
    end
    comparisons_row = comparisons_row + 1;
    end
    
    if Samples{i}{12}(Fe) ~= 0
    % fill export array for trace metals figures
    trace_metals_txt{metals_row,1} = Samples{i}{12}(d13C_PDB);
    trace_metals_txt{metals_row,2} = Samples{i}{12}(sterr_d13C);
    trace_metals_txt{metals_row,3} = Samples{i}{12}(d18O_min_PDB);
    trace_metals_txt{metals_row,4} = Samples{i}{12}(sterr_d18O_min);
    trace_metals_txt{metals_row,5} = Samples{i}{12}(temp);
    trace_metals_txt{metals_row,6} = Samples{i}{12}(sterr_temp);
    trace_metals_txt{metals_row,7} = Samples{i}{12}(sample_type);
    trace_metals_txt{metals_row,8} = Samples{i}{12}(bad);
    trace_metals_txt{metals_row,9} = Samples{i}{12}(Fe);
    trace_metals_txt{metals_row,10} = Samples{i}{12}(Fe_sterr);
    trace_metals_txt{metals_row,11} = Samples{i}{12}(Mg);
    trace_metals_txt{metals_row,12} = Samples{i}{12}(Mg_sterr);
    trace_metals_txt{metals_row,13} = Samples{i}{12}(Mn);
    trace_metals_txt{metals_row,14} = Samples{i}{12}(Mn_sterr);
    trace_metals_txt{metals_row,15} = Samples{i}{12}(Sr);
    trace_metals_txt{metals_row,16} = Samples{i}{12}(Sr_sterr);
    trace_metals_txt{metals_row,17} = Samples{i}{11}{1,ID};
    metals_row = metals_row + 1;
    end
    
    % fill export array for Brand figure
    Brand_csv(i,1) = Samples{i}{12}(sample_type);
    Brand_csv(i,2) = Samples{i}{12}(bad);
    Brand_csv(i,3) = Samples{i}{12}(Mn);
    Brand_csv(i,4) = Samples{i}{12}(Mn_sterr);
    Brand_csv(i,5) = Samples{i}{12}(Sr);
    Brand_csv(i,6) = Samples{i}{12}(Sr_sterr);
    
    % fill export summary array
    summary_txt{i}   = Samples{i}{11}{1,ID};
    summary_csv(i,:) = Samples{i}{12};
    
end

% export csv data to make figures in igor
csvwrite('fig_data_section.csv',section_csv);
%trace_metals_csv = trace_metals_csv(any(trace_metals_csv,2),:);
%csvwrite('fig_data_trace_metals.csv',trace_metals_csv);
csvwrite('fig_data_Brand.csv',Brand_csv);

% export cell array data to make comparisons figure in igor
[nrows,~] = size(comparisons_txt);
names_array = cell(nrows,1);
for i = 1:nrows
    names_array{i,1} = comparisons_txt{i,1};
end
[sorted, index] = sort(names_array);
comparisons_txt = comparisons_txt(index,:);
fid = fopen('fig_data_comparisons', 'w');
label = 1;
for row = 1:nrows
    if mod(row,2)==1 || row==8
        fprintf(fid,'%s,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d\n',comparisons_txt{row,:},label);
        label = label + 1;
    elseif mod(row,2)==0
        fprintf(fid,'%s,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d\n',comparisons_txt{row,:},label);
        fprintf(fid,',,,,,,,,,,%d\n',label+1);
        label = label + 2;
    end
end
fclose(fid);

% export cell array data to make trace metals figure in igor
trace_metals_export = reshape(trace_metals_txt(~cellfun('isempty',trace_metals_txt)),[],17);
[nrows,~] = size(trace_metals_export);
fid = fopen('fig_data_trace_metals', 'w');
for row = 1:nrows
    fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s\n',trace_metals_export{row,:});
end
fclose(fid);

% save Samples variable so I can use it to make the ellipses figure
save('fig_data_ellipses','Samples');

% export summary files
[nrows,~] = size(summary_txt);
fid = fopen('summary_names', 'w');
for row = 1:nrows
    fprintf(fid,'%s\n',summary_txt{row});
end
fclose(fid);
csvwrite('summary_data.csv',summary_csv);



