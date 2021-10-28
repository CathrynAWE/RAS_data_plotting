%%%% Plot RAS ALK and salinity data for quick look plots,
%%%% including all data, with their associated flags
%%% and CTD data (CWE)
clear;
close all;

% add gsw functions for conversion of umol/l to umol/kg (CTD only)
folder = fileparts('C:\Users\cawynn\OneDrive - University of Tasmania\Documents\MATLAB\gsw_matlab_v3_06_13'); 

addpath(genpath(folder));

% set the voyage IDs to compare with 
%voyageID = {'in2020_V09'};% 'in2021_V02'};
voyageID = {'in2019_V02'; 'in2020_V09'};
% RAS data path
search_path = 'C:\Users\cawynn\cloudstor\Shared\RAS share';
cd(search_path)
%files = dir('IMOS_DWM-SOTS_KPSRT_20200906_SOFS_FV01_SOFS-9-2020-RAS-3-48-500-5.3m_END-20200924_C-20211001.nc');
files = dir('IMOS_DWM-SOTS_KPSRT_20190331_SOFS_FV01_SOFS-8-2019-RAS3-48-500FH-5.3m_END-20200410_C-20210922.nc');


% set the CTD bottle depth to compare with RAS data
comp_depth = 10;
% set the tolerances for CTD bottle depth
comp_depth_t = 10;
% set the tolerances for lat and lon of the cast as compared to the mooring
lat_t = 1;
lon_t = 1;
% set the limits for plotting axes
% plotting limits set to 2260 to 2360 umol/kg for Talk (based on QC
% report) and 34 to 35.5 for PSAL.
% Therefore all values higher than that will be set to the
% maximum in order to have all axes in all plots going forward
% the same.
Talk_min = 2260;
Talk_max = 2360;
Sal_min = 34;
Sal_max = 35.5;

% set the lab temperature for density calculations
lab_temp = 20;

% plotting the RAS data first
fig = figure();
cl = get(groot,'defaultAxesColorOrder');

for k = 1:length(files)
    
    fn = files(k);
    file = [fn.folder '/' fn.name];
    
    plotVar={'PSAL','TALK'};
    lat_m = ncread(file, 'LATITUDE');
    lon_m = ncread(file,'LONGITUDE');
    
    try
        for i = 1:2
        var = ncread(file, plotVar{i});
        
        var_unit = ncreadatt(file, plotVar{i}, 'units');
        var_name = ncreadatt(file, plotVar{i}, 'long_name');

        varQCname = strsplit(ncreadatt(file, plotVar{i}, 'ancillary_variables'), ' ');
        varQC = ncread(file, varQCname{1});

        time = ncread(file, 'TIME') + datetime(1950,1,1);
        deployment_code = ncreadatt(file, '/', 'deployment_code');
        nominal_depth = ncread(file, 'NOMINAL_DEPTH');

        timestart = datetime(ncreadatt(file, '/', 'time_deployment_start'), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
        timeend = datetime(ncreadatt(file, '/', 'time_deployment_end'), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
                
        % plot the data with line through good data
        if i==1
            yyaxis left % salinity
            plottingVar = var;
            plottingmsk_max = var >= Sal_max; 
            plottingVar(plottingmsk_max,1) = Sal_max;
            plottingmsk_min = var <= Sal_min; 
            plottingVar(plottingmsk_min,1) = Sal_min;
            
        elseif i>1
            yyaxis right % alkalinity
            plottingVar = var;
            plottingmsk_max = var >= Talk_max; 
            plottingVar(plottingmsk_max,1) = Talk_max;
            plottingmsk_min = var <= Talk_min; 
            plottingVar(plottingmsk_min,1) = Talk_min;
         end
        
        p1 = plot(time(varQC <= 2), plottingVar(varQC <= 2),'.-','MarkerSize',5,'LineWidth',1.2,'DisplayName',plotVar{i}, 'Color', cl(i,:));
       
        hold on
        
        % plot the QC info
             
        plot(time(varQC == 2), plottingVar(varQC == 2), 'o', 'MarkerSize',6, 'MarkerFaceColor', cl(i,:), 'Color', cl(i, :));

        plot(time(varQC == 3), plottingVar(varQC == 3),'d','MarkerSize',6,'MarkerFaceColor', cl(i,:), 'Color', cl(i, :));

        plot(time(varQC == 4), plottingVar(varQC == 4),'s','MarkerSize', 6, 'MarkerFaceColor', cl(i,:), 'Color', cl(i, :))
        
          %i = i + 1;
           
        end
        
      catch
    end
    %grid on
t = title({deployment_code ;'Salinity and Alkalinity'; 'QC flag 1 = line, QC flag 2 = circle, QC flag 3 = diamond, QC flag 4 = square'},'Interpreter','tex');
xlim([timestart-10 timeend+10]);
yyaxis left
ylim([Sal_min Sal_max])
ylabel('Salinity','Interpreter', 'tex')

yyaxis right
ylim([Talk_min Talk_max])
ylabel('Alk \mumol kg^{-1}', 'Interpreter', 'tex')
%hold off
%print( fig, '-dpng', [files(1).folder '/' deployment_code '-RAS-nuts-report_plot.png'])
figurename = [files(1).folder '/' deployment_code '-RAS-Alk-report_plot.png'];
end

% next add the CTD data
% first get relevant Alk / salinity data
Alk_file = 'C:\Users\cawynn\cloudstor\Shared\CTD\DIC_Alk\DIC_Alk_compilation.xlsx';
Alk_data = readtable(Alk_file);
data=[];
for i = 1:length(voyageID)
    data = [data; Alk_data(strcmp(Alk_data.VoyageID,voyageID{i}),:)]; 
end

% now add salinity and alkalinity to the RAS plot
cl = get(groot,'defaultAxesColorOrder');
  
% take out the CTD casta there not close enough in space to the mooring
data(data.Lat-lat_m >= lat_t,:) = [];
data(data.Long-lon_m >= lon_t,:) = [];

if isempty(data)
    disp('no CTD cast close enough')
else
    % use SOMMA salinity when sensor salinity not available
    data.Psal(isnan(data.Psal)) = data.SOMMAPsal(isnan(data.Psal));
         
    % only use the bottle closest to 10m
    for s = 1:length(voyageID)
        data_pre = data(strcmp(data.VoyageID,voyageID{s}),:);
        [idx,d] = knnsearch(data_pre.Pres, comp_depth);
    % limit the search to +- 10m
    if d <= comp_depth_t
                    ind1 = find(var(:) == var_N(idx));
                    % calculate seawater density to convert nutrients
                    % in umol/l to umol/kg
                    SA = gsw_SA_from_SP(psal,press,lon,lat);
                    % use lab temperature instead of insitu temp
                    lab_t=ones(length(press),1);
                    lab_t=lab_t*lab_temp;
                    CT = gsw_CT_from_t(SA,lab_t,press);
                    sigma0 = gsw_sigma0(SA, CT);
                    % see if this results in NaN because data is missing
                    if ~isnan(sigma0(ind1))
                        % convert the value at 10m for plotting
                        var_conv = ((var(ind1)/(sigma0(ind1) + 1000))*1000);
                        disp('density calculated from bottle data');
                    else
                        % double checking the files' deployments match
                        if deployment == deploymentSensor 
                            [ds, ixs] = min(abs(pressSensor-press(ind1)));
                            % calculate seawater density to convert nutrients
                            % in umol/l to umol/kg
                            SASensor = gsw_SA_from_SP(psalSensor(ixs),pressSensor(ixs),lonSensor,latSensor) ;
                            % use lab temperature instead of insitu temp
                            %CTSensor = gsw_CT_from_t(SA,tempSensor(ixs),pressSensor(ixs));
                            CTSensor = gsw_CT_from_t(SA,lab_temp,pressSensor(ixs));
                            sigma0Sensor = gsw_sigma0(SA, CT);
                            var_conv = ((var(ind1)/(sigma0Sensor + 1000))*1000);
                            dips('density calculated from sensor data');
                        end
                    end


                    % plot the data with line through good data
                    if i==1
                        yyaxis left
                        plottingVar = var_conv;
                        plottingmsk = var_conv > NOx_lim; 
                        plottingVar(plottingmsk,1) = NOx_lim;
                    elseif i>1
                        yyaxis right
                        plottingVar = var_conv;
                        plottingmsk = var_conv > PO4_Si_lim; 
                        plottingVar(plottingmsk,1) = PO4_Si_lim;
                     end
                    if varQC(ind1) ==0

                        plot(time, plottingVar,'*-','MarkerSize',5,'LineWidth',1.2,'DisplayName',plotVar{i}, 'Color', cl(i, :));
                        hold on
                    else
                        disp('nearest bottle nutrient data has bad QC');
                    end


                else
                    disp('nearest bottle more than 10m away')
                end
                %i = i + 1;
            end

         catch
                  %i = i + 1;

         end
     else
         disp('CTD cast too far away')
     end
end
grid on
yyaxis left
ylim([0 , NOx_lim]);
yyaxis right
ylim([0 , PO4_Si_lim]);

end

%print( fig, '-dpng', figurename);

