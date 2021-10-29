%%%% Plot RAS nutrient data for quick look plots,
%%%% including all data, with their associated flags
%%% and CTD data (CWE addition)
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
NOx_lim = 20;
PO4_Si_lim = 6;
% set the lab temperature for density calculations
lab_temp = 20;

% plotting the RAS data first
fig = figure();
cl = get(groot,'defaultAxesColorOrder');

for k = 1:length(files)
    
    fn = files(k);
    file = [fn.folder '/' fn.name]
    
    plotVar={'NTRI','PHOS' ,'SLCA'};
    
    lat_m = ncread(file, 'LATITUDE');
    lon_m = ncread(file,'LONGITUDE');
    
    try
        for i = 1:3
        var = ncread(file, plotVar{i});
        
        var_unit = ncreadatt(file, plotVar{i}, 'units');
        var_name = ncreadatt(file, plotVar{i}, 'long_name')

        varQCname = strsplit(ncreadatt(file, plotVar{i}, 'ancillary_variables'), ' ');
        varQC = ncread(file, varQCname{1});

        time = ncread(file, 'TIME') + datetime(1950,1,1);
        deployment_code = ncreadatt(file, '/', 'deployment_code');
        nominal_depth = ncread(file, 'NOMINAL_DEPTH');

        timestart = datetime(ncreadatt(file, '/', 'time_deployment_start'), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
        timeend = datetime(ncreadatt(file, '/', 'time_deployment_end'), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');

        % plotting limits set to 20 for NOx and 6 for PO4 and Si
        % therefore all values higher than that will be set to the
        % maximum in order to have all axes in all plots going forward
        % the same.
               
        % plot the data with line through good data
        if i==1
            yyaxis left
            plottingVar = var;
            plottingmsk = var > NOx_lim; 
            plottingVar(plottingmsk,1) = NOx_lim;
        elseif i>1
            yyaxis right
            plottingVar = var;
            plottingmsk = var>PO4_Si_lim; 
            plottingVar(plottingmsk,1) = PO4_Si_lim;
            
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
t = title({deployment_code ;'\rm \color[rgb]{0 0.4470 0.7410}NO_{x}   \color[rgb]{0.8500 0.3250 0.0980}Phosphate   \color[rgb]{0.9290 0.6940 0.1250}Silicate'; '\rm \color[rgb]{black}QC flag 1 = line, QC flag 2 = circle, QC flag 3 = diamond, QC flag 4 = square'},'Interpreter','tex');
%xlim([timestart-10 timeend+10]);
yyaxis left
ylabel('NO_{x} \mumol kg^{-1}', 'Interpreter', 'tex')
yyaxis right
ylabel('Phosphate \color[rgb]{black}and \color[rgb]{0.9290 0.6940 0.1250}silicate \color[rgb]{black}\mumol kg^{-1}','Interpreter', 'tex')
%print( fig, '-dpng', [files(1).folder '/' deployment_code '-RAS-nuts-report_plot.png'])
figurename = [files(1).folder '/' deployment_code '-RAS-nuts-report_plot.png'];
end


% next add the CTD data

% first get all the CTD files together, hyrdo and sensor data
for i = 1:size(voyageID,1)
    search_path = ['C:\Users\cawynn\cloudstor\Shared\CTD\' voyageID{i}];
    cd(search_path)
    f = dir('*Hydro*.nc');
    CTDhydro{1,i} = f;
    s = dir('*CtdAvg*.nc');
    CTDsensor{1,i} = s;
    
end


cl = get(groot,'defaultAxesColorOrder');
for l = 1:length(CTDhydro)
    CTDhydrofiles = CTDhydro{1,l};
    CTDsensorfiles = CTDsensor{1,l};


    for k = 1:length(CTDhydrofiles)

        fn = CTDhydrofiles(k);
        file = [fn.folder '/' fn.name]

         plotVar={'nox','phosphate' ,'silicate'};
         press = ncread(file, 'pressure');
         psal = ncread(file,'salinity')';
         temp = ncread(file,'temperature')';
         lat = ncread(file,'latitude');
         lon = ncread(file,'longitude');
         deployment = ncreadatt(file, '/', 'Deployment');
         
         % check whether the CTD was close enough in space to the mooring
         if (lat-lat_m) <= lat_t && (lon-lon_m) <= lon_t
             
             % find the matching Sensor file by deployment number
             filesensor=[];
             for m = 1:length(CTDsensorfiles)
                 fnsensor = CTDsensorfiles(m);
                 sensor = [fnsensor.folder '/' fnsensor.name];
                 deploymentSensor = ncreadatt(sensor, '/', 'Deployment');
                 if str2double(deploymentSensor) == deployment
                     filesensor = sensor;
                 end
             end
             if size(filesensor,1) == 0
                 disp('did not find matching sensor file');
             end

             pressSensor = ncread(filesensor, 'pressure');
             psalSensor = ncread(filesensor,'salinity')';
             tempSensor = ncread(filesensor,'temperature')';
             latSensor = ncread(filesensor,'latitude');
             lonSensor = ncread(filesensor,'longitude');
                                         
             try
                for i = 1:3
                    var = ncread(file, plotVar{i})';

                    var_unit = ncreadatt(file, plotVar{i}, 'units');
                    var_name = ncreadatt(file, plotVar{i}, 'long_name')

                    varQCname = [plotVar{i} 'Flag'];
                    varQC = ncread(file, varQCname)';

                    time = datetime(ncread(file, 'woce_date'), 'ConvertFrom' , 'yyyymmdd');% + datetime(1950,1,1);

                    % only use the bottle closest to 10m
                    [row,col] = find(~isnan(var));
                    var_N = var(row);
                    press_N = press(row);
                    [idx,d] = knnsearch(press_N, comp_depth);
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

                            plot(time, var_conv,'*-','MarkerSize',5,'LineWidth',0.8,'DisplayName',plotVar{i}, 'Color', cl(i, :));
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
ylim([0 , 20]);
yyaxis right
ylim([0 , 6]);
end
figurename = [deployment_code '-RAS-nuts-report_plot.png'];
figure_path = 'C:\Users\cawynn\OneDrive - University of Tasmania\Documents\GitHub\RAS_data_plotting\RAS_data_plotting\figures'
cd(figure_path)
print( fig, '-dpng', figurename);
