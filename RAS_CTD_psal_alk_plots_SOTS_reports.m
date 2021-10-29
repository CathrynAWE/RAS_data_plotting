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
search_path = 'C:\Users\cawynn\OneDrive - University of Tasmania\Documents\GitHub\RAS_data_plotting\RAS_data_plotting\netCDFs_for_plotting';
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
%xlim([timestart-10 timeend+10]);
yyaxis left
ylim([Sal_min Sal_max])
ylabel('Salinity','Interpreter', 'tex')

yyaxis right
ylim([Talk_min Talk_max])
ylabel('Alk \mumol kg^{-1}', 'Interpreter', 'tex')
end

% next add the CTD data
% first get relevant Alk / salinity data
Alk_file = 'C:\Users\cawynn\cloudstor\Shared\CTD\DIC_Alk\DIC_Alk_compilation.xlsx';
Alk_data = readtable(Alk_file);
data=[];
for i = 1:length(voyageID)
    data = [data; Alk_data(strcmp(Alk_data.VoyageID,voyageID{i}),:)]; 
end
  
% take out the CTD casta there not close enough in space to the mooring
data(data.Lat-lat_m >= lat_t,:) = [];
data(data.Long-lon_m >= lon_t,:) = [];

if isempty(data)
    disp('no CTD cast close enough')
else
    % use SOMMA salinity when sensor salinity not available
    data.Psal(isnan(data.Psal)) = data.SOMMAPsal(isnan(data.Psal));
    date=datetime(data.Date,'InputFormat','yyyy-mm-dd HH:MM:SS');
    data.Date_2 = date;
    % now split the data frame into a structure, one per cast
    % a complicated way of getting the number of rows per cast
   
    [uv,idx] = unique(date);
    for i = 1:length(idx)-1
        r(i) = idx(i+1)-idx(i);
    end
    r=r';
    r(end+1) = length(date)+1-idx(end);
    % so that I can feed that information into the function mat2cell
    data_n = mat2cell(data, r, size(data,2));
        
% now add salinity and alkalinity to the RAS plot
% cl = get(groot,'defaultAxesColorOrder');
    % only use the bottle closest to 10m
    for s = 1:size(data_n,1)
           [idx,d] = knnsearch(data_n{s}.Pres, comp_depth);
        % limit the search to +- 10m
            if d <= comp_depth_t
                plottingDate = data_n{s}.Date_2(idx);
         % plot the data with line through good data
                yyaxis left % salinity
                PSAL = data_n{s}.Psal(idx);
                PSALmsk_max = PSAL > Sal_max; 
                PSAL(PSALmsk_max,1) = Sal_max;
                PSALmsk_min = PSAL < Sal_min; 
                PSAL(PSALmsk_min,1) = Sal_min;
                plot(plottingDate, PSAL,'*-','MarkerSize',5,'LineWidth',1.2,'DisplayName','Salinity');%, 'Color', cl(i, :));
                hold on

                yyaxis right % alkalinity
                ALK = data_n{s}.Alkalinity(idx);
                ALKmsk_max = ALK > Talk_max; 
                ALK(ALKmsk_max,1) = Talk_max;
                ALKmsk_min = ALK < Talk_min; 
                ALK(ALKmsk_min,1) = Talk_min;
                ALK_flag = data_n{s}.TAFlag(idx);  

                if ALK_flag == 2
                    plot(plottingDate, ALK,'*-','MarkerSize',5,'LineWidth',1.2,'DisplayName','Alkalinity');%, 'Color', cl(i, :));
                    hold on
                else
                    disp('nearest bottle alkalinity data has bad QC');
                end

            else
                disp('nearest bottle more than 10m away')
     end
end
end
grid on
yyaxis left
ylim([Sal_min , Sal_max]);
yyaxis right
ylim([Talk_min, Talk_max]);
t = title({deployment_code ;'\rm \color[rgb]{0 0.4470 0.7410}Salinity  \color[rgb]{0.8500 0.3250 0.0980}Alkalinity'; '\rm \color[rgb]{black}QC flag 1 = line, QC flag 2 = circle, QC flag 3 = diamond, QC flag 4 = square'},'Interpreter','tex');
%t = title({deployment_code ;'\rm \color[rgb]{0 0.4470 0.7410}NO_{x}   \color[rgb]{0.8500 0.3250 0.0980}Phosphate   \color[rgb]{0.9290 0.6940 0.1250}Silicate'; '\rm \color[rgb]{black}QC flag 1 = line, QC flag 2 = circle, QC flag 3 = diamond, QC flag 4 = square'},'Interpreter','tex');

figurename = [deployment_code '-RAS-Alk-sal-report_plot.png'];
figure_path = 'C:\Users\cawynn\OneDrive - University of Tasmania\Documents\GitHub\RAS_data_plotting\RAS_data_plotting\figures'
cd(figure_path)
print(fig, '-dpng', figurename)


