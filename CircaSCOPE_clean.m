%%%%CircaSCOPE main script%%%%
%Last update 4.3.2021, developed on MATLAB R2020b by Gal Manella

%script input should be csv output files from the CellProfiler pipeline of
%CircaSCOPE, one file per condition/well
%The script analyses the data to reconstruct PTCs and additional parameters,
%as well as produces plots of the results.

%% script body
%required additional functions from FileExchange:
addpath('CircStat2012a');
addpath('costomcolormap');
addpath('rsquared');
addpath('harmfit');
addpath('raacampbell-sigstar-c1927a6');

%define the time constants for analysis:
FirstTP = 1;
LastTP = 216; 
bef_window = [24 96]; 
aft_window = [120 192];
cueT = 108;

well_labels = {'UT','Dex 100'}; %user-defined label for each well. should be the same order as PATH
PATH = { 'D:\galman.WISMAIN\Documents\CellProfiler\exp324\output_exp324_I_B6.csv';
    'D:\galman.WISMAIN\Documents\CellProfiler\exp324\output_exp324_II_B6.csv' }; %user-defined location of the input files (CellProfiler .csv outputs, one per well). 
OutDir = 'C:\Users\galman.WISMAIN\OneDrive\M.Sc\MSc Project\exp324\'; %where to save output figures?
PREFIX = 'exp324_Dex_Demo'; %for output filenames

tic
for i = 1:length(PATH) 
DataStruct(i) = LoadDATA(PATH{i},FirstTP,LastTP,cueT,0);
end
toc

%% PTC construction
%define concentration of treatment in each well, same order as PATH
Conc = [0 100];
UNIT = 'nM'; % '%', '\muM', 'nM', 'mM'
%Conc and UNIT are defined for plotting purposes

%Run the main PTC constructing function per well:
for i = 1:length(PATH)
PTCStruct(i) = PTC_construct(DataStruct(i).DataTable_Z,DataStruct(i).DataTable_unstacked,bef_window,aft_window,cueT,1);
end

%save the output structure with additional metadata:
for i = 1:length(PATH)
    PTCStruct(i).PATH = PATH{i};
    PTCStruct(i).well_label = well_labels{i};
    PTCStruct(i).OutDir = OutDir;
    PTCStruct(i).Conc = Conc(i);
    PTCStruct(i).UNIT = UNIT;
end 
save([PREFIX,'PTCStruct.mat'],'PTCStruct','-v7.3');

%define number of columns and rows in subplots throughout the script:
subcol = 2;
subrow = 1;
%% load existing PTCStruct
%this enable to reload previous analyses without running LoadData and
%PTC_construct again

load([PREFIX,'PTCStruct.mat']);
for i = 1:length(PATH)
    Conc(i)= PTCStruct(i).Conc;
end
UNIT = PTCStruct(1).UNIT;

subcol = 2;
subrow = 1;
%% Plot all PTCs

r_th = 0.5; %threshold for r-squared of fits before or after, a measure of rhythmicity robusteness
per_low = 22.5; %lower period length limit (in hours)
per_high = 28.5; %upper period length limit (in hours)

figure;
for i = 1:length(PATH)
    conds = find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th & PTCStruct(i).Periods_aft > per_low & PTCStruct(i).Periods_aft < per_high);
    subplot(subrow,subcol,i);
    plotDoublePTC(PTCStruct(i).TIPA_null_scaled(conds),...
    PTCStruct(i).fit_peaks_aft_scaled(conds),1,...
    [0,0,0],well_labels{i},1)
end

%% Plot all PRCs
r_th = 0.5;
per_low = 22; 
per_high = 29;

figure;
for i = 1:length(PATH)
       
    conds = find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th & PTCStruct(i).Periods_aft > per_low & PTCStruct(i).Periods_aft < per_high);
    subplot(subrow,subcol,i);
    plotDoublePRC(PTCStruct(i).TIPA_null_scaled(conds),...
        PTCStruct(i).ph_sh_TIPA_scaled(conds),...
        1,[0,0,0],well_labels{i},1)

end


%% Bootstrap with random sampling from control
%this module test for significant phase-shift using bootstrapping approach.
%the sampled parameter is the maximal absolute phase-shift (denoted as Amp
%in script), and it is resampled from a control well (ci is the index of
%control well) n_iter times. then the measured value is compared to the
%distribution of the sampled values and a p_val is calculated per well.

ci = 1;
conds1 = find(PTCStruct(ci).rsquare_bef > r_th & PTCStruct(ci).rsquare_aft > r_th & PTCStruct(ci).Periods_aft > per_low & PTCStruct(ci).Periods_aft < per_high);
n_iter = 1000;
    figure;

    p_val = NaN(length(PATH),1);
    
   
for i = 1:length(PATH)
    conds = find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th & PTCStruct(i).Periods_aft > per_low & PTCStruct(i).Periods_aft < per_high);
    k = length(conds);
    
    ph_bef1 = [PTCStruct(i).TIPA_null_scaled(conds)-1; PTCStruct(i).TIPA_null_scaled(conds); PTCStruct(i).TIPA_null_scaled(conds)+1];
    ph_sh1 = [PTCStruct(i).ph_sh_TIPA_scaled(conds);PTCStruct(i).ph_sh_TIPA_scaled(conds);PTCStruct(i).ph_sh_TIPA_scaled(conds)];
    ph_aft1 = [PTCStruct(i).fit_peaks_aft_scaled(conds);PTCStruct(i).fit_peaks_aft_scaled(conds);PTCStruct(i).fit_peaks_aft_scaled(conds)];
    
    ph_bef1 = ph_bef1(~isnan(ph_bef1) & ~isnan(ph_sh1));
    ph_sh1 = ph_sh1(~isnan(ph_bef1) & ~isnan(ph_sh1));
    ph_aft1 = ph_aft1(~isnan(ph_bef1) & ~isnan(ph_sh1));
    [ph_bef, IX] = sort(ph_bef1);
    ph_sh = ph_sh1(IX);
    ph_aft = ph_aft1(IX);      
    
    Amp(i) = fitfour1(ph_bef,ph_sh); 
    
    
     for j = 1:n_iter
    y = datasample([PTCStruct(ci).TIPA_null_scaled(conds1),PTCStruct(ci).ph_sh_TIPA_scaled(conds1)],k);
    ph_bef1 = [y(:,1)-1; y(:,1); y(:,1)+1];
    ph_sh1 = [y(:,2);y(:,2);y(:,2)];
    
    ph_bef1 = ph_bef1(~isnan(ph_bef1) & ~isnan(ph_sh1));
    ph_sh1 = ph_sh1(~isnan(ph_bef1) & ~isnan(ph_sh1));
    [ph_bef, IX] = sort(ph_bef1);
    ph_sh = ph_sh1(IX);       
        
     Amps_sampled(j) = fitfour1(ph_bef,ph_sh);     
        
     end
    
    
     Amps_sampled_table{i} = Amps_sampled;
     
    subplot(subrow,subcol,i);
    histogram(Amps_sampled, 'Normalization', 'probability','FaceColor',[0.5 0.5 0.5], 'BinWidth', 0.005,'EdgeColor', 'none');
    hold on
    xline(Amp(i),'Color','r','LineWidth',1.5,'LineStyle','--');
    p_val(i) = length(find(Amps_sampled >= Amp(i)))./n_iter;
    title([well_labels{i},'; p=',num2str(p_val(i))],'FontSize',7,'FontName','Arial');
    h = gca; h.YAxis.Visible = 'off'; 
    set(gca,'FontSize',7,'FontName','Arial');
    xlabel('max|\Delta\phi|');
    box off
end

%% PTC Classification by fourier fits

mypink = [208,28,139]./255;
mygreen = [77,172,38]./255;

    rmse0 = NaN(length(Conc),1);
    rmse1 = NaN(length(Conc),1);
    ph_bef = cell(length(PATH),1);
    ph_sh = cell(length(PATH),1);
    ph_aft = cell(length(PATH),1);
    Y1 = cell(length(PATH),1);
    Y_absph = cell(length(PATH),1);
    ph_aft_calc = cell(length(PATH),1);
    Y0 = cell(length(PATH),1);
    Y0_ph_sh = cell(length(PATH),1);
    ci1_lo = cell(length(PATH),1);
    ci1_up = cell(length(PATH),1);
    ci0_lo = cell(length(PATH),1);
    ci0_up = cell(length(PATH),1);
    
for i = 1:length(PATH)
    conds = find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th & PTCStruct(i).Periods_aft > per_low & PTCStruct(i).Periods_aft < per_high);

    ph_bef1 = [PTCStruct(i).TIPA_null_scaled(conds)-1; PTCStruct(i).TIPA_null_scaled(conds); PTCStruct(i).TIPA_null_scaled(conds)+1];
    ph_sh1 = [PTCStruct(i).ph_sh_TIPA_scaled(conds);PTCStruct(i).ph_sh_TIPA_scaled(conds);PTCStruct(i).ph_sh_TIPA_scaled(conds)];
    ph_aft1 = [PTCStruct(i).fit_peaks_aft_scaled(conds);PTCStruct(i).fit_peaks_aft_scaled(conds);PTCStruct(i).fit_peaks_aft_scaled(conds)];
    
    ph_bef1 = ph_bef1(~isnan(ph_bef1) & ~isnan(ph_sh1));
    ph_sh1 = ph_sh1(~isnan(ph_bef1) & ~isnan(ph_sh1));
    ph_aft1 = ph_aft1(~isnan(ph_bef1) & ~isnan(ph_sh1));
    [ph_bef{i}, IX] = sort(ph_bef1);
    ph_sh{i} = ph_sh1(IX);
    ph_aft{i} = ph_aft1(IX);
    
    if length(ph_aft1) < 7 || length(ph_bef1) < 7
        continue
    end
    
   
    mean_ph_aft = mod(circ_mean(PTCStruct(i).fit_peaks_aft_scaled(conds)*2*pi)/(2*pi),1);
    ph_aft_adj = mod(ph_aft{i} - mean_ph_aft + 0.5,1);
    
    [fitres_type1,gof1] = fit(ph_bef{i}, ph_sh{i},'fourier2','Upper',[Inf,Inf,Inf,Inf,Inf,2*pi]);
    
    [ci1,Y1{i}] = predint(fitres_type1,ph_bef{i},0.95,'observation','off');
    Y_absph{i} = Shift2phaseAft(ph_bef{i},Y1{i},1,0);
    ph_aft_calc{i} = Shift2phaseAft(ph_bef{i},ph_sh{i},1,0);
    ci1_lo{i} = ci1(:,1);
    ci1_up{i} = ci1(:,2);
    
    [fitres_type0,gof0] = fit(ph_bef{i}, ph_aft_adj,'fourier2','Upper',[Inf,Inf,Inf,Inf,Inf,2*pi]);   
    
    [ci0,Y0_temp] = predint(fitres_type0,ph_bef{i},0.95,'observation','off');
    Y0{i} = Y0_temp + mean_ph_aft - 0.5; %convert back to real ph_aft units
    
    Y0_ph_sh{i} = phaseShift(ph_bef{i},Y0{i},1,1,1);
    ci0_lo{i} = ci0(:,1) + mean_ph_aft - 0.5;
    ci0_up{i} = ci0(:,2) + mean_ph_aft - 0.5;
    
    rmse0(i) = gof0.rmse;
    rmse1(i) = gof1.rmse;
end
    


figure('Units','inches','Position',[4,4,7.5,5.625]);
for i = 1:length(PATH)
    if length(ph_aft{i}) < 7 || length(ph_bef{i}) < 7
        continue
    end
    
    subplot(subrow,subcol,i);
    plotDoublePTC(ph_bef{i}(ph_bef{i} >0 & ph_bef{i} <= 1),ph_aft_calc{i}(ph_bef{i} >0 & ph_bef{i} <= 1),...
    1,[0.7 0.7 0.7],[well_labels{i}],1)
    set(gca,'FontSize',8,'FontName','Arial');
    hold on;
    
    %plot type-1 fit
    patch([ph_bef{i};flipud(ph_bef{i})],[ci1_lo{i}+ph_bef{i};flipud(ci1_up{i}+ph_bef{i})],mypink,'FaceAlpha',0.25,'EdgeColor','none');
    plot(ph_bef{i},Y1{i}+ph_bef{i},'Color',mypink,'LineWidth',1.5);   
    %plot type-0 fit
    patch([ph_bef{i};flipud(ph_bef{i})],[ci0_lo{i};flipud(ci0_up{i})],mygreen,'FaceAlpha',0.25,'EdgeColor','none');
    plot(ph_bef{i},Y0{i},'Color',mygreen,'LineWidth',1.5);    

    
end
print([OutDir,PREFIX,'PTCs_BOTH.emf'],'-dmeta','-painters');

figure('Units','inches','Position',[4,4,7.5,5.625]);
for i = 1:length(PATH)
    
    if length(ph_aft{i}) < 7 || length(ph_bef{i}) < 7
        continue
    end
    
    subplot(subrow,subcol,i);
    plotDoublePRC(ph_bef{i}(ph_bef{i} >0 & ph_bef{i} <= 1),ph_sh{i}(ph_bef{i} >0 & ph_bef{i} <= 1),...
    1,[0.7 0.7 0.7],[well_labels{i}],1)
    set(gca,'FontSize',8,'FontName','Arial');
    hold on;
    
    %plot type-1 fit
    patch([ph_bef{i};flipud(ph_bef{i})],[ci1_lo{i};flipud(ci1_up{i})],mypink,'FaceAlpha',0.25,'EdgeColor','none');
    plot(ph_bef{i},Y1{i},'Color',mypink,'LineWidth',1.5);   
    %plot type-0 fit
    patch([ph_bef{i};flipud(ph_bef{i})],[ci0_lo{i}-ph_bef{i};flipud(ci0_up{i}-ph_bef{i})],mygreen,'FaceAlpha',0.25,'EdgeColor','none');
    plot(ph_bef{i},Y0{i}-ph_bef{i},'Color',mygreen,'LineWidth',1.5);    
   
end
print([OutDir,PREFIX,'PRCs_BOTH.emf'],'-dmeta','-painters');


%% Bootsrap RMSE (model selection between type-1 and type-0)
%This module test the significance of the difference between RMSE of type-1
%and type-0 models. this is done with the inner function bootstrap_rmse, 
%based on resampling of the PTC data and each time fitting both models and
%retrieving the RMSE. p-value of type-0 (RMSE_p0) represents the 
%probability that RMSE0-RMSE1 is larger than 0. A similar p-value for
%type-1 represents the probablilty that RMSE0-RMSE1 is smaller than 0. 

n_iter = 1000;
bstrp_rmse = boostrap_rmse(PTCStruct,well_labels,r_th,per_low,per_high,n_iter,subrow,subcol);

RMSE_p0 = bstrp_rmse.RMSE_p0;
p0_th = 0.05;

%95% confidence interval for delta-RMSE:
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
rmse_dlt_obs = rmse0 - rmse1;
rmse_CI = cellfun(@(x) CIFcn(x,95),bstrp_rmse.rmse_dlt_table','UniformOutput',false);
rmse_CI = cell2mat(rmse_CI);

%%plot delta-RMSE graph:
x_vals = 1:length(PATH);
figure('Units','inches','Position',[5,5,3.5,1.5]);
yline(0,'-k');
hold on;
errorbar(x_vals,rmse_dlt_obs,rmse_dlt_obs-rmse_CI(:,1),rmse_CI(:,2)-rmse_dlt_obs,'Color','k','LineStyle','none');
scatter(x_vals(RMSE_p0 < p0_th), rmse_dlt_obs(RMSE_p0 < p0_th),20,mygreen,'filled');
scatter(x_vals(RMSE_p0 >= p0_th),rmse_dlt_obs(RMSE_p0 >= p0_th),20,mypink,'filled');
xlim([0 max(x_vals)+1]);
ylim([min(rmse_CI(:,1))-0.05,max(rmse_CI(:,2))+0.05]);
ylabel('RMSE_{0} - RMSE_{1}');
set(gca,'XTick',x_vals,'XTickLabel',well_labels,'XTickLabelRotation',0,'FontSize',8,'FontName','Arial');
print([OutDir,PREFIX,'rmse_dlt.emf'],'-dmeta');

%% plot PRCs according to RMSE analysis

figure('Units','inches','Position',[4,4,7.5,5.625]);
for i = 1:length(PATH)
     if length(ph_aft{i}) < 7 || length(ph_bef{i}) < 7
        continue
    end
    
    subplot(subrow,subcol,i);
    plotDoublePRC(ph_bef{i}(ph_bef{i} >0 & ph_bef{i} <= 1),ph_sh{i}(ph_bef{i} >0 & ph_bef{i} <= 1),...
    1,[0.7 0.7 0.7],[well_labels{i}],1)
    hold on;
    set(gca,'FontSize',8,'FontName','Arial');
    
    if RMSE_p0(i) >= p0_th
     patch([ph_bef{i};flipud(ph_bef{i})],[ci1_lo{i};flipud(ci1_up{i})],mypink,'FaceAlpha',0.25,'EdgeColor','none');
     plot(ph_bef{i},Y1{i},'Color',mypink,'LineWidth',1.5);
    else
     patch([ph_bef{i};flipud(ph_bef{i})],[ci0_lo{i}-ph_bef{i};flipud(ci0_up{i}-ph_bef{i})],mygreen,'FaceAlpha',0.25,'EdgeColor','none');   
     plot(ph_bef{i},Y0{i}-ph_bef{i},'Color',mygreen,'LineWidth',1.5);
    end
end

print([OutDir,PREFIX,'PRCs_',num2str(r_th),'.emf'],'-dmeta','-painters');

%% plot PTCs according to RMSE analysis

figure('Units','inches','Position',[4,4,7.5,5.625]);
for i = 1:length(PATH)
    
     if length(ph_aft{i}) < 7 || length(ph_bef{i}) < 7
        continue
    end
    
    subplot(subrow,subcol,i);
    plotDoublePTC(ph_bef{i}(ph_bef{i} >0 & ph_bef{i} <= 1),ph_aft_calc{i}(ph_bef{i} >0 & ph_bef{i} <= 1),...
    1,[0.7 0.7 0.7],[well_labels{i}],1)
    set(gca,'FontSize',8,'FontName','Arial');
    hold on;
    
    if RMSE_p0(i) >= p0_th
     patch([ph_bef{i};flipud(ph_bef{i})],[ci1_lo{i}+ph_bef{i};flipud(ci1_up{i}+ph_bef{i})],mypink,'FaceAlpha',0.25,'EdgeColor','none');
     plot(ph_bef{i},Y1{i}+ph_bef{i},'Color',mypink,'LineWidth',1.5);
    else
     patch([ph_bef{i};flipud(ph_bef{i})],[ci0_lo{i};flipud(ci0_up{i})],mygreen,'FaceAlpha',0.25,'EdgeColor','none');   
     plot(ph_bef{i},Y0{i},'Color',mygreen,'LineWidth',1.5);
    end
end

print([OutDir,PREFIX,'PTCs_',num2str(r_th),'.emf'],'-dmeta','-painters');


%% Feature extraction and plotting

%max/min phase shift, type 1 or 0
[max_ph_sh1,I_max_ph_sh1] = cellfun(@(x,y) max(x(y>=0 & y<=1)),Y1,ph_bef);
[min_ph_sh1,I_min_ph_sh1] = cellfun(@(x,y) min(x(y>=0 & y<=1)),Y1,ph_bef);
ph_bef_subs1 = cellfun(@(x) x(x>=0 & x<=1),ph_bef,'UniformOutput',false);

% [max_ph_sh0,I_max_ph_sh0] = cellfun(@max,Y0_ph_sh);
% [min_ph_sh0,I_min_ph_sh0] = cellfun(@min,Y0_ph_sh);
[max_ph_sh0,I_max_ph_sh0] = cellfun(@(x,y) max(x(x-y<=0.5 & x-y>=-0.5)-y(x-y<=0.5 & x-y>=-0.5)),Y0,ph_bef);
[min_ph_sh0,I_min_ph_sh0] = cellfun(@(x,y) min(x(x-y<=0.5 & x-y>=-0.5)-y(x-y<=0.5 & x-y>=-0.5)),Y0,ph_bef);
ph_bef_subs0 = cellfun(@(x,y) y(x-y<=0.5 & x-y>=-0.5),Y0,ph_bef,'UniformOutput',false);

max_ph_sh(rmse0 >= rmse1) = max_ph_sh1(RMSE_p0 >= p0_th);
max_ph_sh(rmse0 < rmse1) = max_ph_sh0(RMSE_p0 < p0_th);
min_ph_sh(rmse0 >= rmse1) = min_ph_sh1(RMSE_p0 >= p0_th);
min_ph_sh(rmse0 < rmse1) = min_ph_sh0(RMSE_p0 < p0_th);

%phase of maximal phase-shift, type1
ph_bef_max1 = cellfun(@(x,inx) mod(x(inx),1),ph_bef_subs1,num2cell(I_max_ph_sh1));
ph_bef_min1 = cellfun(@(x,inx) mod(x(inx),1),ph_bef_subs1,num2cell(I_min_ph_sh1));
ph_bef_max0 = cellfun(@(x,inx) mod(x(inx),1),ph_bef_subs0,num2cell(I_max_ph_sh0));
ph_bef_min0 = cellfun(@(x,inx) mod(x(inx),1),ph_bef_subs0,num2cell(I_min_ph_sh0));
ph_bef_max(rmse0 >= rmse1) =  ph_bef_max1(RMSE_p0 >= p0_th);
ph_bef_max(rmse0 < rmse1) =  ph_bef_max0(RMSE_p0 < p0_th);
ph_bef_min(rmse0 >= rmse1) =  ph_bef_min1(RMSE_p0 >= p0_th);
ph_bef_min(rmse0 < rmse1) =  ph_bef_min0(RMSE_p0 < p0_th);

%mean final phase of type 0

% [mean_Y0] = cellfun(@(x) mod(circ_mean(x*2*pi)/(2*pi),1),Y0);
[mean_Y0,cil_Y0,ciu_Y0] = cellfun(@(x) circ_mean(x*2*pi),Y0);
mean_Y0 = mod(mean_Y0/(2*pi),1);
cil_Y0 = mod(cil_Y0/(2*pi),1);
ciu_Y0 = mod(ciu_Y0/(2*pi),1);
max_Y0 = cellfun(@(x) max(x),Y0);
min_Y0 = cellfun(@(x) min(x),Y0);
std_Y0 = cellfun(@(x) circ_std(x*2*pi)/2*pi,Y0);


first_i = 1;
last_i = length(Conc);
x_vals = 1:(last_i - first_i + 1);

figure('Units','inches','Position',[5,5,3.5,1.5]);
plot(x_vals,max_ph_sh(first_i:last_i),'^-',x_vals,abs(min_ph_sh(first_i:last_i)),'v-','MarkerFaceColor','auto');
ylabel('Max. shift');
xlabel(['Conc. (',UNIT,')']); %(\muM)
set(gca,'XLim',[0.5 x_vals(end)+0.5],'XTick',x_vals,'XTickLabel',Conc(first_i:last_i),'FontSize',7,'FontName','Arial');
box off;
%legend('|max\Delta\phi|','|min\Delta\phi|','Location','southeast');
print([OutDir,PREFIX,'minmax_ph_sh.emf'],'-dmeta');

figure('Units','inches','Position',[5,5,3.5,1.5]);
plot(x_vals,ph_bef_max(first_i:last_i),'^-',x_vals,ph_bef_min(first_i:last_i),'v-','MarkerFaceColor','auto');
ylabel('Phase of max shift');
xlabel(['Conc. (',UNIT,')']); %(\muM)
set(gca,'XLim',[0.5 x_vals(end)+0.5],'XTick',x_vals,'XTickLabel',Conc(first_i:last_i),'FontSize',7,'FontName','Arial');
box off;
%legend('|max\Delta\phi|','|min\Delta\phi|','Location','southeast');
print([OutDir,PREFIX,'phase_of_max_sh.emf'],'-dmeta');

figure('Units','inches','Position',[5,5,3.5,1.5]);
bar([max_ph_sh;abs(min_ph_sh)]');
ylabel('Max. shift');
set(gca,'XTickLabel',well_labels,'XTickLabelRotation',0,'FontSize',8,'FontName','Arial');
box off;
%legend('|max\Delta\phi|','|min\Delta\phi|','Location','northwest');
print([OutDir,PREFIX,'minmax_ph_sh_BAR.emf'],'-dmeta');

figure('Units','inches','Position',[5,5,2.5,1.5]);
bar([max_ph_sh;abs(min_ph_sh)]');
ylabel('Max. shift');
set(gca,'XTickLabel',well_labels,'XTickLabelRotation',0,'FontSize',8,'FontName','Arial');
box off;
%legend('|max\Delta\phi|','|min\Delta\phi|','Location','northwest');
xlim([1.5,x_vals(end)+0.5]);
print([OutDir,PREFIX,'minmax_ph_sh_BAR_compact.emf'],'-dmeta');

figure('Units','inches','Position',[5,5,3.5,1.5]);
bar([ph_bef_max;ph_bef_min]');
ylabel('Phase of max shift');
set(gca,'XTickLabel',well_labels,'XTickLabelRotation',0,'FontSize',8,'FontName','Arial');
box off;
%legend('|max\Delta\phi|','|min\Delta\phi|','Location','northwest');
print([OutDir,PREFIX,'phase_of_max_sh_BAR.emf'],'-dmeta');

figure('Units','inches','Position',[5,5,2.5,1.5]);
bar([ph_bef_max;ph_bef_min]');
ylabel('Phase of max shift');
set(gca,'XTickLabel',well_labels,'XTickLabelRotation',0,'FontSize',8,'FontName','Arial');
box off;
%legend('|max\Delta\phi|','|min\Delta\phi|','Location','northwest');
xlim([1.5,x_vals(end)+0.5]);
print([OutDir,PREFIX,'phase_of_max_sh_BAR_compact.emf'],'-dmeta');


%% type-specific plots
first_i = 1;
last_i = 1; %length(Conc);
x_vals = 1:(last_i - first_i + 1);

figure('Units','inches','Position',[5,5,2,1.5]);
plot(x_vals,ph_bef_max(first_i:last_i),'^-',x_vals,ph_bef_min(first_i:last_i),'v-','MarkerFaceColor','auto');
ylim([0 1]);
ylabel('Phase of max shift');
xlabel(['Conc. (',UNIT,')']); %(\muM)
set(gca,'XLim',[0.5 x_vals(end)+0.5],'XTick',x_vals,'XTickLabel',Conc(first_i:last_i),'FontSize',7,'FontName','Arial');
box off;
ax = gca;
print([OutDir,PREFIX,'phase_of_max_sh_type1.emf'],'-dmeta');


first_i = 2;
last_i = 2;
x_vals = 1:(last_i - first_i + 1);

figure('Units','inches','Position',[5,5,1.5,1.5]);
errorbar(x_vals,mean_Y0(first_i:last_i),mean_Y0(first_i:last_i)-min_Y0(first_i:last_i),max_Y0(first_i:last_i)-mean_Y0(first_i:last_i),'.-','MarkerSize',10);
ylabel('Mean final phase');
xlabel(['Conc. (',UNIT,')']); %(\muM)
ylim([0 1]);
set(gca,'XLim',[0.5 x_vals(end)+0.5],'XTick',x_vals,'XTickLabel',Conc(first_i:last_i),'FontSize',7,'FontName','Arial');
box off;
print([OutDir,PREFIX,'mean_ph_aft_type0.emf'],'-dmeta');


%% 3D PTC by binning
%relevant for concentration gradients. plot the PTCs as a phase plane, with
%color representing the New Phase and y-axis the cue concentration/intensity.

clear('binned_aft');
clear('binned_bef');
intvl = 0.05; %bin width
bef_times = -0.5:intvl:1.5+intvl; %define the x-axis

for i = 1:length(Conc)
    for j = 1:length(bef_times)-1
        binned_aft(i,j) = mod(circ_median(ph_aft_calc{i}(ph_bef{i} > bef_times(j+1) - (bef_times(j+1)-bef_times(j)) & ph_bef{i} <= bef_times(j+1))*2*pi)/(2*pi),1);
        binned_bef(i,j) = (bef_times(j+1)-bef_times(j))./2;

    end
end


figure('Units','inches','Position',[5,5,3.5,1.5]);
imagesc(flip(binned_aft));
caxis([0,1]);
colormap(hsv);
xlabel('Old Phase');
ylabel('Conc. (\muM)');
ylim([1.5 length(Conc)+0.5]);
set(gca,'YTick',2:length(Conc),'YTickLabel',flip(Conc(2:end)),...
    'XTick',1:2:length(bef_times),'XTickLabel',-0.5:intvl*2:1.5, 'FontName','Arial','FontSize',7,'XTickLabelRotation',90);
c = colorbar;
c.Label.String = 'New Phase';
print([OutDir,PREFIX,'3d_binned.emf'],'-dmeta');

    
%% sorted heatmap
%requires DateStruct in the workspace

%plot a heatmap where each row correspond to one cell, sorted by initial
%phase.

    figure;

for i = 1:length(PATH) 
    subplot(2,3,i)
    conds = find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th);

    datatable = DataStruct(i).DataTable_norm';
    datafiltered = datatable(conds,:);
    [phasesorted,IX] = sort(PTCStruct(i).fit_peaks_bef(conds));
    datasorted = datafiltered(IX,:);

    mycmap = customcolormap_preset('pasteljet');
  
    imagesc(log2(datasorted));
    colormap(mycmap);
    caxis([-1 1]);
    title(well_labels{i});
end

%% survival curves:
%requires DateStruct in the workspace

%%calculate and plot different statistics regarding tracking efficiency,
%%cell survival and cell rhythmicity

figure;
for i = 1:length(PATH)
    obj_per_timepoint = [];
    for j = 1:max(DataStruct(i).TIME)
        obj_per_timepoint(j) = length(find(DataStruct(i).TIME == j)); 
    end
    
    subplot(subrow,subcol,i);
    plot(1:max(DataStruct(i).TIME),obj_per_timepoint);
    title(well_labels{i});
    xlabel('Time');
    xlim([0 max(DataStruct(i).TIME)]);
    %ylim([4000 8000]);
    disp([num2str(i),'. full track obj: ',num2str(length(find(DataStruct(i).Obj_label_life >= LastTP)))]);
    disp([num2str(i),'. full track obj/first TP: ',num2str(length(find(DataStruct(i).Obj_label_life >= LastTP))),'/',num2str(length(find(DataStruct(i).Obj_label_first == FirstTP))),'(',num2str(length(find(DataStruct(i).Obj_label_life >= LastTP))/length(find(DataStruct(i).Obj_label_first == FirstTP))),')']);
    disp([num2str(i),'. full track obj/last TP: ',num2str(length(find(DataStruct(i).Obj_label_life >= LastTP))),'/',num2str(length(find(DataStruct(i).Obj_label_last >= LastTP))),'(',num2str(length(find(DataStruct(i).Obj_label_life >= LastTP))/length(find(DataStruct(i).Obj_label_last >= LastTP))),')']);

end
suptitle('Cells per timepoint');
print([OutDir,PREFIX,'Cells per timepoint.emf'],'-dmeta');

figure;
for i = 1:length(PATH)
    subplot(subrow,subcol,i);
    histogram(DataStruct(i).Obj_label_life);
    title(well_labels{i});
    xlabel('Age (hours)');    
end
suptitle('Cell lifespan');

figure
for i = 1:length(PATH)
    subplot(subrow,subcol,i);
    ecdf(DataStruct(i).Obj_label_life(DataStruct(i).Obj_label_first == FirstTP),'Function','survivor');
    title(well_labels{i});
    xlabel('Age (hours)');
end
suptitle('Cumulative survival curve');
print([OutDir,PREFIX,'survival.emf'],'-dmeta');


%summary table

for i = 1:length(PATH)
   cells_1st_tp(i) = length(find(DataStruct(i).Obj_label_first == FirstTP));
   cells_full_track(i) = length(find(DataStruct(i).Obj_label_life >= LastTP));
   mean_lifespan_per_well(i) = mean(DataStruct(i).Obj_label_life);
   sd_lifespan_per_well(i) = std(DataStruct(i).Obj_label_life);
   mean_lifespan_1st_tp_per_well(i) = mean(DataStruct(i).Obj_label_life(DataStruct(i).Obj_label_first == FirstTP));
   sd_lifespan_1st_tp_per_well(i) = std(DataStruct(i).Obj_label_life(DataStruct(i).Obj_label_first == FirstTP));
   cells_rhythmic(i) = length(find(PTCStruct(i).rsquare_bef > r_th)); %& PTCStruct(i).rsquare_aft > r_th);
   cells_rhythmic_period(i) = length(find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th & PTCStruct(i).Periods_aft > per_low & PTCStruct(i).Periods_aft < per_high));
end

mean_cells_1st_tp = mean(cells_1st_tp);
sd_cells_1st_tp = std(cells_1st_tp);
mean_cells_full_track = mean(cells_full_track);
sd_cells_full_track = std(cells_full_track);
mean_cells_rhythmic = mean(cells_rhythmic);
sd_cells_rhythmic = std(cells_rhythmic);

frac_full_outof_1st = cells_full_track ./ cells_1st_tp;
mean_frac_full_outof_1st = mean(frac_full_outof_1st);
sd_frac_full_outof_1st = std(frac_full_outof_1st);

frac_rhythmic_outof_1st = cells_rhythmic ./ cells_1st_tp;
mean_frac_rhythmic_outof_1st = mean(frac_rhythmic_outof_1st);
sd_frac_rhythmic_outof_1st = std(frac_rhythmic_outof_1st);

sprintf("Cells in 1st timepoint: %.1f%s%.1f",mean_cells_1st_tp,'±',sd_cells_1st_tp)
sprintf("Fully trackable cells: %.1f%s%.1f (%.2g%s%.2g)",mean_cells_full_track,'±',sd_cells_full_track,mean_frac_full_outof_1st,'±',sd_frac_full_outof_1st)
sprintf("Rhythmic cells: %.1f%s%.1f (%.2g%s%.2g)",mean_cells_rhythmic,'±',sd_cells_rhythmic,mean_frac_rhythmic_outof_1st,'±',sd_frac_rhythmic_outof_1st)

sprintf(" Cells in 1st timepoint: \t %.1f%s%.1f \n Fully trackable cells: \t %.1f%s%.1f (%.2g%s%.2g) \n Rhythmic cells: \t\t\t %.1f%s%.1f (%.2g%s%.2g) \n" ...
    ,mean_cells_1st_tp,'±',sd_cells_1st_tp ...
    ,mean_cells_full_track,'±',sd_cells_full_track,mean_frac_full_outof_1st,'±',sd_frac_full_outof_1st ...
    ,mean_cells_rhythmic,'±',sd_cells_rhythmic,mean_frac_rhythmic_outof_1st,'±',sd_frac_rhythmic_outof_1st)


%% periods before and after
%Compare the period lengths before and after for each well

myred = [203,24,29]./255;
myblue = [8,69,148]./255;

lb_rot = 0; % rotation of x-label if necessary

figure;
for i = 1:length(PATH)
    subplot(subrow,subcol,i)
    conds = find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th);
    histogram(PTCStruct(i).Periods_bef(conds),12);
    hold on;
    histogram(PTCStruct(i).Periods_aft(conds),12);
    title(well_labels{i});
    ylabel('Freq.');
    xlabel('Period (h)');
    colororder([myblue;myred]);
    set(gca,'FontSize',8,'FontName','Arial');
end
legend('Before','After');

conds = find(PTCStruct(1).rsquare_bef > r_th & PTCStruct(1).rsquare_aft > r_th);
periods = [PTCStruct(1).Periods_bef(conds);PTCStruct(1).Periods_aft(conds)];
Groups1 = repmat(1,length(PTCStruct(1).Periods_bef(conds))*2,1);
Groups2 = [repmat(1,length(PTCStruct(1).Periods_bef(conds)),1);repmat(2,length(PTCStruct(1).Periods_bef(conds)),1)];


for i = 2:length(PATH)
   conds = find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th);
   periods = [periods;PTCStruct(i).Periods_bef(conds);PTCStruct(i).Periods_aft(conds)];
   Groups1 = [Groups1;repmat(i,length(PTCStruct(i).Periods_bef(conds))*2,1)];
   Groups2 = [Groups2;repmat(1,length(PTCStruct(i).Periods_bef(conds)),1);repmat(2,length(PTCStruct(i).Periods_bef(conds)),1)];
 
end

figure('Units','inches','Position',[5,5,3.5,2]);
boxplot(periods,[Groups1,Groups2],'ColorGroup',Groups2,'PlotStyle','compact','FactorGap',0,'Colors', [myblue;myred], 'FactorSeparator',1);
ylabel('Period (h)');
%legend('Before','After');
box off;
set(gca,'XTick', 1.5:2:length(PATH)*2,'XTickLabel',well_labels,'XTickLabelRotation',lb_rot,'FontSize',8,'FontName','Arial');
set(gca,'XTickLabel',well_labels);
print([OutDir,PREFIX,'periods_box.emf'],'-dmeta');


figure('Units','inches','Position',[5,5,3.5,2]);
boxplot(periods,[Groups1,Groups2],'ColorGroup',Groups2,'PlotStyle','compact','FactorGap',0,'Colors', [myblue;myred], 'FactorSeparator',1);
ylabel('Period (h)');
%legend('Before','After');
box off;
set(gca,'XTick', 1.5:2:length(PATH)*2,'XTickLabel',Conc,'XTickLabelRotation',0,'FontSize',7,'FontName','Arial');
xlabel(['Conc. (',UNIT,')']);
%xlim([2.5,length(PATH)*2+0.5]); %remove if u want the first condition
print([OutDir,PREFIX,'periods_box_cont.emf'],'-dmeta');


dlt_periods = NaN(length(periods),length(PATH));
periods_aft = NaN(length(periods),length(PATH));
for i= 1:length(PATH)
    conds = find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th);
    dlt_periods(1:length(PTCStruct(i).Periods_bef(conds)),i) = PTCStruct(i).Periods_aft(conds) - PTCStruct(i).Periods_bef(conds);
    periods_aft(1:length(PTCStruct(i).Periods_aft(conds)),i) = PTCStruct(i).Periods_aft(conds);
end

p_periods = ones(length(PATH),1);
p_periods_aft = ones(length(PATH),1);
for i= 2:length(PATH)
[~,p_periods(i)] = ttest2(dlt_periods(:,1),dlt_periods(:,i));
[~,p_periods_aft(i)] = ttest2(periods_aft(:,1),periods_aft(:,i));
end

figure('Units','inches','Position',[5,5,3.5,2]);
yline(0,'Color','k');
hold on
boxplot(dlt_periods,'PlotStyle','compact','FactorGap',0,'Colors',[0.5 0.5 0.5]);
set(gca,'XTick', 1:length(PATH),'XTickLabel',well_labels,'XTickLabelRotation',lb_rot,'FontSize',8,'FontName','Arial');
ylabel('Period_{aft} - Period_{bef}');
%xlim([1.5,length(PATH)+0.5]); %remove if u want the first condition
box off
G = {[1,2]};
sigstar(G,p_periods(2:end)');
print([OutDir,PREFIX,'dlt_periods_box_sig.emf'],'-dmeta');


figure('Units','inches','Position',[5,5,3.5,2]);
yline(0,'Color','k');
hold on
boxplot(dlt_periods,'PlotStyle','compact','FactorGap',0,'Colors',[0.5 0.5 0.5]);
set(gca,'XTick', 1:length(PATH),'XTickLabel',Conc,'XTickLabelRotation',0,'FontSize',7,'FontName','Arial');
ylabel('Period_{aft} - Period_{bef}');
xlabel(['Conc. (',UNIT,')']); 
%xlim([1.5,length(PATH)+0.5]); %remove if u want the first condition
box off


print([OutDir,PREFIX,'dlt_periods_box_cont.emf'],'-dmeta');


figure('Units','inches','Position',[5,5,3.5,2]);
boxplot(periods,[Groups1,Groups2],'ColorGroup',Groups2,'PlotStyle','compact','FactorGap',0,'Colors', [myblue;myred], 'FactorSeparator',1);
ylabel('Period (h)');
%legend('Before','After');
box off;
set(gca,'XTick', 1.5:2:length(PATH)*2,'XTickLabel',well_labels,'XTickLabelRotation',0,'FontSize',8,'FontName','Arial');
set(gca,'XTickLabel',well_labels);
G = {[2,4]};
sigstar(G,p_periods_aft(2:end)');
print([OutDir,PREFIX,'periods_box_sig.emf'],'-dmeta');


%% inner functions

function OUT = LoadDATA(PATH,FirstTP,LastTP,CueT,OmitCueTFlag)
%%%%%load data:%%%%%
 
    DATA = readtable(PATH);

    ImageNumber = DATA.ImageNumber;
    ObjectNumber = DATA.ObjectNumber;
    
    TrackObjects_Label = DATA.TrackObjects_Label_50;
    TrackObjects_ParentImageNumber =DATA.TrackObjects_ParentImageNumber_50;
    TrackObjects_ParentObjectNumber =DATA.TrackObjects_ParentObjectNumber_50;

    Intensity_IntegratedIntensity_CorrGreen =DATA.Intensity_IntegratedIntensity_CorrGreen; %N(:,24);
    Intensity_MeanIntensity_CorrGreen =DATA.Intensity_MeanIntensity_CorrGreen; %N(:,45);
    Intensity_MeanIntensity_CorrGreenBgSub =DATA.Intensity_MeanIntensity_CorrGreenBgSub; 

    Children = ImageNumber*10000 + ObjectNumber;
    Parents = TrackObjects_ParentImageNumber*10000 + TrackObjects_ParentObjectNumber;

    %group/frame
    Metadata_Vessel = DATA.Metadata_Vessel;
    Metadata_Site = DATA.Metadata_Site;
    Metadata_Well = DATA.Metadata_Well;
    GROUP = string([num2str(Metadata_Site),char(Metadata_Well)]);
    GROUP_list = unique(GROUP,'stable');

    GROUP_n = NaN(length(GROUP),1);

    for i = 1:length(GROUP)
       GROUP_n(i) = find(strcmp(GROUP_list,GROUP(i))); 
    end

    Well_list = unique(Metadata_Well,'stable');
    Well_n = NaN(length(Metadata_Well),1);
    for i = 1:length(GROUP)
       Well_n(i) = find(strcmp(Well_list,Metadata_Well(i))); 
    end

    %timestamp
    Metadata_Year = DATA.Metadata_Year;
    Metadata_Month = DATA.Metadata_Month;
    Metadata_Day = DATA.Metadata_Day;
    Metadata_Hour = DATA.Metadata_Hour;
    Metadata_Min = DATA.Metadata_Min;
    Metatext_Year =  num2str(Metadata_Year);
    Metatext_Month = num2str(Metadata_Month);
    if size(Metatext_Month,2) == 1 %for the case all Month<10 for example. add a similar condition to other metadata as needed
        clear Metatext_Month;
        Metatext_Month(:,2) = num2str(Metadata_Month);
    end
    Metatext_Month(Metadata_Month <10,1) = '0';
    Metatext_Day = num2str(Metadata_Day);
    Metatext_Day(Metadata_Day <10,1) = '0';
    Metatext_Hour = num2str(Metadata_Hour);
    Metatext_Hour(Metadata_Hour <10,1) = '0';
    Metatext_Min = num2str(Metadata_Min);
    if size(Metatext_Min,2) == 1 %for the case all Min==0 for example. add a similar condition to other metadata as needed
        clear Metatext_Min;
        Metatext_Min(:,2) = num2str(Metadata_Min);
    end
    Metatext_Min(Metadata_Min <10,1) = '0';

    Location_Center_X = DATA.Location_Center_X;
    Location_Center_Y = DATA.Location_Center_Y;

    Metatext_TIME = string([Metatext_Year,Metatext_Month,Metatext_Day,Metatext_Hour,Metatext_Min]);

    TIME_list = unique(Metatext_TIME,'stable');
    TIME = NaN(length(Metatext_TIME),1);

    %added 10/2020: CP skip the CueT timepoint because of focus problems,
    %so the MATLAB script need to cope with that:
    if OmitCueTFlag
    TIME_list = {TIME_list{1:CueT},'XXX',TIME_list{CueT+1:end}}; %to allow skipping the CueT, without having a lag as a result
    end 
    
    for i = 1:length(Metatext_TIME)
       TIME(i) = find(strcmp(TIME_list,Metatext_TIME(i))); 
    end
%%%%%Track Labels%%%%%

    Obj_label = TrackObjects_Label + GROUP_n*10000;
    unique_table = table(unique(Obj_label(~isnan(Obj_label))),ones(length(unique(Obj_label(~isnan(Obj_label)))),1)); % doesn't have to be table, it is for compitability with the current script

    Obj_label_life = NaN(length(unique_table{:,1}),1);
    Obj_label_first = NaN(length(unique_table{:,1}),1);
    Obj_label_length = NaN(length(unique_table{:,1}),1);
    Obj_label_last = NaN(length(unique_table{:,1}),1);
    Obj_well = cell(length(unique_table{:,1}),1);
    Obj_well_n = NaN(length(unique_table{:,1}),1);
    Obj_GROUP = NaN(length(unique_table{:,1}),1);
    for i = 1:length(unique_table{:,1})
       
       temp_TIME = TIME(Obj_label == unique_table{i,1});
       Obj_label_first(i) = temp_TIME(1);
       Obj_label_last(i) = temp_TIME(end);
       Obj_label_life(i) = Obj_label_last(i)- Obj_label_first(i)+ 1;
       Obj_label_length(i) = length(find(Obj_label == unique_table{i,1}));

       Obj_well(i) = Metadata_Well(find(Obj_label == unique_table{i,1},1));
       Obj_well_n(i) = Well_n(find(Obj_label == unique_table{i,1},1));
       Obj_GROUP(i) = GROUP_n(find(Obj_label == unique_table{i,1},1));
    end

%%%%%Subset and normalize%%%%%

    ind = find(ismember(Obj_label,unique_table{Obj_label_first <=FirstTP & Obj_label_last >=LastTP,1}));

    DataTable = table(TIME(ind),Obj_label(ind),Intensity_MeanIntensity_CorrGreenBgSub(ind));

    DataTable_unstacked = unstack(DataTable,3,2);
    DataTable_norm = DataTable_unstacked{:,2:end} ./ mean(DataTable_unstacked{:,2:end},1,'omitnan');
    DataTable_Z = (DataTable_unstacked{:,2:end} - mean(DataTable_unstacked{:,2:end},1,'omitnan')) ./ std(DataTable_unstacked{:,2:end},[],1,'omitnan');

    %%%if Z should be detrended:
    DataTable_trend = zeros(size(DataTable_norm));
    for i = 1:size(DataTable_Z,2)
    DataTable_trend(:,i) = smooth(DataTable_unstacked{:,i+1},24);
    end
    DataTable_Z = (DataTable_unstacked{:,2:end} - DataTable_trend) ./ std(DataTable_unstacked{:,2:end},[],1,'omitnan');
       
    
%%%%%OUTPUT%%%%%

    OUT.ImageNumber = ImageNumber; %?
    OUT.ObjectNumber = ObjectNumber; %?
    OUT.TrackObjects_Label = TrackObjects_Label; %?
    OUT.TrackObjects_ParentImageNumber =TrackObjects_ParentImageNumber; %?
    OUT.TrackObjects_ParentObjectNumber =TrackObjects_ParentObjectNumber; %?
    OUT.Intensity_IntegratedIntensity_CorrGreen =Intensity_IntegratedIntensity_CorrGreen; %?
    OUT.Intensity_MeanIntensity_CorrGreen =Intensity_MeanIntensity_CorrGreen; %?
    OUT.Intensity_MeanIntensity_CorrGreenBgSub =Intensity_MeanIntensity_CorrGreenBgSub; 
    OUT.Children = Children; %?
    OUT.Parents = Parents; %?
    OUT.Metadata_Vessel = Metadata_Vessel;
    OUT.Metadata_Site = Metadata_Site;
    OUT.Metadata_Well = Metadata_Well;
    OUT.GROUP = GROUP;
    OUT.GROUP_list = GROUP_list;
    OUT.GROUP_n = GROUP_n;
    OUT.Metadata_Year = Metadata_Year;
    OUT.Metadata_Month = Metadata_Month;
    OUT.Metadata_Day = Metadata_Day;
    OUT.Metadata_Hour = Metadata_Hour;
    OUT.Metadata_Min = Metadata_Min;
    OUT.Metatext_Year =  Metatext_Year;
    OUT.Metatext_Month = Metatext_Month;
    OUT.Metatext_Hour =  Metatext_Hour;
    OUT.Metatext_Day = Metatext_Day;
    OUT.Metatext_Min =  Metatext_Min;
    OUT.Location_Center_X = Location_Center_X;
    OUT.Location_Center_Y = Location_Center_Y;
    OUT.Metatext_TIME = Metatext_TIME;
    OUT.TIME_list = TIME_list;
    OUT.TIME = TIME;
    

    OUT.Obj_label = Obj_label;
    OUT.unique_table = unique_table;
    OUT.Obj_label_life = Obj_label_life;
    OUT.Obj_label_first = Obj_label_first;
    OUT.Obj_label_length = Obj_label_length;
    OUT.Obj_label_last = Obj_label_last;
    OUT.Obj_well = Obj_well;
    OUT.Obj_well_n = Obj_well_n;
    OUT.Obj_GROUP = Obj_GROUP;

    OUT.DataTable_Z = DataTable_Z;
    OUT.DataTable_unstacked = DataTable_unstacked;
    OUT.DataTable_norm = DataTable_norm;
end

function OUT = PTC_construct(DataTable_Z,DataTable_unstacked,bef_window,aft_window,cueT,RepressFlag)
%%RepressFlag - to repress the appearance of the individual subplotting QC
    
    harmfun = @(a,b,c,x) a*cos(x*2*pi./c + b); %a=amp, b=phase(in rad), c=period (in h), x=time (in h)
    tp = DataTable_unstacked{:,1};
    tp_ind_bef = find(tp>=bef_window(1) & tp<bef_window(2));
    tp_ind_aft = find(tp>=aft_window(1) & tp<aft_window(2));
    
    tp_high_res = tp(1):0.1:tp(end);

    fit_bef = NaN(length(tp_high_res),size(DataTable_Z,2));
    fit_aft = NaN(length(tp_high_res),size(DataTable_Z,2));
    Amp_bef = NaN(size(DataTable_Z,2),1);
    Amp_aft = NaN(size(DataTable_Z,2),1);
    Phases_bef = NaN(size(DataTable_Z,2),1);
    Phases_aft = NaN(size(DataTable_Z,2),1);
    Periods_bef = NaN(size(DataTable_Z,2),1);
    Periods_aft = NaN(size(DataTable_Z,2),1);
    fit_peaks_bef = NaN(size(DataTable_Z,2),1);
    fit_peaks_aft = NaN(size(DataTable_Z,2),1);
    fit_peaks_bef_last = NaN(size(DataTable_Z,2),1);
    rsquare_bef = NaN(size(DataTable_Z,2),1);
    rsquare_aft = NaN(size(DataTable_Z,2),1);
    
    period_range = 22.5:0.1:28.5;
    
    k=1;
    
    if ~RepressFlag
        figure;
    end
    
    for ix = 1:size(DataTable_Z,2)
                 
        temp_rsquare_bef = NaN(length(period_range),1);
        temp_Amp_bef     = NaN(length(period_range),1);
        temp_Phase_bef   = NaN(length(period_range),1);
        temp_rsquare_aft = NaN(length(period_range),1);
        temp_Amp_aft     = NaN(length(period_range),1);
        temp_Phase_aft   = NaN(length(period_range),1);
        
        for i = 1:length(period_range)

           h_bef = harmfit(tp(tp_ind_bef)*2*pi/period_range(i),DataTable_Z(tp_ind_bef,ix),1);
           temp_rsquare_bef(i) = rsquared(DataTable_Z(tp_ind_bef,ix),harmfun(h_bef(2),h_bef(3),period_range(i),tp(tp_ind_bef))); 
           temp_Amp_bef(i) = h_bef(2);  
           temp_Phase_bef(i) = h_bef(3);
           
           h_aft = harmfit(tp(tp_ind_aft)*2*pi/period_range(i),DataTable_Z(tp_ind_aft,ix),1);
           temp_rsquare_aft(i) = rsquared(DataTable_Z(tp_ind_aft,ix),harmfun(h_aft(2),h_aft(3),period_range(i),tp(tp_ind_aft))); 
           temp_Amp_aft(i) = h_aft(2);  
           temp_Phase_aft(i) = h_aft(3);

        end
        
        [~,max_rsq_bef] = max(temp_rsquare_bef);
        [~,max_rsq_aft] = max(temp_rsquare_aft);
          
        rsquare_bef(ix) = temp_rsquare_bef(max_rsq_bef); 
        rsquare_aft(ix) = temp_rsquare_aft(max_rsq_aft);

        if rsquare_bef(ix) < 0.5 || rsquare_aft(ix) < 0.5 %|| abs(Periods_bef(ix) - Periods_aft(ix))>0.5
            continue
        end
   
        Amp_bef(ix)     = temp_Amp_bef(max_rsq_bef);
        Amp_aft(ix)     = temp_Amp_aft(max_rsq_aft);
        Phases_bef(ix)  = temp_Phase_bef(max_rsq_bef);
        Phases_aft(ix)  = temp_Phase_aft(max_rsq_aft);
        Periods_bef(ix) = period_range(max_rsq_bef);
        Periods_aft(ix) = period_range(max_rsq_aft);
    
        fit_bef(1:length(tp_high_res),ix) = harmfun(Amp_bef(ix),Phases_bef(ix),Periods_bef(ix),tp_high_res); 
        fit_aft(1:length(tp_high_res),ix) = harmfun(Amp_aft(ix),Phases_aft(ix),Periods_aft(ix),tp_high_res); 

        [~,temppeaks_bef] = findpeaks(fit_bef(tp_high_res >=aft_window(1),ix),10); 
        [~,temppeaks_aft] = findpeaks(fit_aft(tp_high_res >=aft_window(1),ix),10); 
        
        if length(temppeaks_bef) <1 || length(temppeaks_aft) <1
            continue
        end
        fit_peaks_bef(ix) = temppeaks_bef(1);
        fit_peaks_aft(ix) = temppeaks_aft(1);

        %for TIPA
        [~,temppeaks_bef_last] = findpeaks(fit_bef(tp_high_res <=cueT,ix),10); 
        fit_peaks_bef_last(ix) = temppeaks_bef_last(end);

        if k > 100 || RepressFlag
            continue
        end
        subplot(10,10,k);
        scatter(tp,DataTable_Z(:,ix),4,'k','filled');
        hold on
        plot(tp_high_res,fit_bef(:,ix));
        plot(tp_high_res,fit_aft(:,ix));
        title([num2str(Periods_bef(ix),3),' ',num2str(Periods_aft(ix),3)],'FontSize',6)
        set(gca,'XTick',0:24:max(tp),'XLim',[0 max(tp)],'FontSize',6);
        k=k+1;      
      
    end
    
    
    %TIPA
    
    frac_aft = 1 - (cueT - fit_peaks_bef_last)./Periods_bef;
    delt_t_aft = frac_aft.*Periods_aft;
    DELTA = aft_window(1) - (cueT + delt_t_aft); %distance from first peak of null after, to phase determination window
    TIPA_null = (Periods_aft - mod(DELTA,Periods_aft)+1); %TIPA null should be used as an alternative to 'fit_peaks_bef'
    
    TIPA_null_scaled = TIPA_null./Periods_aft;
    fit_peaks_aft_scaled = fit_peaks_aft./Periods_aft;
    
    ph_sh_scaled = phaseShift(fit_peaks_bef,fit_peaks_aft,Periods_bef,Periods_aft,1);
    ph_sh_TIPA = phaseShift(TIPA_null,fit_peaks_aft,Periods_aft,Periods_aft,0);
    ph_sh_TIPA_scaled = phaseShift(TIPA_null,fit_peaks_aft,Periods_aft,Periods_aft,1);
    
   

    %%%%OUTPUT%%%%
    OUT.fit_bef = fit_bef;
    OUT.fit_aft = fit_aft;
    OUT.Amp_bef = Amp_bef;
    OUT.Amp_aft = Amp_aft;
    OUT.Phases_bef = Phases_bef;
    OUT.Phases_aft = Phases_aft;
    OUT.Periods_bef = Periods_bef;
    OUT.Periods_aft = Periods_aft;
    OUT.fit_peaks_bef =fit_peaks_bef;
    OUT.fit_peaks_aft =fit_peaks_aft;
    OUT.rsquare_bef =rsquare_bef;
    OUT.rsquare_aft = rsquare_aft;
    OUT.fit_peaks_bef_last = fit_peaks_bef_last;
    OUT.TIPA_null = TIPA_null;
    OUT.ph_sh_scaled = ph_sh_scaled;
    OUT.ph_sh_TIPA = ph_sh_TIPA;
    OUT.ph_sh_TIPA_scaled = ph_sh_TIPA_scaled;
    OUT.TIPA_null_scaled = TIPA_null_scaled;
    OUT.fit_peaks_aft_scaled = fit_peaks_aft_scaled;

end

function plotDoublePTC(pre_peak,post_peak,PER,colvar,TITLE,NormFlag)
    %PER: the period to plot the PTC according to
    %colvar: a variable for the color axis. to have a constant color, put a
    %single color value.
    %NormFlag: if true, data will be normalized to period length
       
    shaded = [0.7 0.7 0.7];
    sz = 15;
   
    if NormFlag
       pre_peak = pre_peak ./ PER;
       post_peak = post_peak ./ PER;
       PER = PER ./ PER;
    end
    
    hold on;
    
    line([0 0],[-PER PER*2],'Color',[0.5 0.5 0.5]);
    line([-PER PER*2],[0 0],'Color',[0.5 0.5 0.5]);
    line([PER PER],[-PER PER*2],'Color',[0.5 0.5 0.5]);
    line([-PER PER*2],[PER PER],'Color',[0.5 0.5 0.5]);
    line([-PER PER],[0 PER*2],'Color',[0.5 0.5 0.5]);
    line([0 PER*2],[-PER PER],'Color',[0.5 0.5 0.5]);

    scatter(pre_peak-PER,post_peak,sz,shaded,'filled');
    scatter(pre_peak+PER,post_peak,sz,shaded,'filled');
    scatter(pre_peak,post_peak-PER,sz,shaded,'filled');
    scatter(pre_peak,post_peak+PER,sz,shaded,'filled');
    scatter(pre_peak-PER,post_peak-PER,sz,shaded,'filled');
    scatter(pre_peak+PER,post_peak+PER,sz,shaded,'filled');
    scatter(pre_peak+PER,post_peak-PER,sz,shaded,'filled');
    scatter(pre_peak-PER,post_peak+PER,sz,shaded,'filled');

    scatter(pre_peak,post_peak,sz,colvar,'filled'); 

    ylim([-PER*0.5,+PER*1.5]);
    xlim([-PER*0.5,+PER*1.5]);
    set(gca,'XTick',-PER*0.5:PER*0.25:PER*1.5,'XGrid','off','YTick',-PER*0.5:PER*0.25:PER*1.5, ...
        'XTickLabel',{' .5','.75',' 0','.25',' .5','.75',' 0','.25',' .5'},...%[PER*0.5,PER*0.75,0,PER*0.25,PER*0.5,PER*0.75,0,PER*0.25,PER*0.5],...
        'YTickLabel',{' .5','.75',' 0','.25',' .5','.75',' 0','.25',' .5'},...%[PER*0.5,PER*0.75,0,PER*0.25,PER*0.5,PER*0.75,0,PER*0.25,PER*0.5],...
        'TickLength',[0.03 0.03],'FontSize',7,'FontName','Arial',...
        'DataAspectRatio',[1 1 1],'XTickLabelRotation',90);

    title(TITLE,'FontSize',8,'FontName','Arial');
    xlabel('Old Phase');
    ylabel('New Phase');
    box on
    line([-PER PER*2],[-PER PER*2],'Color',[0.3 0.3 0.3]);
    
end

function ph_sh = phaseShift(pre_peak,post_peak,pre_period,post_period,normFlag)
%The function convert PTC data to PRC data
    if normFlag %scale to period
        ph_sh = (post_peak./post_period - pre_peak./pre_period);
        ph_sh(ph_sh < -0.5) = ph_sh(ph_sh < -0.5) +1;
        ph_sh(ph_sh > 0.5) = ph_sh(ph_sh > 0.5) -1;
    else %no scaling, notice that pre_period is unused here
        ph_sh = (post_peak - pre_peak);
        ph_sh(ph_sh < -0.5*post_period) = ph_sh(ph_sh < -0.5*post_period) +post_period(ph_sh < -0.5*post_period);
        ph_sh(ph_sh > 0.5*post_period) = ph_sh(ph_sh > 0.5*post_period) -post_period(ph_sh > 0.5*post_period);
    end
end

function post_peak = Shift2phaseAft(pre_peak,ph_sh,post_period,normFlag)
%The function convert PTC data to PRC data

    post_peak = pre_peak + ph_sh;
    post_peak(post_peak < 0) = post_peak(post_peak < 0) +post_period;
    post_peak(post_peak > post_period) = post_peak(post_peak > post_period) - post_period;
    if normFlag %scale to period 
        post_peak = post_peak./post_period; 
    end
end

function plotDoublePRC(pre_peak,ph_sh,PER,colvar,TITLE,NormFlag)
%PER: the period to plot the PRC according to
    %colvar: a variable for the color axis. to have a constant color, put a
    %single color value.
    %NormFlag: if true, data will be normalized to period length
       
    shaded = [0.7 0.7 0.7];
    sz = 15;
   
    if NormFlag
       pre_peak = pre_peak ./ PER;
       ph_sh = ph_sh ./ PER;
       PER = PER ./ PER;
    end
    
    hold on;
    line([0 0],[-0.5*PER PER*0.5],'Color',[0.5 0.5 0.5]);
    line([PER PER],[-0.5*PER PER*0.5],'Color',[0.5 0.5 0.5]);
    
    scatter(pre_peak-PER,ph_sh,sz,shaded,'filled');
    scatter(pre_peak+PER,ph_sh,sz,shaded,'filled');
    scatter(pre_peak,ph_sh,sz,colvar,'filled'); 

    ylim([-PER*0.5,+PER*0.5]);
    xlim([-PER*0.5,+PER*1.5]);
    set(gca,'XTick',-PER*0.5:PER*0.25:PER*1.5,'XGrid','off','YTick',-PER*0.5:PER*0.25:PER*0.5, ...
        'XTickLabel',{' .5','.75',' 0','.25',' .5','.75',' 0','.25',' .5'},...%[PER*0.5,PER*0.75,0,PER*0.25,PER*0.5,PER*0.75,0,PER*0.25,PER*0.5],...
        'YTickLabel',{'-.5','-.25','  0',' .25',' .5'},...
        'TickLength',[0.03 0.03],'FontSize',7,'FontName','Arial',...
        'DataAspectRatio',[1 1 1],'XTickLabelRotation',90);

    title(TITLE);
    xlabel('Old Phase');
    ylabel('Phase Shift');
    box on
    line([-0.5*PER PER*1.5],[0 0],'Color',[0.3 0.3 0.3]);
    
end

function maxHeight = fitfour1(phBef, phSh)
        [fitres_type1,~] = fit(phBef, phSh,'fourier2','Upper',[Inf,Inf,Inf,Inf,Inf,2*pi]);
        %[fitres_type1,gof1] = fit(ph_bef, ph_sh,cosfun);
        [~,Y1] = predint(fitres_type1,phBef,0.95,'observation','off');
        maxHeight = max(abs(Y1));       
end


function OUT = boostrap_rmse(PTCStruct,well_labels,r_th,per_low,per_high,n_iter,subrow,subcol)
        rmse0 = NaN(length(PTCStruct),1);
        rmse1 = NaN(length(PTCStruct),1);

        mypink = [208,28,139]./255;
        mygreen = [77,172,38]./255;

        clear('rmse_dlt_table');  
tic 
        
     figure;   
    for i = 1:length(PTCStruct)
        conds = find(PTCStruct(i).rsquare_bef > r_th & PTCStruct(i).rsquare_aft > r_th & PTCStruct(i).Periods_aft > per_low & PTCStruct(i).Periods_aft < per_high);

        ph_bef1 = [PTCStruct(i).TIPA_null_scaled(conds)-1; PTCStruct(i).TIPA_null_scaled(conds); PTCStruct(i).TIPA_null_scaled(conds)+1];
        ph_sh1 = [PTCStruct(i).ph_sh_TIPA_scaled(conds);PTCStruct(i).ph_sh_TIPA_scaled(conds);PTCStruct(i).ph_sh_TIPA_scaled(conds)];
        ph_aft1 = [PTCStruct(i).fit_peaks_aft_scaled(conds);PTCStruct(i).fit_peaks_aft_scaled(conds);PTCStruct(i).fit_peaks_aft_scaled(conds)];

        ph_bef1 = ph_bef1(~isnan(ph_bef1) & ~isnan(ph_sh1));
        ph_sh1 = ph_sh1(~isnan(ph_bef1) & ~isnan(ph_sh1));
        ph_aft1 = ph_aft1(~isnan(ph_bef1) & ~isnan(ph_sh1));
        [ph_bef, IX] = sort(ph_bef1);
        ph_sh = ph_sh1(IX);
        ph_aft = ph_aft1(IX);

        [~,gof1] = fit(ph_bef, ph_sh,'fourier2','Upper',[Inf,Inf,Inf,Inf,Inf,2*pi]);

        mean_ph_aft = circ_mean(ph_aft*2*pi)/(2*pi);
        ph_aft_adj = mod(ph_aft - mean_ph_aft + 0.5,1);
        [~,gof0] = fit(ph_bef, ph_aft_adj,'fourier2','Upper',[Inf,Inf,Inf,Inf,Inf,2*pi]);

        rmse0(i) = gof0.rmse;
        rmse1(i) = gof1.rmse;

       rmse0_smp = NaN(n_iter,1);
       rmse1_smp = NaN(n_iter,1);

        for j = 1:n_iter
            y = datasample([PTCStruct(i).TIPA_null_scaled(conds),PTCStruct(i).ph_sh_TIPA_scaled(conds),PTCStruct(i).fit_peaks_aft_scaled(conds)],length(conds));
            ph_bef1 = [y(:,1)-1; y(:,1); y(:,1)+1];
            ph_sh1 = [y(:,2);y(:,2);y(:,2)];
            ph_aft1 = [y(:,3);y(:,3);y(:,3)];

            ph_bef1 = ph_bef1(~isnan(ph_bef1) & ~isnan(ph_sh1));
            ph_sh1 = ph_sh1(~isnan(ph_bef1) & ~isnan(ph_sh1));
            ph_aft1 = ph_aft1(~isnan(ph_bef1) & ~isnan(ph_sh1));
            [ph_bef, IX] = sort(ph_bef1);
            ph_sh = ph_sh1(IX);
            ph_aft = ph_aft1(IX);

            [~,gof1] = fit(ph_bef, ph_sh,'fourier2','Upper',[Inf,Inf,Inf,Inf,Inf,2*pi]);

            mean_ph_aft = circ_mean(ph_aft*2*pi)/(2*pi);
            ph_aft_adj = mod(ph_aft - mean_ph_aft + 0.5,1);
            [~,gof0] = fit(ph_bef, ph_aft_adj,'fourier2','Upper',[Inf,Inf,Inf,Inf,Inf,2*pi]);

            rmse0_smp(j) = gof0.rmse;
            rmse1_smp(j) = gof1.rmse;        

        end

        rmse_dlt_table{i} = (rmse0_smp - rmse1_smp);

        RMSE_p0(i) =  length(find(rmse0_smp > rmse1_smp))./n_iter; %propaply this is the one to use
        RMSE_p1(i) =  length(find(rmse0_smp < rmse1_smp))./n_iter; % this will give the complementary value. they cannot both be significant, but there is possibility that neither will be.
        RMSE_p_diff(i) = length(find((rmse0_smp - rmse1_smp) >= (rmse0(i) - rmse1(i))))./n_iter;

    subplot(subrow,subcol,i);


    rmse_dlt = rmse_dlt_table{i};
    rmse_dlt_neg = rmse_dlt(rmse_dlt <0);
    rmse_dlt_pos = rmse_dlt(rmse_dlt >=0);


    histogram(rmse_dlt_neg,'FaceColor',mygreen, 'BinWidth', 0.005,'EdgeColor', 'none');
    hold on
    histogram(rmse_dlt_pos,'FaceColor',mypink, 'BinWidth', 0.005,'EdgeColor', 'none');
    XLIM = get(gca,'XLim');
    XLIM = [-max(abs(XLIM)), max(abs(XLIM))];
    set(gca,'XLim',XLIM,'FontSize',7,'FontName','Arial');
    title([well_labels{i},'; p_0=',num2str(RMSE_p0(i)),'; p_1=',num2str(RMSE_p1(i))],'FontSize',7,'FontName','Arial');
    xlabel('RMSE_{type-0} - RMSE_{type-1}');
    box off
    h = gca; h.YAxis.Visible = 'off'; 
    xline(0,'Color','k','LineWidth',1);

    end

toc

OUT.rmse_dlt_table = rmse_dlt_table;
OUT.RMSE_p0 = RMSE_p0;
OUT.RMSE_p1 = RMSE_p1;
OUT.per_low = per_low;
OUT.per_high = per_high;
OUT.r_th = r_th;
OUT.n_iter = n_iter;
OUT.well_labels = well_labels;

end
