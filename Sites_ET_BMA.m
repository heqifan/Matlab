close all; clear all; clc;

%% Read models
% ET observations on EC tower, mm/m2/8day --> mm/m2/day
fxls = 'D:\Sanjy250m\ET_BMA\dat\hgc_site_8days_20190606.xls';
[dat, txt, raw] = xlsread(fxls, 'ET');
dt_obs = datenum(dat(:,1), 1, 1) + dat(:,2);
yr1_obs = min(year(dt_obs));yr2_obs = max(year(dt_obs));
et_obs = dat(:,3)/8;

% PT-JPL模拟的海北灌丛站蒸散，单位mm/m2/day,时间分辨率8d
fxls = 'D:\Sanjy250m\ET_BMA\dat\站点提取HBGC_ET.xlsx';
[dat, txt, raw] = xlsread(fxls);
dt_jpl = datenum(txt(2:1611));
et_jpl = dat;
% plot(dt_jpl, et_jpl); datetick('x'); title('PT-JPL');

% ARTIS: bi-weekly ET, mm/m2/day,
fxls = 'D:\Sanjy250m\ET_BMA\dat\ET_HB_ARTIS.xls';
[dat, txt, raw] = xlsread(fxls, 'ARTIS2');
jd = dat(:,1) - dat(:,2) * 1000;
dt_artis = datenum(dat(:,2), 1, 1) + jd;
et_artis = dat(:, 3);
yr = year(dt_artis);
yr1 = min(yr); yr2 = max(yr); xt = [];i = 1;
for yt = yr1: yr2, for d = 1 : 8 : 361,xt(i,1) = datenum(yt,1,1)+d-1;i = i + 1;end;end
et = spline(dt_artis, et_artis, xt);
dt_artis = xt;
et_artis = et;

% MOD16: mm/m2/8dayS
fxls = 'D:\Sanjy250m\ET_BMA\dat\ET_HB-month.xls';
[dat, txt, raw] = xlsread(fxls, 'Sheet1', '$b1:$c706');
yr = fix(dat(:,1)/1000);
jd = dat(:,1) - yr * 1000;
dt_mod16 = datenum(yr, 1, 1) + jd - 1;
et_mod16 = dat(:, 2)/8;
yr = year(dt_mod16);
yr1 = min(yr); yr2 = max(yr); xt = [];i = 1; j = 1; yp = [];
for yt = yr1: yr2,
    for d = 1 : 8 : 361,
        t0 = datenum(yt,1,1)+d-1;
        ii = find(dt_mod16==t0);
        if isempty(ii)
            et(i,1) = -999;
        else
            et(i,1) = et_mod16(ii);
        end
        xt(i,1) = t0;
        i = i + 1;
    end
end
for j = 1 : length(xt)
    if et(j) < -900
        j1 = max(1, j - 5);
        j2 = min(length(xt), j + 5);
        
        x1 = xt(j1 :j2);
        y1 = et(j1 :j2);
        
        x2 = x1(y1 > -900);
        y2 = y1(y1 > -900);
        
        if length(x2) > 3
            [b, bs, r, rs, p] = regress(y2, [ones([length(x2),1]), x2]);
            yp(j, 1) = b(1) + b(2) * xt(j);
        else
            yp(j, 1) = -9999;
        end
    else
        yp(j,1) = et(j);
    end
end
dt_mod16 = xt;
et_mod16 = yp;

% SSEBop: The operational simplified Surface Energy Balance (SSEBop)
% Senay et al. 2012
% monthly ET
fxls = 'D:\Sanjy250m\ET_BMA\dat\ET_HB-month.xls';
[dat, ~, raw] = xlsread(fxls, 'SSE');  %, '$A1:$B191'
yr = fix(dat(:,1)/1000);
jd = dat(:,1) - yr * 1000;
dt_sse = datenum(yr, 1, 1) + jd;
et_sse = dat(:, 6);
yr = year(dt_sse);
yr1 = min(yr); yr2 = max(yr); xt = [];i = 1;
for yt = yr1: yr2, for d = 1 : 8 : 361,xt(i,1) = datenum(yt,1,1)+d-1;i = i + 1;end;end
et = spline(dt_sse, et_sse, xt); et(et < 0) = 0;
dt_sse = xt;
et_sse = et;

% GLCV: Bi-weekly ET in 8km of resolution
fxls = 'D:\Sanjy250m\ET_BMA\dat\ET_HB_GC_8km.xls';
[dat, txt, raw] = xlsread(fxls, 'GC_8km');
yr = fix(dat(:,1)/1000);
jd = dat(:,1) - yr * 1000;
dt_glcv8km = datenum(yr, 1, 1) + jd;
et_glcv8km = dat(:, 3);
yr = year(dt_glcv8km);
yr1 = min(yr); yr2 = max(yr); xt = [];i = 1;
for yt = yr1: yr2, for d = 1 : 8 : 361,xt(i,1) = datenum(yt,1,1)+d-1;i = i + 1;end;end
et = spline(dt_glcv8km, et_glcv8km, xt); et(et < 0) = 0;
dt_glcv8km = xt;
et_glcv8km = et;

% GLCV: Bi-weekly ET in 250m of resolution, 8days
fxls = 'D:\Sanjy250m\ET_BMA\dat\ET_HB_GC_250m.xls';
[dat, txt, raw] = xlsread(fxls, 'GC250');
yr = fix(dat(:,1)/1000);
jd = dat(:,1) - yr * 1000;
dt_glcv = datenum(yr, 1, 1) + jd;
et_glcv = dat(:, 3);

yr1_mod = max(year([min(dt_glcv), min(dt_glcv8km), min(dt_sse), min(dt_artis), min(dt_jpl)]));
yr2_mod = min(year([max(dt_glcv), max(dt_glcv8km), max(dt_sse), max(dt_artis), max(dt_jpl)]));

yr1 = max([yr1_mod, yr1_obs]);
yr2 = min([yr2_mod, yr2_obs]);

dt1_obs = dt_obs(year(dt_obs) >= yr1 & year(dt_obs) <= yr2);
et1_obs = et_obs(year(dt_obs) >= yr1 & year(dt_obs) <= yr2);
et1_jpl = et_jpl(year(dt_jpl) >= yr1 & year(dt_jpl) <= yr2);
et1_artis = et_artis(year(dt_artis) >= yr1 & year(dt_artis) <= yr2);
et1_mod16 = et_mod16(year(dt_mod16) >= yr1 & year(dt_mod16) <= yr2);
et1_sse = et_sse(year(dt_sse) >= yr1 & year(dt_sse) <= yr2);
et1_glcv8km = et_glcv8km(year(dt_glcv8km) >= yr1 & year(dt_glcv8km) <= yr2);
et1_glcv = et_glcv(year(dt_glcv) >= yr1 & year(dt_glcv) <= yr2);

plot(dt1_obs, et1_obs, 'x', dt1_obs, et1_jpl, 'b-', dt1_obs, et1_artis, 'r-', dt1_obs, et1_mod16, 'g-');
hold on;
plot(dt1_obs, et1_sse, 'y-');
plot(dt1_obs, et1_glcv8km, 'k-');
plot(dt1_obs, et1_glcv, 'c-');
datetick('x'); legend('Obs','PT-JPL','ARTIS','MOD16','SSE','GLCV8km', 'GLCV250m');
hold off;

% dat = [et1_obs, et1_jpl, et1_artis, et1_mod16, et1_sse, et1_glcv8km, et1_glcv];
dat = [et1_obs, et1_jpl, et1_artis, et1_mod16, et1_sse,  et1_glcv];
yo = dat(:,1);
rr = corr(dat);
n = 368; m = 5;
%% Performance of each model
B = [];
for i = 1 : m
    yp = dat(:,i + 1);
    [b, bs, r, rs, p] = regress(yo, [ones([n,1]), yp]);
    bias = sum((yo - yp))/n;
    rmse = (sum((yo - yp).^2)/n)^0.5;
    re = rmse / mean(yo);
    nse = 1 - sum((yo - yp).^2)/sum((yo - mean(yo)).^2);
    B(i,:) = [b(1), b(2), p(1), p(3), bias, rmse, re, nse, mean(yp), std(yp),mean(yo), std(yo)];
end

%% BMA
Tainingdata =dat;
[wk0, sigmak0, loglikelihood0] = codeBMA(Tainingdata);
figure;
wk_model = [];
for i = 1 : m
    subplot(2,3,i);
    hist(wk0(:,i));
    w = quantile(wk0(:,i), [0.5, 0.95, 0.99]); % disp(w);
    c = hist(wk0(:,i), w);
    hold on;
    plot(w, c, 'rx')
    wm = mean(wk0(:,i));
    ws = std(wk0(:,i));
    wm1 = [wm-ws*2/3, wm-ws/3, wm, wm+ws/3, wm+ws*2/3];
    cm1 = hist(wk0(:,i), wm1);
    plot(wm1, cm1, 'bo');
    plot(wm, cm1(3), 'kd')
    plot(w(1),c(1), 'gx');
    wk1(1,i) = w(1);
    wk_model(i,:) = [wm, ws, w];
end
disp(wk_model);
% mean, std, median, x_0.95, x_0.99
csvwrite('..\Tables\wk_model.csv', wk_model);

loc = [1 50 100];
for j = 1 : 100
    wk = wk0(j,:);
    sigmak = sigmak0(j,:);
    
    n=size(Tainingdata,1);
    
    for i=1:n
        BMAV(i,j)=sum(wk.*Tainingdata(i,2:end));
        temp1=sum(wk.*(Tainingdata(i,2:end)-sum(wk.*Tainingdata(i,2:end))).^2);
        temp2=sum(wk.*sigmak.^2);
        VarV(i,1)=sqrt(temp1+temp2);
    end
    
end

for i = 1 : n
    et_bma(i,1) = sum(wk1 .* Tainingdata(i,2:end));
end
figure; plot(et_bma, et_obs, '.');
et_max = max([et_obs; et_bma]);
hold on;
plot([0, et_max], [0, et_max], 'k:');
[b, bs, r, rs, p] = regress(et_obs, [ones(length(et_bma),1), et_bma]);
plot(et_bma, b(1) + b(2) * et_bma, '-');
xlabel('ET_{BMA}'); ylabel('ET_{OBS}');
title(['y = ', num2str(b(1), '%.2f'), ' + ', num2str(b(2), '%.2f'), ' x ( R^2 = ', num2str(p(1), '%.2f'), ', p = ', num2str(p(3), '%.2f'), ')']);
csvwrite('..\Tables\ET_BMA_Haibei.csv', [dt1_obs, et_bma, et_obs]);

%% dat = [et1_obs, et1_jpl, et1_artis, et1_mod16, et1_sse, et1_glcv8km, et1_glcv];
figure;

num_times = length(et_bma) ; station_years = '2003-2010' ;

station_name = 'HaiBeiShrub' ; data_unit = 'mm·m^(-2)·day^(-1)' ;

% model_name = {'PT-JPL','ARTIS','MOD16','SSEop','GLCV-AVHRR', 'GLCV-MODIS'} ;
model_name = {'PT-JPL','ARTIS','MOD16','SSEop', 'GLO-AT'} ;

num_models = length(model_name);

chinflux_et = dat(:,1);
BMA_et      = et_bma;
all_model_et = (dat(:, 2:end))';
%--------------------------------------------------------------------------
% To combine all kinds of ET in one array:
%  et(1,:) the ET from flux obeservation.
%  et(2,:) the ET from BMA method.
%  et(3:num_models,:) the ET from different models.
%--------------------------------------------------------------------------
all_et = zeros(num_models+2,num_times) ;

all_et(1,:) = chinflux_et' ;
all_et(2,:) = BMA_et' ;
for k = 3 : num_models + 2
    all_et(k,:) =  all_model_et(k-2,:) ;
end

et = all_et ;

%%
%--------------------------------------------------------------------------
% To calculate statistics for Taylor diagram
% https://github.com/PeterRochford/SkillMetricsToolbox/wiki
%--------------------------------------------------------------------------
struct_et(1).name = station_name ;
struct_et(1).years = station_years ;
struct_et(1).data = et(1,:) ;
struct_et(1).unit = data_unit ;

struct_et(2).name = 'BMA' ;
struct_et(2).years = station_years ;
struct_et(2).data = et(2,:) ;
struct_et(2).unit = data_unit ;

for k = 1 : num_models
    
    struct_et(k+2).name = model_name{k} ;
    struct_et(k+2).years = station_years ;
    struct_et(k+2).data = et(k+2,:) ;
    struct_et(k+2).unit = data_unit ;
    
end

for k = 2 : num_models + 2
    
    et_stats(k-1) = taylor_statistics(struct_et(k),struct_et(1),'data') ;
    
end



%--------------------------------------------------------------------------
% To store statistics in arrays.
%--------------------------------------------------------------------------
sdev(1) = et_stats(1).sdev(1) ;
crmsd(1) = et_stats(1).crmsd(1) ;
ccoef(1) = et_stats(1).ccoef(1) ;

for k = 2 : num_models + 2
    
    sdev(k) = et_stats(k-1).sdev(2) ;
    crmsd(k) = et_stats(k-1).crmsd(2) ;
    ccoef(k) = et_stats(k-1).ccoef(2) ;
    
end
figure
taylordiag(sdev,crmsd,ccoef);

% %-------------------------------------------------------------------------
% %    To normalize for Standard deviations (sdev) and Centered Root Mean
% % Square Difference (crmsd) by obsevered standard deviations.
% %-------------------------------------------------------------------------
% sdev(:) =  sdev(:) ./ sdev(1) ;
% crmsd(:) = crmsd(:) ./ sdev(1) ;
% 
% 
% label = {'Non-Dimensional Observation',''} ;
% 
% % colorm = [ [0 245 255];[0 255 0];[238 180 34];[0 134 139];[139 101 8]; ...
% %     [0 0 0];[0 0 128];[176 48 96];[100 149 237];[205 181 205]; ...
% %     [0 0 255];[255 0 0];[0 238 0];[255 127 80]; ...
% %     [238 238 0];[255 0 255];[34 139 34]] / 255 ;
% colorm = 'rgbcmyk';
% sdev2 = sdev(1:2) ;
% crmsd2 = crmsd(1:2) ;
% ccoef2 = ccoef(1:2) ;
% 
% 
% [hp(1), ht(1), ax] = taylor_diagram(sdev2,crmsd2,ccoef2, ...
%     'checkStats','on','styleOBS','-.','colOBS','r', ...
%     'widthOBS',1.0,'markerobs','*','titleOBS','', ...
%     'markerLabel',label,'markerLabelColor','black', ...
%     'markerColor', colorm(1), ...
%     'tickRMS',0.2:0.2:1.0,'tickRMSangle',130, ...
%     'showlabelsRMS','off','styleRMS','--','colRMS', ...
%     'black','widthRMS',1.0,'titleRMS','off', ...
%     'tickSTD',0.1:0.2:1.1,'limSTD',1.1,'colSTD','black', ...
%     'styleSTD','-.','widthSTD',0.2, 'titleSTD','on', ...
%     'tickCOR',[0.2:0.1:0.9 0.95 0.98 0.99], ...
%     'colCOR','black','styleCOR',':','widthCOR',0.2, ...
%     'titleCOR','on') ;
% 
% for k = 1 : num_models
%     disp(k);
%     
%     label{2} = num2str(k) ;
%     ccoef(isnan(ccoef)) = 0;
%     sdev2(2) = sdev(k+2) ;
%     crmsd2(2) = crmsd(k+2) ;
%     ccoef2(2) = ccoef(k+2) ;
%     
%     hp(k+1) =    taylor_diagram(sdev2,crmsd2,ccoef2, ...
%         'checkStats','on','overlay','on')
% %         'markerLabel',label, 'markerLabelColor','black', ...
% %         'markerColor',colorm(k+1)) ;
% end
% 
% data_label{1} = 'BMA' ;
% for k = 1 : num_models
%     data_label{k+1} = model_name{k} ;
% end
% 
% title('HaiBeiShiDi','position',[0.1,1.15])
% gtext('RMSD','position',[0.5,0])
% h_legend = legend(hp,data_label) ;
% 
% set(h_legend,'Position',[0.83 0.143 0.1 0.7], ...
%     'FontName','Times New Roman','FontSize',8) ;

outputfile = ['..\Figures\', station_name,'_TaylorDiagram'] ;

print(outputfile,'-djpeg') ;


%% Performance of each model
BS = [];
dat = all_et';
[~,m] = size(dat);
yo = dat(:, 1);
for i = 1 : m-1
    yp = dat(:,i+1);
    [b, bs, r, rs, p] = regress(yo, [ones([n,1]), yp]);
    bias = sum((yo - yp))/n;
    rmse = (sum((yo - yp).^2)/n)^0.5;
    re = rmse / mean(yo);
    nse = 1 - sum((yo - yp).^2)/sum((yo - mean(yo)).^2);
    BS(i,:) = [b(1), b(2), p(1), p(3), bias, rmse, re, nse, mean(yp), std(yp),mean(yo), std(yo)];
end
csvwrite('..\Tables\ET_allmodel_Haibei.csv', [year(dt1_obs),month(dt1_obs),day(dt1_obs), dat]);

disp('program is over') ;
