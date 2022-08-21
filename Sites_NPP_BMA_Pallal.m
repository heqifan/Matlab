close all; clear all; clc;

styear = 2003; edyear = 2017;
wks = 'G:\BMA_text';
parfor year = styear : edyear
    fprintf('value of year: %d\n', year);
    %% Read models
    %LAI
    ftif = [wks, '\Normal_LAI\','Normal_LAI_', num2str(year), '.tif'];    
    [npp_lai, R, cc] = geotiffread(ftif);  
    npp_lai = npp_lai(:);
    % Geodata
    ftif = [wks, '\Normal_Geodata\','Normal_Geodata_', num2str(year), '.tif'];    
    [npp_geo, R, cc] = geotiffread(ftif);  
    info = geotiffinfo(ftif);   
    npp_geo = npp_geo(:);
    % Glass
    ftif = [wks, '\Normal_GLASS\','Normal_GLASS_', num2str(year), '.tif'];    
    [npp_glass, R, cc] = geotiffread(ftif);  
    npp_glass = npp_glass(:);
    %MODIS,
    ftif = [wks, '\Normal_MODIS\','Normal_MODIS_', num2str(year), '.tif'];    
    [npp_modis, R, cc] = geotiffread(ftif);  
    npp_modis = npp_modis(:);
    %TPDC
    ftif = [wks, '\Normal_TPDC\','Normal_TPDC_', num2str(year), '.tif'];    
    [npp_tpdc, R, cc] = geotiffread(ftif);  
    npp_tpdc = npp_tpdc(:);
    %W
    ftif = [wks, '\Normal_W\','Normal_W_', num2str(year), '.tif'];    
    [npp_w, R, cc] = geotiffread(ftif);  
    npp_w = npp_modis(:);
    
    
    dat = [npp_lai,npp_geo, npp_glass, npp_modis, npp_tpdc, npp_w];
    m = 5;  
    b = dat(all(~isnan(dat),2),:); %É¾³ýº¬ÓÐNANµÄÐÐ

    %% BMA
    Tainingdata =b;
    [wk0, sigmak0, loglikelihood0] = codeBMA(Tainingdata);
    figure;
    wk_model = [];
    for i = 1 : m
        subplot(2,3,i);
        hist(wk0(:,i));
        w = quantile(wk0(:,i), [0.5, 0.95, 0.99]); 
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
        wk_model(i,:) = [wm, ws, w];
    end
    disp(wk_model);
    % mean, std, median, x_0.95, x_0.99
    csvpath = ['..\wk_model_',  num2str(year), '.csv'];  
    csvwrite(csvpath, wk_model);
    disp('program is over') ;
end
delete(p);

