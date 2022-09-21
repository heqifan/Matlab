% Spatial_NPP_BMA.m
close all; clear all; clc;

% mean, std, median, x_0.95, x_0.99
wk_model = csvread('K:\HeQiFan\BMA\\wk_model_npp.csv');   %读取csv数据
wk = wk_model(:, 3);  %读取csv数据的第三列

% model_name = {'Geodata','GLASS','MODIS','TPDC', 'W'} ;    %每种模型的名字


wks = 'K:\HeQiFan\1Y';
yr1 = 2000;
yr2 = 2017;



a = [1 2 3 4 5 6 7 8 9 10;2 2 3 4 5 6 7 8 9 10;3 2 3 4 5 6 7 8 9 10;4 2 3 4 5 6 7 8 9 10;5 2 3 4 5 6 7 8 9 10;6 2 3 4 5 6 7 8 9 10];
bnd = quantile(a, [0, 0.5, 0.5,1]);    %定义百分位
%% Read Annual ET from each model
for yr = yr1 : yr2   %年份循环
    figure(1);   %定义第一个图窗
    % Geodata
    ftif = [wks, '\MuSyQ_1981_2018_1y_chinese\', num2str(yr),'\Mask_Resample_Mul_',num2str(yr),'.tif'];    
    [npp_geodata, R, cc] = geotiffread(ftif);  
    npp_geodata(npp_geodata <0) = nan;
    info = geotiffinfo(ftif);   
    bnd = quantile(npp_geodata(npp_geodata > -900), [0, 0.05, 0.5, 0.95, 1]);  %返回数据中大于-900的数据的百分位，0，0.05，0.95，1
    disp(bnd)   
    subplot(2,3,1);   %定义画布的位置
    imagesc(npp_geodata, [0, 1600]);  %展示图像  0-1600之间的一个映射
    title('MuSyQ'); colorbar('horizaonal'); axis off; axis equal;      %title：标题，
    %colbar； 在当前坐标区或图的右侧显示一个垂直颜色栏。颜色栏显示当前颜色图并指示数据值到颜色图的映射。
    
    % GLASS
    ftif = [wks, '\GLASS_1982_2018_1y_chinese\', num2str(yr), '\Mask_Mul_',num2str(yr),'.tif'];
    [npp_glass, rr, cc] = geotiffread(ftif);
    npp_glass(npp_glass <0) = nan;
    bnd = quantile(npp_glass(npp_glass> -900), [0, 0.05, 0.5, 0.95, 1]);
    subplot(2,3,2);
    disp(bnd)
    imagesc(npp_glass, [0, 1600]);
    title('GLASS'); colorbar('horizaonal'); axis off; axis equal;
    
    % MODIS
    ftif = [wks, '\MODIS_2000_2017_1y_chinese\', num2str(yr), '\Mask_Mul_',num2str(yr),'.tif'];
    [npp_modis, rr, cc] = geotiffread(ftif);
    npp_modis(npp_modis <0) = nan;
    %npp_modis(npp_modis <0) = nan;
    bnd = quantile(npp_modis(npp_modis > -900), [0, 0.05, 0.5, 0.95, 1]);    
    disp(bnd)
    subplot(2,3,3);
    imagesc(npp_modis, [0, 1700]);
    title('MODIS'); colorbar('horizaonal'); axis off; axis equal;
    
%     %TPDC
%     ftif = [wks, '\TPDC_2000_2017_1y\', num2str(yr), '\Mask_Reproject_Mul_',num2str(yr),'.tif'];
%     [npp_tpdc, rr, cc] = geotiffread(ftif);
%     %npp_tpdc(npp_tpdc >1700) = nan;
%     bnd = quantile(npp_tpdc(npp_tpdc > -900), [0, 0.05, 0.5, 0.95, 1]);
%     disp(bnd)
%     subplot(2,3,4);
%     imagesc(npp_tpdc, [0, 1700]);
%     title('CASA'); colorbar('horizaonal'); axis off; axis equal;
    
    % W
    ftif = [wks, '\GLOPEM-CEVSA_1980_2020_1y_chinese\', num2str(yr), '\RNPP_',num2str(yr),'.flt'];
%     [npp_w, rr, cc] = geotiffread(ftif);
    fip = fopen(ftif, 'r');     %读取flt数据时候用，先打开
    npp_w = (fread(fip, [4998, 4088], 'float32'))';    %再读取
    fclose(fip);   %再关闭
    npp_w(npp_w <0) = nan;
%     npp_w(npp_w ==3200) = nan;
    bnd = quantile(npp_w(npp_w > -900), [0, 0.05, 0.5, 0.95, 1]);
    disp(bnd)
    subplot(2,3,4);
    imagesc(npp_w, [0, 2000]);
    title('GLOPEM-CEVSA'); colorbar('horizaonal'); axis off; axis equal;
    
%     fip = fopen(ftif, 'r');     %读取flt数据时候用，先打开
%     dat = fread(fip, [5231, 3201], 'float32');    %再读取
%     fclose(fip);   %再关闭
    
    % ET_BMA 
    % model_name = {'PT-JPL','ARTIS','MOD16','SSEop', 'GLO-AT'} ;
    npp_bma = wk(1) * npp_modis + wk(2) * npp_geodata + wk(3) *npp_w + wk(4) * npp_glass;  %对应的模型数据乘以对应的权重
    bnd = quantile(npp_bma(npp_bma > -900), [0, 0.05, 0.5, 0.95, 1]);
    npp_bma(npp_bma <0) = nan;
    disp(bnd)
    subplot(2,3,6);
    imagesc(npp_bma, [0, 2000]);
    title('BMA');colorbar('horizaonal'); axis off; axis equal;
    file_png = ['K:\HeQiFan\BMA\\',  num2str(yr), '.png'];
    print(file_png, '-r600', '-dpng');

    
    figure(2)
    imagesc(npp_bma, [0, 2000]); colorbar('horizaonal'); axis off; axis equal; title(['BMA_NPP ', num2str(yr)]);
    file_png = ['K:\HeQiFan\BMA\\BMA',  num2str(yr), '.png'];
    print(file_png, '-r600', '-dpng');
    file_out = ['K:\HeQiFan\BMA\\BMA_',   num2str(yr), '.tif'];
    geotiffwrite(file_out,npp_bma, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);  %写出数据
    disp(wk')
    disp(num2str(yr)) ;
end

