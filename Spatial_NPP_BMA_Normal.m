% Spatial_NPP_BMA.m
close all; clear all; clc;

% mean, std, median, x_0.95, x_0.99
wk_model = csvread('K:\HeQiFan\BMA\wk_model_.csv');   %读取csv数据
wk = wk_model(:, 3);  %读取csv数据的第三列

model_name = {'Geodata','GLASS','MODIS','TPDC', 'W'} ;    %每种模型的名字

wks = 'K:\HeQiFan\Out';   
yr1 = 2003;
yr2 = 2017;



a = [1 2 3 4 5 6 7 8 9 10;2 2 3 4 5 6 7 8 9 10;3 2 3 4 5 6 7 8 9 10;4 2 3 4 5 6 7 8 9 10;5 2 3 4 5 6 7 8 9 10;6 2 3 4 5 6 7 8 9 10];
bnd = quantile(a, [0, 0.5, 0.5,1]);    %定义百分位
%% Read Annual ET from each model
for yr = yr1 : yr2   %年份循环
    figure(1);   %定义第一个图窗
    % Geodata
    ftif = [wks, '\Normal_Geodata\',  'Mask_Normal_Geodata_',num2str(yr),'.tif'];    
    [npp_geodata, R, cc] = geotiffread(ftif);  
    info = geotiffinfo(ftif);   
    bnd = quantile(npp_geodata(npp_geodata > -900), [0, 0.05, 0.5, 0.95, 1]);  %返回数据中大于-900的数据的百分位，0，0.05，0.95，1
    disp(bnd)   
    subplot(2,3,1);   %定义画布的位置
    imagesc(npp_geodata, [-1, 1]);  %展示图像  0-1600之间的一个映射
    title('MuSyQ'); colorbar('horizaonal'); axis off; axis equal;      %title：标题，
    %colbar； 在当前坐标区或图的右侧显示一个垂直颜色栏。颜色栏显示当前颜色图并指示数据值到颜色图的映射。
    
    % GLASS
    ftif = [wks, '\Normal_GLASS\','Mask_Normal_GLASS_',num2str(yr),'.tif'];
    [npp_glass, rr, cc] = geotiffread(ftif);
    bnd = quantile(npp_glass(npp_glass> -900), [0, 0.05, 0.5, 0.95, 1]);
    subplot(2,3,2);
    disp(bnd)
    imagesc(npp_glass, [-1, 1]);
    title('GLASS'); colorbar('horizaonal'); axis off; axis equal;
    
    % MODIS
    ftif = [wks, '\Normal_MODIS\', 'Mask_Normal_MODIS_',num2str(yr),'.tif'];
    [npp_modis, rr, cc] = geotiffread(ftif);
    bnd = quantile(npp_modis(npp_modis > -900), [0, 0.05, 0.5, 0.95, 1]);    
    disp(bnd)
    subplot(2,3,3);
    imagesc(npp_modis, [-1, 1]);
    title('MODIS'); colorbar('horizaonal'); axis off; axis equal;
    
    %TPDC
    ftif = [wks, '\Normal_TPDC\','Mask_Normal_TPDC_',num2str(yr),'.tif'];
    [npp_tpdc, rr, cc] = geotiffread(ftif);
    bnd = quantile(npp_tpdc(npp_tpdc > -900), [0, 0.05, 0.5, 0.95, 1]);
    disp(bnd)
    subplot(2,3,4);
    imagesc(npp_tpdc, [-1, 1]);
    title('CASA'); colorbar('horizaonal'); axis off; axis equal;
    
    % W
    ftif = [wks, '\Normal_W\', 'Mask_Normal_W_',num2str(yr),'.tif'];
    [npp_w, rr, cc] = geotiffread(ftif);
    bnd = quantile(npp_w(npp_w > -900), [0, 0.05, 0.5, 0.95, 1]);
    disp(bnd)
    subplot(2,3,4);
    imagesc(npp_w, [-1, 1]);
    title('GLOPEM-CEVSA'); colorbar('horizaonal'); axis off; axis equal;
    
    % ET_BMA
    % model_name = {'PT-JPL','ARTIS','MOD16','SSEop', 'GLO-AT'} ;
    npp_bma = wk(1) * npp_geodata + wk(2) * npp_glass + wk(3) *npp_modis + wk(4) * npp_tpdc + wk(5) * npp_w;  %对应的模型数据乘以对应的权重
    bnd = quantile(npp_bma(npp_bma > -900), [0, 0.05, 0.5, 0.95, 1]);
    disp(bnd)
    subplot(2,3,6);
    imagesc(npp_bma, [-1, 1]);
    title('BMA');colorbar('horizaonal'); axis off; axis equal;
    file_png = ['K:\HeQiFan\BMA\Normal',  num2str(yr), '.png'];
    print(file_png, '-r600', '-dpng');

    
    figure(2)
    imagesc(npp_bma, [-1, 1]); colorbar('horizaonal'); axis off; axis equal; title(['BMA_NPP ', num2str(yr)]);
    file_png = ['K:\HeQiFan\BMA\Normal',  num2str(yr), '.png'];
    print(file_png, '-r600', '-dpng');
    file_out = ['K:\HeQiFan\BMA\Normal',   num2str(yr), '.tif'];
    geotiffwrite(file_out,npp_bma, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);  %写出数据
    disp(wk')
    disp(num2str(yr)) ;
end

