% Spatial_ET_BMA.m
close all; clear all; clc;

% mean, std, median, x_0.95, x_0.99
wk_model = csvread('wk_model.csv');   %读取csv数据
wk = wk_model(:, 3);  %读取csv数据的第三列

model_name = {'PT-JPL','ARTIS','MOD16','SSEop', 'GLO-AT'} ;    %每种模型的名字

wks = 'D:\Sanjy250m\ET_BMA\Annually';   
yr1 = 2003;
yr2 = 2017;



a = [1 2 3 4 5 6 7 8 9 10;2 2 3 4 5 6 7 8 9 10;3 2 3 4 5 6 7 8 9 10;4 2 3 4 5 6 7 8 9 10;5 2 3 4 5 6 7 8 9 10;6 2 3 4 5 6 7 8 9 10];
bnd = quantile(a, [0, 0.5, 0.5,1]);    %定义百分位
%% Read Annual ET from each model
for yr = yr1 : yr2   %年份循环
    figure(1);   %定义第一个图窗
    % PT-JPL模拟的蒸散，单位mm/m2/a
    ftif = [wks, '\JPL\JPL_ET', num2str(yr), '.tif'];    % PT-JPL数据的路径
    [et_jpl, R, cc] = geotiffread(ftif);  %读取PT-JPL数据
    info = geotiffinfo(ftif);   %获取PT-JPL数据的基本信息
    
    bnd = quantile(et_jpl(et_jpl > -900), [0, 0.05, 0.5, 0.95, 1]);  %返回pt-jpl数据中大于-900的数据的百分位，0，0.05，0.95，1
    disp(bnd)   
    subplot(2,3,1);   %定义画布的位置
    imagesc(et_jpl, [0, 600]);  %展示图像  0-600之间的一个映射
    title('PT-JPL'); colorbar('horizaonal'); axis off; axis equal;      %title：标题，
    %colbar； 在当前坐标区或图的右侧显示一个垂直颜色栏。颜色栏显示当前颜色图并指示数据值到颜色图的映射。
    
    % ARTIS: bi-weekly ET, mm/m2/day,
    ftif = [wks, '\ARTS\ET_YH', num2str(yr), '.tif'];
    [et_artis, rr, cc] = geotiffread(ftif);
    bnd = quantile(et_artis(et_artis > -900), [0, 0.05, 0.5, 0.95, 1]);
    subplot(2,3,2);
    disp(bnd)
    imagesc(et_artis, [0, 600]);
    title('ARTS'); colorbar('horizaonal'); axis off; axis equal;
    
    % MOD16: mm/m2/8days
    ftif = [wks, '\MOD16\MOD_ET', num2str(yr), '.tif'];
    [dat, rr, cc] = geotiffread(ftif);
    et_mod16 = double(dat);    %将数据类型转为双精度
    et_mod16(dat > 60000) = -9999;
    bnd = quantile(et_mod16(et_mod16 > -900), [0, 0.05, 0.5, 0.95, 1]);    
    disp(bnd)
    subplot(2,3,3);
    imagesc(et_mod16, [0, 600]);
    title('MOD16'); colorbar('horizaonal'); axis off; axis equal;
    
    % SSEBop: The operational simplified Surface Energy Balance (SSEBop)
    % Senay et al. 2012
    % monthly ET
    ftif = [wks, '\SSEBop\SSEBop_ET_', num2str(yr), '.tif'];
    [et_sseb, rr, cc] = geotiffread(ftif);
    bnd = quantile(et_sseb(et_sseb > -900), [0, 0.05, 0.5, 0.95, 1]);
    disp(bnd)
    subplot(2,3,4);
    imagesc(et_sseb, [0, 600]);
    title('SSEBop'); colorbar('horizaonal'); axis off; axis equal;
    
    % GLCV: Bi-weekly ET in 8km of resolution
    
    % GLCV: Bi-weekly ET in 250m of resolution, 8days
    ftif = [wks, '\GLO-AT_ET\ET_', num2str(yr), '.flt'];
    fip = fopen(ftif, 'r');     %读取flt数据时候用，先打开
    dat = fread(fip, [5231, 3201], 'float32');    %再读取
    fclose(fip);   %再关闭
    et_glcv = dat';     %转置
    
    % file_out = [wks, '\GLO-AT_ET\GLO-AT_ET_',  num2str(yr), '.tif'];
    % geotiffwrite2(file_out,et_glcv, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

    %[et_glcv, rr, cc] = geotiffread(ftif);
    bnd = quantile(et_glcv(et_glcv > -900), [0, 0.05, 0.5, 0.95, 1]);
    disp(bnd)
    subplot(2,3,5);
    imagesc(et_glcv, [0, 600]);
    title('ARTS_{opt}'); colorbar('horizaonal'); axis off; axis equal;
    
    % ET_BMA
    % model_name = {'PT-JPL','ARTIS','MOD16','SSEop', 'GLO-AT'} ;
    et_bma = wk(1) * et_jpl + wk(2) * et_artis + wk(3) * et_mod16 + wk(4) * et_sseb + wk(5) * et_glcv;  %对应的模型数据乘以对应的权重
    bnd = quantile(et_bma(et_bma > -900), [0, 0.05, 0.5, 0.95, 1]);
    disp(bnd)
    subplot(2,3,6);
    imagesc(et_bma, [0, 600]);
    title('BMA');colorbar('horizaonal'); axis off; axis equal;
    file_png = [wks, '\BMA\Models_ET_',  num2str(yr), '.png'];
    print(file_png, '-r600', '-dpng');

    
    figure(2)
    imagesc(et_bma, [0, 600]); colorbar('horizaonal'); axis off; axis equal; title(['BMA ET ', num2str(yr)]);
    file_png = [wks, '\BMA\BMA_ET_',  num2str(yr), '.png'];
    print(file_png, '-r600', '-dpng');
    file_out = [wks, '\BMA\BMA_ET_',  num2str(yr), '.tif'];
    geotiffwrite2(file_out,et_bma, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);  %写出数据
    
    
    disp(wk')
end

