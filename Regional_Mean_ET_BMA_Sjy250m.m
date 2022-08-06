close all; clear all; clc
yr1 = 2003; yr2 = 2015;
x = (yr1: yr2)';
xt = fix(yr1 + (yr2 - yr1) / 2);

owk= 'D:\Sanjy250m\ET_BMA\Trend';
owk2 = 'D:\Sanjy250m\ET_BMA\Annually\BMA_Park';

yr1 = 1982;
yr2 = 2018;
nyr = yr2 - yr1 + 1;
bk = (yr1+4:yr2-4);

viway = 'D:\Sanjy250m\ET_BMA\Annually';

v = 1;
d = 1;
% rg = {10,'OP'; 11, 'CJ'; 12, 'LC'; 13, 'HH'};
fmsk = 'E:\Sanjy250m\Parameters\Sanjy_Park.tif';
info = geotiffinfo(fmsk);
[reg, R, bbox] = geotiffread(fmsk);
figure; imagesc(reg, [10, 14]);
sid = unique(reg(reg < 100));
nsd = length(sid);
model = {'BMA','ARTS','JPL','GLO-AT_ET','MOD16','SSEBop'};
file_name = {'BMA_ET_', 'ET_YH', 'JPL_ET', 'GLO-AT_ET_', 'MOD_ET', 'SSEBop_ET_'};
VI = {'ET'};

dn = 5231; dm = 3201;
for m = 5 %: length(model)
    mj = model{m};
    for v = 1 : 1
        tic;
        sub = [viway, '\', model{m}];
        mn = [];
        i = 1; X = []; k = 1;
        tic;
        for yr = yr1 : yr2
            ftif = [sub, '\',file_name{m}, num2str(yr), '.tif'];
            if ~exist(ftif, 'file')
                continue;
            end
            disp(ftif)
            [et, R, cc] = geotiffread(ftif);
            info = geotiffinfo(ftif);
            et = double(et);
            et(et > 3200) = -9999;
            et(isnan(et)) = -9999;
            for j = 1 : nsd + 1
                if j > nsd
                    vi = et(et >=0 & et < 30000 & reg < 100);
                    bnd = quantile(vi, [0.05, 0.95]);
                    mn(i,:) = [i, yr, 250, mean(vi), std(vi), max(vi), min(vi), length(vi)];
                else
                    vi = et(et >=0 & et < 30000 & reg == sid(j));
                    bnd = quantile(vi, [0.05, 0.95]);
                    mn(i,:) = [i, yr, sid(j), mean(vi), std(vi), max(vi), min(vi), length(vi)];
                end
                i = i + 1;
            end
            et(et < 0) = 0;
            if k == 1
                aet = et;
            else
                aet = aet + et;
            end
            
            et(reg <11 | reg > 13) = -9999;
            file_out = [owk2, '\',  VI{v}, '_', model{m}, '_Mean_', num2str(yr), '.tif'];
            geotiffwrite2(file_out, et, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
            disp(yr);
            k = k + 1;
        end
        aet = aet / k;
        figure(1)
        imagesc(aet, [0, 600]);
        colorbar('horizonal'); axis off; axis equal;title(['ET_{',  model{m}, '}']);
        
        file_out = [owk, '\',  VI{v}, '_', model{m}, '_Mean_', num2str(yr1),'-', num2str(yr2), '.tif'];
        geotiffwrite2(file_out, aet, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
        
        % clip for park
        aet(reg < 11 | reg > 13) = -9999;
        file_out = [owk2, '\',  VI{v}, '_', model{m}, '_Mean_', num2str(yr1),'-', num2str(yr2), '.tif'];
        geotiffwrite2(file_out, aet, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
        
        fxls = [owk2, '\ET_Mean_', num2str(yr1),'-', num2str(yr2), '.xls'];
        xlswrite(fxls, mn, model{m});
    end
end

