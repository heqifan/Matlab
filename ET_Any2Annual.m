close all; clear all; clc;
% step = 'week';
% switch step
%     case 'month'
% monthly
yr1 = 2003;
yr2 = 2018;
wks = 'F:\Sanjy250m\ET\SSEBop';

for yr = yr1 : yr2
    for m = 1 : 12
        jd = datenum(yr, m, 1) - datenum(yr, 1, 1) + 1;
        ftif = [wks, '\SSE_ET', num2str(yr * 1000 + jd), '.tif'];
        if exist(ftif, 'file');
            [dat, R, cc] = geotiffread(ftif);
            out_proj = geotiffinfo(ftif);            
            et = double(dat);
            bnd = quantile(et(et > -900), [0, 0.05, 0.5, 0.95, 1]);
            disp([m, bnd])
            % figure(m); imagesc(et, [0, 200]);
            if m == 1
                aet = et;
            else
                aet = aet + et;
            end
        end
    end
    disp(yr)
    bnd = quantile(aet(aet >0), [0, 0.05, 0.5, 0.95, 1]);
    imagesc(aet, [0, 2000]);colorbar
    title(yr);
    fout = ['F:\Sanjy250m\ET\Annually\SSEBop\SSEBop_ET_', num2str(yr), '.tif'];
    geotiffwrite(fout, aet, R, ...
        'GeoKeyDirectoryTag', out_proj.GeoTIFFTags.GeoKeyDirectoryTag);
    
end
%     case 'week'
% weekly
yr1 = 2001;
yr2 = 2018;
wks = 'F:\Sanjy250m\ET\Weekly\MODIS';

for yr = yr1 : yr2
    for m = 1 : 46
        jd = (m - 1) * 8 + 1;
        ftif = [wks, '\MOD_ET', num2str(yr * 1000 + jd), '.tif'];
        if exist(ftif, 'file');
            [dat, R, cc] = geotiffread(ftif);
            out_proj = geotiffinfo(ftif);
            et = double(dat) / 10; et(dat > 3000) = 0;
            bnd = quantile(et(et > -900), [0, 0.05, 0.5, 0.95, 1]);
            disp([m, bnd])
            % figure(m); imagesc(et, [0, 200]);
            if m == 1
                aet = et;
            else
                aet = aet + et;
            end
        end
    end
    disp(yr)
    bnd = quantile(aet(aet >0), [0, 0.05, 0.5, 0.95, 1]);
    disp([m, bnd])
    imagesc(aet, [0, 500]);colorbar
    title(yr);
    
    fout = ['F:\Sanjy250m\ET\Annually\MOD16\MOD16_ET_', num2str(yr), '.tif'];
    geotiffwrite(fout, aet, R, ...
        'GeoKeyDirectoryTag', out_proj.GeoTIFFTags.GeoKeyDirectoryTag);
end
% end
