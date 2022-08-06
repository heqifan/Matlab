close all; clear all; clc;

% A Tutorial on the Piecewise Regression Approach Applied to Bedload Transport Data
warning off;

myfun = @(x0,x)line_pw;
% close all; clear all; clc
yr1 = 2003; yr2 = 2015;
x = (yr1: yr2)';
xt = fix(yr1 + (yr2 - yr1) / 2);

owk= 'D:\Sanjy250m\ET_BMA\Trend';

yr1 = 2003;
yr2 = 2015;
nyr = yr2 - yr1 + 1;
bk = (yr1+4:yr2-4);

viway = 'D:\Sanjy250m\ET_BMA\Annually';

out_name = {'bt', 'bt_{se}', 'R^2_{bt}', 'p', 'b1', 'b2', 'tp', 'b1_{se}', 'b2_{se}', 'tp_{se}',...
    'b_{yp_y}',  'R^2_{yp_y}', 'p_{yp_y}'};


options = optimoptions('lsqcurvefit','Display','none');
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
VI = {'ET'};

fxls = [owk, '\ET_Trend_TP_', num2str(yr1),'-', num2str(yr2), '.xls'];
dn = 5231; dm = 3201; 
for m = 1 : 1
    mj = model{m};
    for v = 1 : 1
        tic;
        sub = [viway, '\', model{m}];
        mn = [];
        i = 1; X = [];
        tic;
        for yr = yr1 : yr2

            ftif = [sub, '\',model{m},  '_', VI{v}, '_', num2str(yr), '.tif'];
            
            disp(ftif)
            [et, R, cc] = geotiffread(ftif);
            info = geotiffinfo(ftif);
            
            %[dn,dm] = size(dat);
            
            % vi = et(et >=0 & dat < 30000);
            % bnd = quantile(vi, [0.05, 0.95]);
            % mn(i,:) = [i, yr, m, mean(vi), std(vi), max(vi), min(vi), bnd];
            % figure(1)
            % imagesc(et, [0, 600]);
            % colorbar('horizonal'); axis off; axis equal
            X(i, :) =  reshape(et, 1, dn * dm);
            i = i + 1;
            disp(yr);
        end
        [mp, np] = size(X);
        disp_pix = fix(np / 100+0.5);
        k = 1;
        % XP = [];
        parfor pix = 1 : np
            xt = 2010;
            
            yr = (yr1:yr2)';
            y1 = X(:, pix);
            
            x = yr(y1 < 10000 & y1 > 0);
            y = y1(y1 < 10000 & y1 > 0);
            
            y2 = unique(y1);
            n2 = length(y2);
            
            n = length(y);
            
            if n < 0.8*nyr || n2 < 0.8*nyr
                xc = -9999 * ones(1, 6);
                xd = -9999 * ones(1, 9);
            else
                %% piece wise linear regression
                % B(i,:) = [b1(1),b1(2),b2(1),b2(2),p1(1),p2(1), r];
                if length(x) == length(y)
                    B = piece_wise_lines([x,y],xt);
                    x0 = [B(1),B(2),B(4), xt];
                end
            
                [xp, resnorm,res,exit_flag,outputs] = lsqcurvefit(@line_pw,x0,x,y, [],[],options);
                yp = line_pw(xp,x);
                
                if length(yp) == length(y)
                    [b, bs, r, rs, p] = regress(y, [ones(length(yp),1), yp]);
                    xc = [xp, p(1), p(3)];
                else
                    xc = [xp, -9999, -9999];
                end
                %% Linearly regress in three periods: 00-10, 10-19, 00-19 with two methods
                % period since 2000 to 2010
                xs = x(x >= 2000 & x <= 2010);
                ys = y(x >= 2000 & x <= 2010);
                ns = length(ys);
                if ns >= 6
                    [b1, b1s, r1, r1s, p1] = regress(ys, [ones(ns,1), xs]);
                else
                    b1 = -9999 * ones([1, 2]);
                    p1 = -9999 * ones([1, 4]);
                end
                
                % period since 2010 to 2019
                xs = x(x >= 2010 & x <= 2019);
                ys = y(x >= 2010 & x <= 2019);
                ns = length(ys);
                if ns >= 6
                    [b2, b2s, r2, r2s, p2] = regress(ys, [ones(ns,1), xs]);
                else
                    b2 = -9999 * ones([1, 2]);
                    p2 = -9999 * ones([1, 4]);
                end
                
                % period since 2000 to 2019
                xs = x(x >= 2000 & x <= 2019);
                ys = y(x >= 2000 & x <= 2019);
                ns = length(ys);
                if ns >= 9
                    [b3, b3s, r3, r3s, p3] = regress(ys, [ones(ns,1), xs]);
                else
                    b3 = -9999 * ones([1, 2]);
                    p3 = -9999 * ones([1, 4]);
                end
                xd = [b1(2), p1(1), p1(3), b2(2), p2(1), p2(3), b3(2), p3(1), p3(3)];
            end
            XP(pix,:) = [xc, xd];
        end
        XP(isnan(XP)) = -9999;
        fout = [owk, '\',  VI{v}, '_', model{m}, '_', num2str(yr1),'-', num2str(yr2), '.flt'];
        fop = fopen(fout, 'w');
        fwrite(fop, XP, 'float32');
        fclose(fop);
        
        % bt, bts, rsq, p; b1, b2, tp, b1_se, b2_se, tp_se; b_yp_y,
        % rsq_yp_y, sig_yp_y
        
        out_name = {'a_0', 'b_1', 'b_2', 'TP', 'R^2', 'p', 'trd1', 'Rsq1','sig1', 'trd2', 'Rsq2','sig2','trd0', 'Rsq0','sig0'};
        for j = 1 : length(out_name)
            x1 = reshape(XP(:, j), dm,  dn);
            file_out = [owk, '\',  VI{v}, '_', model{m}, '_', out_name{j}, '_', num2str(yr1),'-', num2str(yr2), '.tif'];
            geotiffwrite2(file_out, x1, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
        end
        
        owk2 = 'D:\Sanjy250m\ET_BMA\Trend_Park';
        for j = 1 : length(out_name)
            x1 = reshape(XP(:, j), dm,  dn);
            for k = 1 : 5
                if k == 5
                    xc = x1(x1 > -9000 & reg >9 & reg < 14);
                else
                    xc = x1(x1 > -9000 & reg == k + 9);
                end
                disp([j, k, (j - 1) * 5 + k]);
                XR((j - 1) * 5 + k, :) = [j, k, mean(xc), std(xc), max(xc), min(xc), median(xc), length(xc)];
            end
            
            x1(reg < 11 | reg > 13) = -9999;
            file_out = [owk2, '\',  VI{v}, '_', model{m}, '_', out_name{j}, '_', num2str(yr1),'-', num2str(yr2), '.tif'];
            geotiffwrite2(file_out, x1, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
        end
        fxls = [owk2, '\trend_sjy-park_et.xls'];
        xlswrite(fxls, XR);
        
        p = reshape(XP(:, 6), dm,  dn);
        bnd = [-150, 150];
        subplot(2,2,1)
        b1 = reshape(XP(:, 2), dm,  dn);
        b1_bnd = quantile(b1(b1 > -9000), [0.05, 0.95]);
        imagesc(b1, b1_bnd);
        colorbar('horizonal');
        colormap([white(3); parula(250)]);
        title(out_name{2});
        
        subplot(2,2,2)
        b2 = reshape(XP(:, 3),  dn, dm);
        b2_bnd = quantile(b2(b2 > -9000), [0.05, 0.95]);
        imagesc(b2, bnd);
        colorbar('horizonal');
        colormap([white(3); parula(250)]);
        title(out_name{3});
        
        subplot(2,2,3)
        b21 = b2 - b1;
        b21(b2< -9000 | b1 < -9000) = -9999;
        b21_bnd = quantile(b21(b21 > -9000), [0.05, 0.95]);
        imagesc(b21, [-300, 300]);
        colorbar('horizonal');
        colormap([white(2); parula(253)]);
        title('b_2 - b_1');
        
        subplot(2,2,4)
        tp = reshape(XP(:, 4),  dn, dm);
        tp(p > 0.05 | tp < 2001 | tp > 2017) = NaN;
        tp_bnd = quantile(tp(tp > -9000), [0.05, 0.95]);
        td = fix(tp / 5) * 5;
        imagesc(tp, tp_bnd);
        colorbar('horizonal');
        colormap([white(1); parula(10);[0.5,0.5,0.5]]);
        title(out_name{4});
        
        
        fig = [owk, '\', VI{v}, '_', model{m}, '_', num2str(yr1),'-', num2str(yr2), '.png'];
        print(fig, '-dpng', '-r600');
        
        sheet = [VI{v}, '_', model{m}];
        figure
        subplot(2,2, 1);
        [x, nx] = hist(tp(:), 2000:2018);
        xp = x*100/sum(x);
        bar(nx, xp);
        xlswrite(fxls, {['TP_' model{m}], ['num_' model{m}]}, sheet, 'A1');
        xlswrite(fxls, [nx', x'], sheet, 'A2');
        
        subplot(2,2, 2);
        [x1, nx] = hist(b1((b1 > -9000 & p <= 0.05)), -150:150);
        xp1 = x1*100/sum(x1);
        [x2, nx] = hist(b2((b2 > -9000 & p <= 0.05)), -150:150);
        xp2 = x2*100/sum(x2);
        plot(nx, xp1, 'r.', nx, xp2,'b.'); ylim([0, 0.8])
        legend('b_1', 'b_2')
        xlswrite(fxls, {['B_' model{m}], ['B1_' model{m}], ['B2_' model{m}]}, sheet, 'C1');
        xlswrite(fxls, [nx', x1', x2'], sheet, 'C2');
        
        subplot(2,2, 3);
        [x3, nx] = hist(b21((b21 > -9000 & p <= 0.05)), -300:300);
        xp3 = x3*100/sum(x3);
        plot(nx, xp3, 'b.');ylim([0, 0.5])
        xlswrite(fxls, {['B_' model{m}], ['B2-B1_' model{m}]}, sheet, 'F1');
        xlswrite(fxls, [nx', x3'], sheet, 'F2');
        
        subplot(2,2, 4);
        [x3, nx] = hist(b21((b21 > -9000 & p <= 0.05)), -10000:100:10000);
        xp3 = x3*100/sum(x3);
        plot(nx, xp3, 'b.'); % ylim([0, 0.5])
        xlswrite(fxls, {['b_' model{m}], ['b2-b1_' model{m}]}, sheet, 'H1');
        xlswrite(fxls, [nx', x3'], sheet, 'H2');
        
        fig = [owk, '\', VI{v}, '_', model{m}, '_', num2str(yr1),'-', num2str(yr2), '_hist.png'];
        print(fig, '-dpng', '-r600');
    end
end
