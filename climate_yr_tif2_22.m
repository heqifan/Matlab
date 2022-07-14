clear all; clc; close all;
setenv('BLAS_VERSION','')  %定义环境变量名  变量名为“BLAS_VERSION”  值为空，等于setenv('BLAS_VERSION')
wks = 'F:/Cheng';
owk = 'F:/tongjiresult';  %定义输出路径

if ~exist(owk, 'dir')  %判断输出路径是否存在，不存在就创建
    mkdir(owk);
end

scenario = {'RCP45','RCP85'};   %定义输入文件夹名
scenario2 = {'Result_RCP45','Result_RCP85'};   %定义对应输出文件夹名

vflx = {'prc','tas', 'hum','clo'};    %定义变量名
prc=[];tas=[];hum=[];clo=[];
vflxxx={prc,tas,hum,clo};
vflxx = {'降水','气温', '相对湿度','云量'};    %对应的中文名，用于设置图例

uflx = {'mm', '^oC', '%', '%'};   %定义单位
exp_num = 5;   %定义分区的数量
yr1 = 2006; yr2 = 2099;   %定义开始年份和结束年份

ftif = 'F:/Climate4region/Climate_4R_CEVSA.tif';   %分区的tif数据
info = geotiffinfo(ftif);   %读取分区tif数据的基本信息
[msk, R, bbox] = geotiffread(ftif);  %msk为二维数组，R为空间参考,包括数据的大小，经纬度范围与起始经纬度，数据单位、方向等），bbox为边界信息

dat = load([wks, filesep, '2005']);  %读取数据
carea = dat(:,1);   %读取数据的第一列：好像是经度
dat = load([wks, '/locat']);  %读取数据
[n, m] = size(dat);  %获取数组的行数和列数
%carea2 = dat(:,1);  %读取数据的第一列：纬度
ll = dat(:,1:2);  %读取数组的一二两列，纬度和经度
lat = round((53.5-ll(:,1))*10+1);  %计算中国的纬度
lon = round((ll(:,2)-73.3)*10+1);   %计算中国的经度
np = max(lon); nl = max(lat);   %计算最大的行数和列数
disp([np nl]);
dat = [];
xsum = [];

carea_img= -9999*ones([nl np]);  %定义一个全是-9999的数组，数组大小为行数为  nl:354 ，列数为np:617

for k = 1 : n
    carea_img(lat(k),lon(k)) = carea(k);   %这一步好像是赋值          有疑惑？？？？？？？？
end


Trend_rcp45 = []; Rsq_rcp45 = [];  pvalue_rcp45 = [];    %定义多个空的要存放数据的数组
Trend_rcp85 = []; Rsq_rcp85 = [];  pvalue_rcp85 = [];

%p=parpool(4);
%prc = -9999*ones([nl np]); tas = -9999*ones([nl np]); hum = -9999*ones([nl np]); clo = -9999*ones([nl np]);   %定义大小一样的各个变量的数组
outdata = -9999*ones([nl np]);
XM = [];
for s = 1 : 2    %循环  {'RCP45','RCP85'}
    xm(:, 1) = yr1:yr2;
    xs(:, 1) = yr1:yr2;
    xc(:, 1) = yr1:yr2;
    
    for yr = yr1 : yr2
        for v = 1 : 4   %循环四个气象因子
            ff = [wks  '/' scenario{s} '/' vflx{v} '/' num2str(yr)];   %获取RCP45或RCP85下的2006-2009年的每个气象因子数据
            disp(ff);
            dat = load(ff);   %导入对应的数据
            [nd, md] = size(dat);    %获取数据(数组)的大小，nd为行，md为列
            xd = dat(:, 3:end);    %获取数据从第三行开始的数据
            
            for k = 1 : n
                if v == 1       %降水是累加值，其他的气温都是求平均值
                    data = sum(xd, 2);     %对行求和
                    outdata(lat(k),lon(k)) = data(k);
                else
                    data = mean(xd, 2);   %对行求平均
                    outdata(lat(k),lon(k)) = data(k);
                end
            end
            %outf = [owk '/' scenario2{s} '/' vflx{v} num2str(yr), '.tif'];
            %geotiffwrite(outf,single(outdata), R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
            xdat =outdata;
            for i = 1 : 5
                if i < 5
                    x1 = xdat(msk == i & xdat > -900);    %xdat是218418个像元，去除中国区域外面的像元，中国面积内的只有95483个?
                    sub_area = carea_img(msk == i & xdat > -900);    %%%msk是气候分区，carea_imag是像元面积的空间数据
                else
                    x1 = xdat(msk < 100 & xdat > -900);    %非中国区域是255，求全国的是1.2.3.4加在一起
                    sub_area = carea_img(msk < 100 & xdat > -900);
                end
                xm(yr - yr1 + 1, i+1) = sum(sum(x1 .* sub_area))/ sum(sub_area);%%%gc m-2a-1
                %xs(yr - yr1 + 1, i+1) = std(x1(:));
                
            end%fenqujieshuxunhuan
            %             if v==1
            %                 xm1(yr - yr1 + 1,2:6)=xm(yr - yr1 + 1,2:6);
            %                 xs1(yr - yr1 + 1, 2:6) = xs(yr - yr1 + 1,2:6);
            %             elseif v==2
            %                 xm2(yr - yr1 + 1,2:6)=xm(yr - yr1 + 1,2:6);
            %                 xs2(yr - yr1 + 1, 2:6) = xs(yr - yr1 + 1,2:6);
            %             elseif v==3
            %                 xm3(yr - yr1 + 1,2:6)=xm(yr - yr1 + 1,2:6);
            %                 xs3(yr - yr1 + 1, 2:6) =xs(yr - yr1 + 1,2:6);
            %             else
            %                 xm4(yr - yr1 + 1,2:6)=xm(yr - yr1 + 1,2:6);
            %                 xs4(yr - yr1 + 1, 2:6) = xs(yr - yr1 + 1,2:6);
            %             end
            
            
        end%%%fenyinzijieshuxunhuan
        %xm(yr - yr1 + 1, i+1) = sum(sum(x1 .* sub_area))/ sum(sub_area);%%%gc m-2a-1
        %xs(yr - yr1 + 1, i+1) = std(x1(:));
        %xc(:,2:6) = xs(:,2:6)./xm(:, 2:6).*100;%cv
        %         disp(yr);
        %         if yr==2099
        %             prc=[xm1,xs1,xs1./xm1];
        %            tas=[xm2,xs2,xs2./xm2];
        %             hum=[xm3,xs3,xs3./xm3];
        %             clo=[xm4,xs4,xs4./xm4];
        %         else
        %         end
    end%nianzhishuxunhuan
    
    %delete(p);
    hdr = {'Year','QT', 'SC', 'NC', 'NW','Whole','Year','SD_QT', 'SD_SC', 'SD_NC', 'SD_NW','SD_Whole','Year','CV_QT_%', 'CV_SC_%', 'CV_NC_%', 'CV_NW_%','CV_Whole_%'};
    hdrr={'青藏高原区','热带-亚热带季风区','温带季风区','温带大陆区','全国陆地'};
    for j=1:4
        xlswrite([owk, filesep, scenario{s}, '_', vflx{j}, '.xls'], vflxxx{j});%均值和标准差.cv共同写出%%%
    end
    
    %
    XM = [XM, xm];%两个情景的均值?
    x0 = XM(:, 1);  % year
    X0 = [ones(length(x0), 1), x0];
    
    for r = 1 : 5
        y1 = XM(:, r + 1);
        y2 = XM(:, r + 8);
        
        disp([vflx{v}, ' ', hdr{r+1}, ': RCP45 = ', num2str(mean(y1)), ' (', num2str(std(y1)), '), RCP85 = ', num2str(mean(y2)), ' (', num2str(std(y2)),')',]);
        
        [b1, bs1, re1, rs1, p1] = regress(y1, X0);
        [b2, bs2, re2, rs2, p2] = regress(y2, X0);
        
        yp1 = b1(1) + b1(2) * x0;
        yp2 = b2(1) + b2(2) * x0;
        
        s1 = ['y = ', num2str(b1(2), '%.4f') ' x + ', num2str(b1(1), '%.2f')];
        r1 = ['R^2 = ', num2str(p1(1), '%.2f'), ' (p = ', num2str(p1(3), '%.4f'), ')'];
        
        s2 = ['y = ', num2str(b2(2), '%.4f') ' x + ', num2str(b2(1), '%.2f')];
        r2 = ['R^2 = ', num2str(p2(1), '%.2f'), ' (p = ', num2str(p2(3), '%.4f'), ')'];
        
        figure(r);
        plot(x0, y1, 'b.', x0, y2, 'r.', 'MarkerSize', 13);
        hold on;
        plot(x0, yp1, 'b-', x0, yp2, 'r-', 'Linewidth', 2);
        hold off
        
        xlabel('\fontname{宋体}年?'); ylabel([vflxx{v}, ' (', uflx{v}, ')']);
        
        legend('\fontname{Times New Roman}RCP45','\fontname{Times New Roman}RCP85');
        legend('boxoff');
        title(hdrr{r});
        set(gca, 'FontSize', 14);
        fig = [owk, '/' vflx{v} '/' scenario{s},'_', vflx{v}, '_', hdr{r + 1} '.jpg'];
        saveas(gca, fig);
        
        Trend_rcp45(r, v) = b1(2);  Rsq_rcp45(r, v) = p1(1);  pvalue_rcp45(r, v) = p1(3);      %趋势 R2 P值?
        Trend_rcp85(r, v) = b2(2);  Rsq_rcp85(r, v) = p2(1);  pvalue_rcp85(r, v) = p2(3);
    end
    
end

T_45 = [Trend_rcp45, Rsq_rcp45, pvalue_rcp45];
T_85 = [Trend_rcp85, Rsq_rcp85, pvalue_rcp85];

fxls = [owk, '/'  scenario{1},'_meteo',  '_trend_r2_p0711', '.xls'];
xlswrite(fxls, T_45);

fxls = [owk, '/'  scenario{2},'_meteo',  '_trend_r2_p0711', '.xls'];
xlswrite(fxls, T_85);
