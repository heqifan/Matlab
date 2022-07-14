clear all; clc; close all;
setenv('BLAS_VERSION','')  %���廷��������  ������Ϊ��BLAS_VERSION��  ֵΪ�գ�����setenv('BLAS_VERSION')
wks = 'F:/Cheng';
owk = 'F:/tongjiresult';  %�������·��

if ~exist(owk, 'dir')  %�ж����·���Ƿ���ڣ������ھʹ���
    mkdir(owk);
end

scenario = {'RCP45','RCP85'};   %���������ļ�����
scenario2 = {'Result_RCP45','Result_RCP85'};   %�����Ӧ����ļ�����

vflx = {'prc','tas', 'hum','clo'};    %���������
prc=[];tas=[];hum=[];clo=[];
vflxxx={prc,tas,hum,clo};
vflxx = {'��ˮ','����', '���ʪ��','����'};    %��Ӧ������������������ͼ��

uflx = {'mm', '^oC', '%', '%'};   %���嵥λ
exp_num = 5;   %�������������
yr1 = 2006; yr2 = 2099;   %���忪ʼ��ݺͽ������

ftif = 'F:/Climate4region/Climate_4R_CEVSA.tif';   %������tif����
info = geotiffinfo(ftif);   %��ȡ����tif���ݵĻ�����Ϣ
[msk, R, bbox] = geotiffread(ftif);  %mskΪ��ά���飬RΪ�ռ�ο�,�������ݵĴ�С����γ�ȷ�Χ����ʼ��γ�ȣ����ݵ�λ������ȣ���bboxΪ�߽���Ϣ

dat = load([wks, filesep, '2005']);  %��ȡ����
carea = dat(:,1);   %��ȡ���ݵĵ�һ�У������Ǿ���
dat = load([wks, '/locat']);  %��ȡ����
[n, m] = size(dat);  %��ȡ���������������
%carea2 = dat(:,1);  %��ȡ���ݵĵ�һ�У�γ��
ll = dat(:,1:2);  %��ȡ�����һ�����У�γ�Ⱥ;���
lat = round((53.5-ll(:,1))*10+1);  %�����й���γ��
lon = round((ll(:,2)-73.3)*10+1);   %�����й��ľ���
np = max(lon); nl = max(lat);   %������������������
disp([np nl]);
dat = [];
xsum = [];

carea_img= -9999*ones([nl np]);  %����һ��ȫ��-9999�����飬�����СΪ����Ϊ  nl:354 ������Ϊnp:617

for k = 1 : n
    carea_img(lat(k),lon(k)) = carea(k);   %��һ�������Ǹ�ֵ          ���ɻ󣿣�������������
end


Trend_rcp45 = []; Rsq_rcp45 = [];  pvalue_rcp45 = [];    %�������յ�Ҫ������ݵ�����
Trend_rcp85 = []; Rsq_rcp85 = [];  pvalue_rcp85 = [];

%p=parpool(4);
%prc = -9999*ones([nl np]); tas = -9999*ones([nl np]); hum = -9999*ones([nl np]); clo = -9999*ones([nl np]);   %�����Сһ���ĸ�������������
outdata = -9999*ones([nl np]);
XM = [];
for s = 1 : 2    %ѭ��  {'RCP45','RCP85'}
    xm(:, 1) = yr1:yr2;
    xs(:, 1) = yr1:yr2;
    xc(:, 1) = yr1:yr2;
    
    for yr = yr1 : yr2
        for v = 1 : 4   %ѭ���ĸ���������
            ff = [wks  '/' scenario{s} '/' vflx{v} '/' num2str(yr)];   %��ȡRCP45��RCP85�µ�2006-2009���ÿ��������������
            disp(ff);
            dat = load(ff);   %�����Ӧ������
            [nd, md] = size(dat);    %��ȡ����(����)�Ĵ�С��ndΪ�У�mdΪ��
            xd = dat(:, 3:end);    %��ȡ���ݴӵ����п�ʼ������
            
            for k = 1 : n
                if v == 1       %��ˮ���ۼ�ֵ�����������¶�����ƽ��ֵ
                    data = sum(xd, 2);     %�������
                    outdata(lat(k),lon(k)) = data(k);
                else
                    data = mean(xd, 2);   %������ƽ��
                    outdata(lat(k),lon(k)) = data(k);
                end
            end
            %outf = [owk '/' scenario2{s} '/' vflx{v} num2str(yr), '.tif'];
            %geotiffwrite(outf,single(outdata), R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
            xdat =outdata;
            for i = 1 : 5
                if i < 5
                    x1 = xdat(msk == i & xdat > -900);    %xdat��218418����Ԫ��ȥ���й������������Ԫ���й�����ڵ�ֻ��95483��?
                    sub_area = carea_img(msk == i & xdat > -900);    %%%msk�����������carea_imag����Ԫ����Ŀռ�����
                else
                    x1 = xdat(msk < 100 & xdat > -900);    %���й�������255����ȫ������1.2.3.4����һ��
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
    hdrr={'��ظ�ԭ��','�ȴ�-���ȴ�������','�´�������','�´���½��','ȫ��½��'};
    for j=1:4
        xlswrite([owk, filesep, scenario{s}, '_', vflx{j}, '.xls'], vflxxx{j});%��ֵ�ͱ�׼��.cv��ͬд��%%%
    end
    
    %
    XM = [XM, xm];%�����龰�ľ�ֵ?
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
        
        xlabel('\fontname{����}��?'); ylabel([vflxx{v}, ' (', uflx{v}, ')']);
        
        legend('\fontname{Times New Roman}RCP45','\fontname{Times New Roman}RCP85');
        legend('boxoff');
        title(hdrr{r});
        set(gca, 'FontSize', 14);
        fig = [owk, '/' vflx{v} '/' scenario{s},'_', vflx{v}, '_', hdr{r + 1} '.jpg'];
        saveas(gca, fig);
        
        Trend_rcp45(r, v) = b1(2);  Rsq_rcp45(r, v) = p1(1);  pvalue_rcp45(r, v) = p1(3);      %���� R2 Pֵ?
        Trend_rcp85(r, v) = b2(2);  Rsq_rcp85(r, v) = p2(1);  pvalue_rcp85(r, v) = p2(3);
    end
    
end

T_45 = [Trend_rcp45, Rsq_rcp45, pvalue_rcp45];
T_85 = [Trend_rcp85, Rsq_rcp85, pvalue_rcp85];

fxls = [owk, '/'  scenario{1},'_meteo',  '_trend_r2_p0711', '.xls'];
xlswrite(fxls, T_45);

fxls = [owk, '/'  scenario{2},'_meteo',  '_trend_r2_p0711', '.xls'];
xlswrite(fxls, T_85);
