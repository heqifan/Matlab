close all;

% ����ͼ�����ԣ�����ͼλ�úͳߴ�
% set(gcf,'units','inches','position',[0,10.0,14.0,10.0]);
% set(gcf,'DefaultAxesFontSize',18); % �����������С
%��ȡ���ݣ�sd rmse �� r��
data=xlsread('Taylor.xlsx','Sheet1','B2:P18');%�ļ�·��  

% sdev = data(:,1);
% crmsd = data(:,2);
% ccoef = data(:,3);
% %mmodel ID���������ֶ���������ΪҪÿ���������ñ�־
ID = {'Obs','MODIS','MuSyQ','GLOCV','GLASS','Mean','Median','Weight',...
'Mul','Bag','Ada','Grad',...
'Stack','RF','LCE','Vote','BMA'};
label = ID;
% %>>���� taylor_diagram

% % �����ļ�
% writepng(gcf,'taylor fig.png');

%%�����������
for iserie = 2 : size(data,1)
 		    S = allstats(data(1,:),data(iserie,:));
 		    MYSTATS(iserie,:) = S(:,2); % We get stats versus reference
end       %for iserie
 MYSTATS(1,:) = S(:,1); % We assign reference stats to the first row
 
%  taylordiag(MYSTATS(:,2),MYSTATS(:,3),MYSTATS(:,4));
stdev=MYSTATS(1,2);

MYSTATS(:,2)=MYSTATS(:,2)/stdev;% ��׼������
MYSTATS(:,3)=MYSTATS(:,3)/stdev;%��׼������

 [hp ht axl] =taylordiag(MYSTATS(:,2),MYSTATS(:,3),MYSTATS(:,4), ...
                 'markerLabel',label, 'markerLegend', 'on','tickrms',[0:.2:1],'titleRMS', 1 ,'showlabelsRMS',1,... 
                 'widthRMS',1,'colRMS','r',...
                 'tickSTD',[0:.25:1.25],'limSTD',1.25,'styleSTD','-',... 
                 'tickCOR',[.1:.1:.9 .95 .99],'showlabelsCOR',1,'titleCOR',1);
% writepng(gcf,'D:\WareHouse\Matlab\Matlab\taylor fig.png');
%  [hp, ht, axl] = taylordiag(MYSTATS(:,2),MYSTATS(:,3),MYSTATS(:,4), ...
%     'markerLabel',label, 'markerLegend', 'on', ...
%     'styleSTD', '-', 'colOBS','r', 'markerObs','o', ...
%     'markerSize',12, 'tickRMS',[0:1:5],'limSTD',5, ...
%     'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
%     'titleRMS','on', 'titleOBS','Observation');

% [hp, ht, axl] = taylor_diagram(MYSTATS(:,2),MYSTATS(:,3),MYSTATS(:,4), ...
%     'markerLabel',label, 'markerLegend', 'on', ...
%     'styleSTD', '-', 'colOBS','r', 'markerObs','o', ...
%     'markerSize',12, 'tickRMS',[0:1:5],'limSTD',5, ...
%     'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
%     'titleRMS','on', 'titleOBS','Observation');

