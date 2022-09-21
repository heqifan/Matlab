function [wk, sigmak, loglikelihood] = codeBMA(Tainingdata)
% clc
% clear
% 
% load Tainingdata
% Tainingdata=Tainingdata4;
% measured  AA	SW	PMv	PM	PT




IterL=100;                     % iteration length
K=size(Tainingdata,2)-1;         % number of model
T=size(Tainingdata,1);         % length of training data

wk=zeros(IterL,K);
sigmak=zeros(IterL,K);
loglikelihood=zeros(IterL,1);


wk(1,:)=1/K;                    % Intialization wk
for i=1:K
    temp=(Tainingdata(:,i+1)-Tainingdata(:,1));
    sigmak(1,i)=std(temp);
end

temp2=0;
for t=1:T
    for i=1:K
        temp(t,i)=normpdf(Tainingdata(t,1),Tainingdata(t,i+1),sigmak(1,i)); 
    end
    temp2=temp2+log(sum(temp(t,:).*wk(1,:)));   
end
loglikelihood(1,1)=temp2;

for Iter=2:IterL
    for t=1:T
        for k=1:K
            temp2(t,k)=normpdf(Tainingdata(t,1),Tainingdata(t,k+1),sigmak(Iter-1,k));
        end
        for k=1:K
            Z(t,k)=temp2(t,k)/sum(temp2(t,:));
        end
    end
    for k=1:K
        wk(Iter,k)=sum(Z(:,k))/T;
        sigmak(Iter,k)=sqrt(sum(Z(:,k).*(Tainingdata(:,1)-Tainingdata(:,k+1)).^2)/sum(Z(:,k)));
    end
    temp2=0;
    for t=1:T
        for i=1:K
            temp(t,i)=normpdf(Tainingdata(t,1),Tainingdata(t,i+1),sigmak(Iter,i)); 
        end        
        temp2=temp2+log(sum(temp(t,:).*wk(Iter,:)));   
    end           
       loglikelihood(Iter,1)=temp2;       
end



        
        
    
    
    
    
    