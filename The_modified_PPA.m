% This program calculate the transportation cost of each supplier in a given graph; 
% ------------------------------------------
% output variables:
% cost_total(i): the cost the Supplier S_i;
% time: the computation time;
% Q0: the total link flow distribution
% ------------------------------------------

clc;clear;close all;
load graph1.mat

% convergence criterion
eps=0.1; % eps is usually set between 0.0001 and 1 depending on the problem

% set-up the effective distance parameters 
m1=5;m2=3;m3=0.3;% effective distance = unit cost = m1+m2*exp(-m3*Q);

% L(:,:,i)stores the link lengths of subgraph G_i
L0(L0==0)=inf;
for i=1:S
    L(:,:,i)=L0; 
end

% A(:,i) stores the the OD pair demand for subgraph G_i
for i=1:S
    A(:,i)=Dmd(:,i);
end

% Initialization
D_ini=ones(N);
D_ini=D_ini-diag(diag(D_ini));
D=zeros(N,N,S);
for i=1:S
    D(:,:,i)=D_ini;
end

for i=1:S
    temp_D(:,:,i)=zeros(N);
end


%% start the iteration
ite=0; % the iteration record

tic;
while (sum(sum(sum( abs(temp_D-D) ))))>eps
    for i=1:S
        temp_D(:,:,i)=D(:,:,i);
    end
    
    for i=1:S
        B=D(:,:,i)./L(:,:,i);
        B=B-diag(sum(B));
        B(:,1)=[];
        P=pinv(B)*A(:,i);
        P=[0;P];
        tempP=repmat(P,1,N)-repmat(P',N,1);
        Q(:,:,i)=(D(:,:,i)./L(:,:,i)).*tempP; % the flow matrix
        Q(:,:,i)=abs(Q(:,:,i)); % undirected graph
        D(:,:,i)=(Q(:,:,i)+D(:,:,i))./2;
    end
   
    % update the total flow on each link
    Q0=zeros(N);
    for i=1:S
        Q0=Q0+Q(:,:,i);
    end
 
    % effective distance on each link
    C_link_unit=m1+m2*exp(-m3*Q0);
    
    %update the SO cost function on each link of each supplier's network
    for i=1:S
        C_SO(:,:,i)=C_link_unit-m3*m2*Q(:,:,i).*exp(-m3*Q0); % convert to SO problem
        L(:,:,i)=C_SO(:,:,i); % update the link length
        temp=L(:,:,i);
        temp(L0==inf)=inf;
        L(:,:,i)=temp;
    end  

    ite=ite+1;
end
% iterations terminate

% calculate the cost for each suppliers
for i=1:S
    C_link(:,:,i)=C_link_unit.*(Q(:,:,i)); % cost of links for each subgraph
    cost_total(i)=sum(sum( C_link(:,:,i) ))/2; % each customer's total cost
end

time=toc;

cost_total
time


