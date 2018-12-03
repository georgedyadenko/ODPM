%
% codename: MasterPiece
%

target_futures = 'RIZ2';

max_distance_from_spot = 10000;

globals
cpt_import

Mo = cpto;
Mf = cptf;

if exist('Mo','var')
    for i = 1:1:length(Mo)
        for j = 1:1:length(Mf)
            if strcmp(Mo(i).baseactive, Mf(j).active)
                Mo(i).baseactive_index = j;
                break
            end;
        end;
    end;
else
    display('error: no options prices');
end;

OPT.matrix = [];

for i = 1:1:length(Mo)    
    if strcmp(Mo(i).baseactive, target_futures)
        OPT.matrix(size(OPT.matrix,1)+1,:) = [i Mo(i).strike strcmp(Mo(i).type, 'Call') Mo(i).bidprice Mo(i).askprice-Mo(i).bidprice Mo(i).bidvol Mo(i).askvol Mo(i).implvol Mo(i).tomat];
    end;
end;

OPT.matrix = sortrows(OPT.matrix, [2 3]);

for i = size(OPT.matrix,1):-1:1
    if abs(OPT.matrix(i,2) - Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice) > max_distance_from_spot
        %OPT.matrix(i,:)=[];
    end;
end;

%здесь появляется ограничение, что матрица опционов должна быть однородной,
%на один и тот же фьючерсный контракт
%индекс контракта в таблице Mf берется исходя из 1го опциона
OPT.matrix(size(OPT.matrix,1)+1,1:9)=[Mo(OPT.matrix(i,1)).baseactive_index Mf(Mo(OPT.matrix(i,1)).baseactive_index).askprice -1 Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice Mf(Mo(OPT.matrix(i,1)).baseactive_index).askprice-Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidvol Mf(Mo(OPT.matrix(i,1)).baseactive_index).askvol 0 Mf(Mo(OPT.matrix(i,1)).baseactive_index).tomat];

for i = 1:1:size(OPT.matrix,1)
    
    if OPT.matrix(i, 3) ~= -1
        [CallDelta, PutDelta] = blsdelta(Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice, OPT.matrix(i,2), 0.0001, OPT.matrix(i,9)/365, OPT.matrix(i,8)/100);
        g = blsgamma(Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice, OPT.matrix(i,2), 0.0001, OPT.matrix(i,9)/365, OPT.matrix(i,8)/100);
        [CallTheta, PutTheta] = blstheta(Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice, OPT.matrix(i,2), 0.0001, OPT.matrix(i,9)/365, OPT.matrix(i,8)/100);
        v = blsvega(Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice, OPT.matrix(i,2), 0.0001, OPT.matrix(i,9)/365, OPT.matrix(i,8)/100);
    end;
        
    switch OPT.matrix(i, 3)
        case 1
            OPT.matrix(i, 10) = CallDelta;
            OPT.matrix(i, 11) = g;            
            OPT.matrix(i, 12) = CallTheta/100;      
            OPT.matrix(i, 13) = v/100;
        case 0
            OPT.matrix(i, 10) = PutDelta;
            OPT.matrix(i, 11) = g;             
            OPT.matrix(i, 12) = PutTheta/100;
            OPT.matrix(i, 13) = v/100;
        case -1
            OPT.matrix(i, 10) = 1;
            OPT.matrix(i, 11) = 0;  
            OPT.matrix(i, 12) = 0; 
            OPT.matrix(i, 13) = 0;             
    end;
    
end;

OPT.my_port = zeros(3, size(OPT.matrix,1));
OPT.algo = -1;

basic_vector = zeros(1, size(OPT.matrix,1));
my_options = optimset('LargeScale','off', 'MaxFunEvals',1000,'Display','off','Algorithm','active-set','DiffMinChange',1,'DiffMaxChange',10,'PlotFcns',@optimplotx);
[solution_vector, z_values] = ...
    fmincon(@objfun, basic_vector, [],[],[],[],[],[], @confun, my_options);

my_plot(solution_vector);

%{
for i = 1:1:size(OPT.matrix,1)
    plot_range = [Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice-50000 ...
        Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice+50000];
    step = (plot_range(2)-plot_range(1))/200;
    x = plot_range(1):step:plot_range(2);
    [call_price, put_price] = blsprice(x, OPT.matrix(i,2), 0.001, 1/365, 0.3);
    switch OPT.matrix(i, 3)
        case 1
            %if you need net result
            option_price = OPT.matrix(i,5).*(call_price - OPT.matrix(i,4));
            
            %if you need option price
            %option_price = OPT.matrix(i,5).*call_price;
        case 0
            %if you need net result
            option_price = OPT.matrix(i,5).*(put_price  - OPT.matrix(i,4));
            
            %if you need option price
            %option_price = OPT.matrix(i,5).*put_price;
        case -1
            option_price = zeros(size(call_price));
    end;
    
    %theoretical bullshit, must be removed in future version
    Z = normpdf(x, Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice, ...
        (Mf(Mo(OPT.matrix(i,1)).baseactive_index).bidprice)*(OPT.matrix(i,6)/100)*sqrt(OPT.matrix(i,7)/365));
    
    OPT.matrix(i, 8) = trapz(Z.*option_price)*step;  
    
    if isnan(OPT.matrix(i, 8))
        OPT.matrix(i, 8) = 0;
    end;
    
    OPT.matrix(i, 9) = min(option_price);
    
    positive_option_price = option_price;
    for j = 1:1:length(option_price)
        if option_price(j) < 0
            positive_option_price(j) = 0;
        end;
    end;
    
    OPT.matrix(i, 10) = trapz(Z.*positive_option_price)*step;  
    
    negative_option_price = option_price;
    for j = 1:1:length(option_price)
        if option_price(j) > 0
            negative_option_price(j) = 0;
        end;
    end;
    
    OPT.matrix(i, 11) = trapz(Z.*negative_option_price)*step;   
    


end;

clear plot_range
clear step
clear x
clear call_price
clear put_price
clear Z
clear theory_price
clear option_price
clear j
clear positive_option_price
clear negative_option_price

my_port = zeros(3, size(OPT.matrix,1));

basic_vector = zeros(1, size(OPT.matrix,1));
my_Mo = optimset('LargeScale','off','MaxFunEvals',1000,'Display','off','Algorithm','active-set');
[solution_vector, z_values] = ...
    fmincon(@optimum_objfun, basic_vector, [],[],[],[],[],[], @optimum_confun, my_Mo);

solution_vector=round(solution_vector);

my_port(3,:) = (sign(my_port(1,:))~=sign(solution_vector)).*(min(abs(my_port(1,:)),abs(solution_vector))).* ...
    (-my_port(2,:)+(OPT.matrix(:,4)'+(solution_vector<0).*spread')).*sign(my_port(1,:));
my_port(2,:) = (sign(my_port(1,:))~=sign(solution_vector)).*(abs(my_port(1,:))<abs(solution_vector)).*(OPT.matrix(:,4)'+(solution_vector<0).*spread')+...
    (sign(my_port(1,:))~=sign(solution_vector)).*(abs(my_port(1,:))>abs(solution_vector)).*my_port(2,:)+...
    (sign(my_port(1,:))==sign(solution_vector)).*(my_port(2,:).*my_port(1,:)+(OPT.matrix(:,4)'+(solution_vector<0).*spread').*solution_vector)./(solution_vector+my_port(1,:));
my_port(1,:) = solution_vector+my_port(1,:);

for i = 1:1:size(my_port,1)
    for j = 1:1:size(my_port,2)
        if isnan(my_port(i,j))
            my_port(i,j)=0;
        end;
    end;
end;

%my_port(2,:) = OPT.matrix(:,4)'+(solution_vector<0).*spread';
solution_vector = basic_vector;

%}

%display(sprintf('%d \t', [OPT.matrix(:,2)]));
%display(sprintf('%d \t\t', [OPT.matrix(:,3)]));
%display(sprintf('%d \t\t', [OPT.matrix(:,5)]));
%display(sprintf('%3.2f \t', [(round(solution_vector'.*100)/100)]));

% ***************************************************************
% HEREBY IS REPRESENTED THE OLD CODE FOR TEST OF Mo STRATEGY
% ***************************************************************
%{
global price
global i
global days
global OM
global p

days = 50;
%{
hourly_returns = randn(8*days,1);
price = 140000;
for i = 1:1:size(hourly_returns, 1)
    price(i+1,1) = price(i,1)*(1+hourly_returns(i,1)/100);
end;
%}

for I = 1:30:days*8-5
    
p = zeros(1,15);
result=0;

for i = 1 : 1 : I%length (price)-3
spot_price = 140000;%round( price(i) / 5000 ) * 5000;

%Mo OPT.matrix
OM = [];
OM(1,1) = spot_price - 15000;
OM(2,1) = spot_price - 10000;
OM(3,1) = spot_price - 5000;
OM(4,1) = spot_price;
OM(5,1) = spot_price + 5000;
OM(6,1) = spot_price + 10000;
OM(7,1) = spot_price + 15000;

OM(8,1) = spot_price - 15000;
OM(9,1) = spot_price - 10000;
OM(10,1) = spot_price - 5000;
OM(11,1) = spot_price;
OM(12,1) = spot_price + 5000;
OM(13,1) = spot_price + 10000;
OM(14,1) = spot_price + 15000;

OM(15,1) = price(i);

for j = 1 : 1 : 7
    OM(j, 2) = blsprice(price(i), OM(j, 1), 0.05, (days-i/8)/365, 0.3036);
end;
for j = 8 : 1 : 14
    [fake, OM(j, 2)] = blsprice(price(i), OM(j, 1), 0.05, (days-i/8)/365, 0.3036);
end;
    OM(15,2) = price(i);

for j = 1 : 1 : 7
    OM(j, 4) = blsdelta(price(i), OM(j, 1), 0.05, (days-i/8)/365, 0.3036);
end;
for j = 8 : 1 : 14
    [fake, OM(j, 4)] = blsdelta(price(i), OM(j, 1), 0.05, (days-i/8)/365, 0.3036);
end;
    OM(15,4) = 1;

for j = 1 : 1 : 7%size(OM,1)
    plot_range = [price(i)-50000 price(i)+50000];
    step = (plot_range(2)-plot_range(1))/200;
    y = plot_range(1):step:plot_range(2);
    [zC, zP] = blsprice(y, OM(j,1), 0.05, 1/365, 0.3036);
    zC = zC - OM(j,2);
    %zP = zP - OM(j,3);
    Z = normpdf(y, price(i), price(i)*0.3036*sqrt((days-(i+1)/8)/365));
    OM(j, 3) = trapz(Z.*zC)*step;
    %OM(j, 5) = trapz(Z.*zP)*step;
end;
for j = 8 : 1 : 14%size(OM,1)
    plot_range = [price(i)-50000 price(i)+50000];
    step = (plot_range(2)-plot_range(1))/200;
    y = plot_range(1):step:plot_range(2);
    [zC, zP] = blsprice(y, OM(j,1), 0.05, 1/365, 0.3036);
    %zC = zC - OM(j,2);
    zP = zP - OM(j,2);
    Z = normpdf(y, price(i), price(i)*0.3036*sqrt((days-(i+1)/8)/365));
    %OM(j, 3) = trapz(Z.*zC)*step;
    OM(j, 3) = trapz(Z.*zP)*step;
end;
    plot_range = [price(i)-50000 price(i)+50000];
    step = (plot_range(2)-plot_range(1))/200;
    y = plot_range(1):step:plot_range(2);
    zP = -price(i)+y;
    Z = normpdf(y, price(i), price(i)*0.3036*sqrt((days-(i+1)/8)/365));
    OM(15, 3) = trapz(Z.*zP)*step;

x=zeros(1,size(OM,1));

Mo = optimset('LargeScale','off','MaxFunEvals',1000,'Display','off','Algorithm','active-set');
[x1, fval1] = fmincon(@optimum_objfun,x,[],[],[],[],[],[],@optimum_confun,Mo);

x1 = round(x1*100);

if isnan(x1)
    x1 = p;
end;
%{
if i < 50
    x1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
else
    x1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
end;
%}
if optimum_objfun(x1/100)~=optimum_objfun(p/100)%result - (x1-p)*OM(:,2) < result
    result = result - (x1-p)*OM(:,2);
    p=x1;
%    display([num2str(round(result)) ' ' num2str(round(p)) ' <' num2str(round(price(i)))]);
%    display([num2str(round(result)) ' ' num2str(round(OM(:,2)'))]);
end;
    display([num2str(i) ' ' num2str(round(result)) ' ' num2str(round(p)) ' <' num2str(round(price(i)))]);
%    display([num2str(round(result)) ' ' num2str(round(OM(:,2)'))]);
end;

port=zeros(size(zC));
pport=0;
for j = 1:1:7
    plot_range = [price(i)-50000 price(i)+50000];
    step = (plot_range(2)-plot_range(1))/200;
    y = plot_range(1):step:plot_range(2);
    [zC, zP] = blsprice(y, OM(j,1), 0.05, (days-i/8)/365, 0.3036);
    port = port + p(j)*(zC);

    [zC, zP] = blsprice(price(i), OM(j,1), 0.05, (days-i/8)/365, 0.3036);
    pport = pport + p(j)*(zC);
end;
for j = 8:1:14
    plot_range = [price(i)-50000 price(i)+50000];
    step = (plot_range(2)-plot_range(1))/200;
    y = plot_range(1):step:plot_range(2);
    [zC, zP] = blsprice(y, OM(j,1), 0.05, (days-i/8)/365, 0.3036);
    port = port + p(j)*(zP);%-OM(j,2));

    [zC, zP] = blsprice(price(i), OM(j,1), 0.05, (days-i/8)/365, 0.3036);
    pport = pport + p(j)*(zP);%-OM(j,2));
end;
    j = 15;
    plot_range = [price(i)-50000 price(i)+50000];
    step = (plot_range(2)-plot_range(1))/200;
    y = plot_range(1):step:plot_range(2);
    zP = y;%-price(i)+y;
    port = port + p(j)*(zP);%-OM(j,2));

    zP = price(i);%-price(i)+price(length(price));
    pport = pport + p(j)*(zP);%-OM(j,2));


%hold off
plot(y,(port+result));
xx=zeros(size(y));
hold all
plot(y,xx)
plot(price(i),pport+result,'*')
text(price(i)+5000,pport+result,[num2str(round(pport+result))])
%round(result + (price(length(price))*7 - sum(p(1:7)*OM(1:7,1))) + (-price(length(price))*7 + sum(p(8:14)*OM(8:14,1))))

end;
%}