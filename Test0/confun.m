function [c, ceq] = confun(x)
globals
solution=x+OPT.my_port(1,:);

plot_range = [min(OPT.matrix(:,2))-50000 max(OPT.matrix(:,2))+50000];
step = (plot_range(2)-plot_range(1))/200;
xx = plot_range(1):step:plot_range(2);
portfolio = zeros(size(xx));    
for i = 1:1:size(OPT.matrix,1)
    switch OPT.matrix(i,3)
        case 1 % if call
            [call_price, put_price] = blsprice(xx, OPT.matrix(i,2), 0.001, 0.1/365, 0.3);
            portfolio = portfolio + solution(i).*OPT.matrix(i,5).*(call_price - OPT.matrix(i,4));
        case 0 % if put
            [call_price, put_price] = blsprice(xx, OPT.matrix(i,2), 0.001, 0.1/365, 0.3);
            portfolio = portfolio + solution(i).*OPT.matrix(i,5).*(put_price - OPT.matrix(i,4));
        case -1
            portfolio = portfolio + solution(i).*OPT.matrix(i,5).*(xx - OPT.matrix(i,4));                
    end;
end;

c = [solution-10, -solution-10, min(portfolio)+300000]';
ceq = [];

%{
% Nonlinear inequality constraints <0

% complex constraints
% -10<x<10
% 0<max_loss<100000

    solution=x+my_port(1,:);
    
    plot_range = [min(matrix(:,2))-50000 max(matrix(:,2))+50000];
    step = (plot_range(2)-plot_range(1))/200;
    xx = plot_range(1):step:plot_range(2);
    portfolio = zeros(size(xx));    
    for i = 1:1:size(matrix,1)
        switch matrix(i,3)
            case 1 % if call
                [call_price, put_price] = blsprice(xx, matrix(i,2), 0.001, 1/365, 0.3);
                portfolio = portfolio + solution(i).*matrix(i,5).*(call_price - matrix(i,4));
            case 0 % if put
                [call_price, put_price] = blsprice(xx, matrix(i,2), 0.001, 1/365, 0.3);
                portfolio = portfolio + solution(i).*matrix(i,5).*(put_price - matrix(i,4));
            case -1
                portfolio = portfolio + solution(i).*matrix(i,5).*(xx - matrix(i,4));                
        end;
    end;


c = [solution-10, -solution-10, -min(portfolio)-300000];%, -min(portfolio)-100000];%, x*matrix(:,9), -x*matrix(:,9)-40000]';

% simple covenant 0<x<10
%c = [x-10, -x]';

% Nonlinear equality constraints
ceq = [];

%uncomment if you want delta-neutral strategy
%optimum_globals
%ceq = [x * matrix(:, ???)];
%}


% ***************************************************************
% HEREBY IS REPRESENTED THE OLD CODE FOR TEST OF OPTIONS STRATEGY
% ***************************************************************
%{
%global ANT
%global DNT
%f = (exp(x*ANT(:,5)+DNT(:,4)'*DNT(:,5))+exp(-x*ANT(:,5)-DNT(:,4)'*DNT(:,5))-2)+...
%    (exp((x*ANT(:,6)+DNT(:,4)'*DNT(:,6))/50)+exp((-x*ANT(:,6)-DNT(:,4)'*DNT(:,6))/50)-2);
%f = (exp(x*ANT(:,5)+DNT(:,4)'*DNT(:,5))+exp(-x*ANT(:,5)-DNT(:,4)'*DNT(:,5))-2)+...
%    (exp((x*ANT(:,6)+DNT(:,4)'*DNT(:,6))/50)+exp((-x*ANT(:,6)-DNT(:,4)'*DNT(:,6))/50)-2);
lowbound=-10000;
%f = exp( (lowbound -p*OM(:,3) - x*OM(:,3)) / 50000 ) + exp( -(lowbound -p*OM(:,3) - x*OM(:,3)) / 50000 )-2;
%f = exp( (lowbound -p*OM(:,3) - x*OM(:,3)) / 50000 ) + exp( -(lowbound -p*OM(:,3) - x*OM(:,3)) / 50000 )-2;

plot_range = [price(i)-50000 price(i)+50000];
step = (plot_range(2)-plot_range(1))/200;
y = plot_range(1):step:plot_range(2);

port=zeros(size(y));
for j = 1:1:7
    plot_range = [price(i)-50000 price(i)+50000];
    step = (plot_range(2)-plot_range(1))/200;
    y = plot_range(1):step:plot_range(2);
    [zC, ~] = blsprice(y, OM(j,1), 0.05, 0.01/365, 0.3036);
    port = port + x(j)*(zC - OM(j,2));
end;
for j = 8:1:14
    plot_range = [price(i)-50000 price(i)+50000];
    step = (plot_range(2)-plot_range(1))/200;
    y = plot_range(1):step:plot_range(2);
    [~, zP] = blsprice(y, OM(j,1), 0.05, 0.01/365, 0.3036);
    port = port + x(j)*(zP - OM(j,2));
end;
    j = 15;
    plot_range = [price(i)-50000 price(i)+50000];
    step = (plot_range(2)-plot_range(1))/200;
    y = plot_range(1):step:plot_range(2);
    zP = -price(i)+y;
    port = port + x(j)*zP;

    %a = max(port);
% Nonlinear inequality constraints <0
c = [x-20 ,-x-20]';%, a-1000]';
% Nonlinear equality constraints
ceq = [];
%}
