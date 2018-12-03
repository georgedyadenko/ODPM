function f = objfun(x)
    globals;

    switch algo
        case -1
            f = (1 /( exp( (x+my_port(1,:))*OPT.matrix(:,10)/2) + exp( -(x+OPT.my_port(1,:))*OPT.matrix(:,10)/2)-2+0.1 )) ...
                * (- (x+OPT.my_port(1,:))*OPT.matrix(:,12)) ;
        case 1
            % delta = 0
            % theta -> max
            f = (1 /( exp( (x+my_port(1,:))*matrix(:,12)/2) + exp( -(x+my_port(1,:))*matrix(:,12)/2)-2+0.1 )) ...
                * (- (x+my_port(1,:))*matrix(:,13)) ;
            %display('**************************');
            %(x)
            %round((x+my_port(1,:))*matrix(:,13))
            %round(f) 
        case 5
            % delta = 0
            % theta -> max
            f =       -(  (x+my_port(1,:))*matrix(:,8));%(1 /( exp( (x+my_port(1,:))*matrix(:,12)/2) + exp( -(x+my_port(1,:))*matrix(:,12)/2)-2+0.1 )) ...
                %* (- (x+my_port(1,:))*matrix(:,13)) ...
           
            %display('**************************');
            %(x)
            %round((x+my_port(1,:))*matrix(:,13))
            %round(f)             
        case 2
            % delta = 0
            % vega  = 0
            % theta -> max


            f = (1 /( exp( (x+my_port(1,:))*matrix(:,12)/2) + exp( -(x+my_port(1,:))*matrix(:,12)/2)-2+0.1 )) ...
                * (- (x+my_port(1,:))*matrix(:,13)) ...
                * (1/( exp((x+my_port(1,:))*matrix(:,14)/1000)+exp(-(x+my_port(1,:))*matrix(:,14)/1000)-2+0.1 ));

        case 3
            f = (1 /( exp( (x+my_port(1,:))*matrix(:,12)/2) + exp( -(x+my_port(1,:))*matrix(:,12)/2)-2+0.1 )) ...
                * (- (x+my_port(1,:))*matrix(:,13)) ...
                * (-sum((sign(my_port(1,:))~=sign(x)).*(min(abs(my_port(1,:)),abs(x))).* ...
                (-my_port(2,:)+(matrix(:,4)'+(x<0).*spread')).*sign(my_port(1,:))));            
        case 4
            plot_range = [futures(options(matrix(1,1)).baseactive_index).bidprice-10000 futures(options(matrix(1,1)).baseactive_index).bidprice+10000];
            step = (plot_range(2)-plot_range(1))/200;
            Ox = plot_range(1):step:plot_range(2);
            
            old_portfolio = zeros(size(Ox));
            for i = 1:1:size(matrix,1)
                [call_price, put_price] = blsprice(Ox, matrix(i,2), 0.001, 0.1/365, matrix(i,6)/100);
                switch matrix(i,3)
                    case 1 % if call
                        old_portfolio = old_portfolio + ...
                            my_port(1,i).*(call_price - my_port(2,i)) ...
                            + my_port(3,i);
                    case 0 % if put
                        old_portfolio = old_portfolio + ...
                            my_port(1,i).*(put_price - my_port(2,i)) ...
                            + my_port(3,i);
                    case -1
                        old_portfolio = old_portfolio + ...
                            my_port(1,i).*(Ox - my_port(2,i)) ...
                            + my_port(3,i);             
                end;
            end;

            solution_vector=x;
            new_port(3,:) = my_port(3,:)+(sign(my_port(1,:))~=sign(solution_vector)).*(min(abs(my_port(1,:)),abs(solution_vector))).* ...
                (-my_port(2,:)+(matrix(:,4)'+(solution_vector<0).*spread')).*sign(my_port(1,:));
            new_port(2,:) = (sign(my_port(1,:))~=sign(solution_vector)).*(abs(my_port(1,:))<abs(solution_vector)).*(matrix(:,4)'+(solution_vector<0).*spread')+...
                (sign(my_port(1,:))~=sign(solution_vector)).*(abs(my_port(1,:))>abs(solution_vector)).*my_port(2,:)+...
                (sign(my_port(1,:))==sign(solution_vector)).*(my_port(2,:).*my_port(1,:)+(matrix(:,4)'+(solution_vector<0).*spread').*solution_vector)./(solution_vector+my_port(1,:)).*((solution_vector+my_port(1,:))~=0);
            new_port(1,:) = solution_vector+my_port(1,:);
            for i = 2:1:2
                for j = 1:1:size(new_port,2)
                    if isnan(new_port(i,j))
                        new_port(i,j)=0;
                    end;
                end;
            end;
            new_portfolio = zeros(size(Ox));
            for i = 1:1:size(matrix,1)
                [call_price, put_price] = blsprice(Ox, matrix(i,2), 0.001, 0.1/365, matrix(i,6)/100);
                switch matrix(i,3)
                    case 1 % if call
                        new_portfolio = new_portfolio + ...
                            new_port(1,i).*(call_price - new_port(2,i)) ...
                            + new_port(3,i);
                    case 0 % if put
                        new_portfolio = new_portfolio + ...
                            new_port(1,i).*(put_price - new_port(2,i)) ...
                            + new_port(3,i);
                    case -1
                        new_portfolio = new_portfolio + ...
                            new_port(1,i).*(Ox - new_port(2,i)) ...
                            + new_port(3,i);             
                end;
            end;
            
            Z = normpdf(Ox, futures(options(matrix(1,1)).baseactive_index).bidprice, ...
                (futures(options(matrix(1,1)).baseactive_index).bidprice)*(matrix(3,6)/100));%*sqrt((matrix(3,7)/365))); %nujno ispravit'

            %hold off
            %plot(Ox,new_portfolio);
            %hold on
            %plot(Ox,old_portfolio);
            
            f=-trapz(Z.*(new_portfolio-old_portfolio))*step;
            
            f = (1 /( exp( (x+my_port(1,:))*matrix(:,12)/2) + exp( -(x+my_port(1,:))*matrix(:,12)/2)-2+0.1 )) ...
                * (- trapz(Z.*(new_portfolio))*step) ;
        case 0

            solution=x+my_port(1,:);

            plot_range = [futures(options(matrix(1,1)).baseactive_index).bidprice-2000 futures(options(matrix(1,1)).baseactive_index).bidprice+2000];
            step = (plot_range(2)-plot_range(1))/200;
            xx = plot_range(1):step:plot_range(2);
            portfolio = zeros(size(xx));

            for i = 1:1:size(matrix,1)
                switch matrix(i,3)
                    case 1 % if call
                        [call_price, put_price] = blsprice(xx, matrix(i,2), 0.001, matrix(i,7)/365, matrix(i,6)/100);
                        portfolio = portfolio + ...
                                x(i).*matrix(i,5).*(call_price - matrix(i,4));
                                %my_port(i).*(call_price - my_price(i)) + ...
                                %(my_port(i)<0)*spread(i)
                                %my_fixedresult + ...
                    case 0 % if put
                        [call_price, put_price] = blsprice(xx, matrix(i,2), 0.001, matrix(i,7)/365, matrix(i,6)/100);
                        portfolio = portfolio + solution(i).*matrix(i,5).*(put_price - matrix(i,4));
                    case -1
                        portfolio = portfolio + solution(i).*matrix(i,5).*(xx - matrix(i,4));
                end;
            end;

            Z = normpdf(x, futures(options(matrix(1,1)).baseactive_index).bidprice, ...
                (futures(options(matrix(1,1)).baseactive_index).bidprice)*(matrix(3,6)/100)*sqrt(matrix(3,7)/365)); %nujno ispravit'

            %f=-trapz(Z.*portfolio)*step;

            %f=trapz(portfolio);

            x=solution;
            %f = (1/( exp((x.*matrix(:,5)')*matrix(:,12)/2)+exp(-(x.*matrix(:,5)')*matrix(:,12)/2)-2+0.1 )) ...
            %    * (-(x.*matrix(:,5)')*matrix(:,13)) ...
            %    * (1/( exp((x.*matrix(:,5)')*matrix(:,14)/1000)+exp(-(x.*matrix(:,5)')*matrix(:,14)/1000)-2+0.1 )); % ...

            f = (1 /( exp( (x+my_port(1,:))*matrix(:,12)/2) + exp( -(x+my_port(1,:))*matrix(:,12)/2)-2+0.1 )) ...
                * (- (x+my_port(1,:))*matrix(:,13)) ;

            %f= (-trapz(Z.*portfolio)*step);
            %* (-trapz(Z.*portfolio)*step) ...        
    end;

            
%{
d   f  Q    z
1   2  -5   -10
0   5  -5   -25
-1  2  -5   -10
    
    
    
    
%}
    
% ***************************************************************
% HEREBY IS REPRESENTED THE OLD CODE FOR TEST OF OPTIONS STRATEGY
% ***************************************************************
%{
global OM
global p
global price
global i
global days

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
  
    Z = normpdf(y, price(i), price(i)*0.3036*sqrt((days-(i+1)/8)/365));
    AM = trapz(Z.*port)*step;

%f = exp( (lowbound - AM) / 1000 ) + exp( -(lowbound - AM) / 1000 ) -2;
f = AM;%(exp( (lowbound - AM) / 1000 ) + exp( -(lowbound - AM) / 1000 ))*...
    %(exp( (lowbound - max(port)) / 1000 ) + exp( -(lowbound - max(port)) / 1000 ));
%}
f = x*OM(:,3);
%}