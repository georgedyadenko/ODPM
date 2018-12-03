function my_plot(solution)
    globals
    hold off
    
    plot_range = [min(matrix(:,2))-50000 max(matrix(:,2))+50000];
    step = (plot_range(2)-plot_range(1))/200;
    xx = plot_range(1):step:plot_range(2);
    portfolio = zeros(size(xx));
    plot(xx, portfolio, 'r');
    hold on
    
    min_time_tomat = min(matrix(:,9));
    
    for i = 1:1:size(matrix, 1)
        [call_price, put_price] = blsprice(xx, matrix(i,2), 0.001, (matrix(i,9)-min_time_tomat+0.1)/365, 0.3);
        switch matrix(i,3)
            case 1 % if call
                portfolio = portfolio + ...
                    my_port(1,i).*(call_price - my_port(2,i)) ...
                    + my_port(3,i) + ...
                    solution(i).*(call_price - (matrix(i,4) + (solution(i)>0).*matrix(i,5)));
            case 0 % if put
                portfolio = portfolio + ...
                    my_port(1,i).*(put_price - my_port(2,i)) ...
                    + my_port(3,i) + ...
                    solution(i).*(put_price - (matrix(i,4) + (solution(i)>0).*matrix(i,5)));
            case -1
                portfolio = portfolio + ...
                    my_port(1,i).*(xx - my_port(2,i)) ...
                    + my_port(3,i) + ...
                    solution(i).*(xx - (matrix(i,4) + (solution(i)>0).*matrix(i,5)));
        end;
    end;
    
    plot(xx, portfolio);
