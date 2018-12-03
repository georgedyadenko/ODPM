globals

logintimeout(20);
if ~exist('conn','var')
    conn = database('EXPORTBASE', '', '');
end;

% ONE-TABLE IMPORT
%{
fields={'active','class','till_maturity','bidprice','bidvol','askprice','askvol','strike','type','baseactive','tick','tick_value','lot','implied_volativity','fut_class','fut_till_maturity','fut_tick','fut_tick_value','fut_lot','fut_bidprice','fut_bidvol','fut_askprice','fut_askvol','fut_baseactive'};
f = fetch(conn, 'select * from OPTandFUT')';
if ~isempty(f)
    cpt = cell2struct(f, fields);
else
    display('ERROR: "OPTandFUT" SQL query is empty');
end;

clear f
clear fields
%}

% TWO-TABLE IMPORT
fields={'desc','active','class','tomat','tick','tickval','lot','bidprice','bidvol','askprice','askvol','baseactive'};
f = fetch(conn, 'select * from FUTPRICES')';
if ~isempty(f)
    cptf = cell2struct(f, fields);  %Mf
end;

if exist('cptf','var')
    fields={'active','class','tomat','bidprice','bidvol','askprice','askvol','strike','type','baseactive','desc','tick','tickval','lot','implvol'};
    f = fetch(conn, 'select * from OPTIONPRICES')';
    if ~isempty(f)
        cpto = cell2struct(f, fields); %Mo
    end;
else
    display('error: no futures prices');
end;

if exist('cpto','var')
    for i = 1:1:length(cpto)
        for j = 1:1:length(cptf)
            if strcmp(cpto(i).baseactive, cptf(j).active)
                cpto(i).baseactive_index = j;
                break
            end;
        end;
    end;
else
    display('error: no options prices');
end;

clear f
clear fields
clear i j