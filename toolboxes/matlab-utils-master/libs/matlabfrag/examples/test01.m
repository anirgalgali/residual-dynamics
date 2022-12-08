hfig = figure;
set(hfig,'units','centimeters','NumberTitle','off','Name','test01');
pos = get(hfig,'position');
set(hfig,'position',[pos(1:2),8,6]);

%%
plot([0,1],rand(6,2));
ylim([0,1]);
hleg=legend('line 1','line 2','3','4','5','6');
h6=findobj(hleg,'string','6');
set(h6,'userdata','matlabfrag:Line6');
%%
matlabfrag('graphics/test01');