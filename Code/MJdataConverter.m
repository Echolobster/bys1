%% MJdataConverter
%
%
%%
tsmin = 600;
tsmax = 1000;

slist = {'CV','GB','LB','PG','SM','SWE','TA','VD','VR'};

for i = 1:length(slist)
    load([slist{i} '_BLSbiasedFitResults20170404'],...
        'tsIn','tpIn','lapseTrials','WM','WP','B');
    
    d.tmin = tsmin; %min of the prior in ms
    d.tmax = tsmax; %max of the prior in ms
    d.ts = tsIn{2}(~lapseTrials{2}); %array of all ts used in RSSG
    d.tp = tpIn{2}(~lapseTrials{2}); %array of all tp used in RSSG
    d.wm1 = mean(WM,1); %wm fitted from RSG
    d.wp = mean(WP,1); %wp fitted from RSG
    d.b = mean(B,1); % b fitted from RSG
    
    save(['~/Projects/RS2G_psychophysics/' slist{i} '_RSSG_4_MJ'],...
        '-struct','d')
    
end
    