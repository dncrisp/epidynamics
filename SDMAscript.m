% sdma script

sdma1 = SDMAprocess('sdma');


pt_ids = {'NVC1001_25_002'}%,...

%pt_ids = {'NVC1001_25_002'}%,...
%     'NVC1001_24_005',...
%     'NVC1001_23_007',...
%     'NVC1001_24_001',...
%     'NVC1001_25_003_2'};

%  'NVC1001_24_002_2',...
%     'NVC1001_24_004',...
%     'NVC1001_23_005',...
%     'NVC1001_23_002',...

chans = 57;
start = 50430;
stop = 53810;

%chans = [15, 13, 9, 6, 2];
%%


for i=1:1
    
    np = sdma1.addnewpatient(pt_ids{i});
    %np.globalchannel = chans(i);
    np.globalchannel = chans;
    np.downloadszfromieeg(chans,start,stop);

    
     np.analyzeszwithisi;
    
end






% clear all
