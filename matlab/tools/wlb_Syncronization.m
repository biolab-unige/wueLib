subj = 'wue10';

pathPcs = ['C:\Users\andrea83\Documents\Projects_data\Data_repository\PC+S-2016\' subj];
pathEmg = ['C:\Users\andrea83\Documents\Projects_data\Data_repository\PC+S-2016\' subj];
pathHdeeg = ['C:\Users\andrea83\Documents\Projects_data\Data_repository\PC+S-2016\' subj];
pathEvent = [''];
outdir = ['C:\Users\andrea83\Documents\Projects_data\Data_repository\PC+S-2016\' subj '\new_hdeeg_pcs_emg_sync'];
close all
status = wlb_TENSSynch_main('pathPcs',pathPcs,'pathHdeeg',pathHdeeg,'pathEmg',pathEmg,'pathEvent',pathEvent,'outdir',outdir,...
    'fNameFilters',{'trial2','walking'},'debug',true,'pcsRefChannel',1,'emgRefChannel','pulse','find_tens_pctg',80,'automaticDetection',true);