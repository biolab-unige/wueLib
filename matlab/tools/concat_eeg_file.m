function [eegHdr, eegData,vmrk_event] = concat_eeg_file(eegFname,eegFname2,ecg_channel)

    [eegHdr, eegData,vmrk_event] = wlb_readBrainvision( eegFname );    
    [eegHdr2, eegData2,vmrk_event2] = wlb_readBrainvision( eegFname2 );
   
    offset = eegHdr.nSamples;
    
    if(not(isempty(vmrk_event2)))
       evt_label2 = unique({vmrk_event2.label});
       
       for iEvtlabl2 = 1:length(evt_label2)
           evt1_ind = find(strcmp({vmrk_event.label},evt_label2{iEvtlabl2}));
           vmrk_event(evt1_ind).type = [vmrk_event(evt1_ind).type; vmrk_event2(iEvtlabl2).type];
           vmrk_event(evt1_ind).epochs = [vmrk_event(evt1_ind).epochs vmrk_event2(iEvtlabl2).epochs];
           vmrk_event(evt1_ind).samples = [vmrk_event(evt1_ind).samples vmrk_event2(iEvtlabl2).samples+offset];
           vmrk_event(evt1_ind).length = [vmrk_event(evt1_ind).length vmrk_event2(iEvtlabl2).length];
           vmrk_event(evt1_ind).chan_num = [vmrk_event(evt1_ind).chan_num vmrk_event2(iEvtlabl2).chan_num];
       end
    end
    
    eegData = cat(2,eegData,eegData2);
    eegHdr.nSamples = size(eegData,2);
    
end