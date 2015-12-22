function [hdr, data,event] = wlb_readBrainvision( eegFname )

[hdr] = read_brainvision_vhdr(eegFname);
[data] = read_brainvision_eeg(eegFname, hdr, 1, hdr.nSamples);
[event] = read_brainvision_vmrk(eegFname);


