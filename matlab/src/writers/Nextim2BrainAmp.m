function Nextim2BrainAmp (filename)

[p,f,x] = fileparts(Filename);

[hdr,data_out] = wlb_read_nextim(filename);

hdr.DataFile = strcat(f,'.eeg');

write_brainvision_eeg(p, hdr, data_out);
write_brainvision_vhdr(p, hdr);
end

