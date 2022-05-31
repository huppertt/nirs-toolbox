function wipe_mat_header(filename)
cmds = ['printf ' '''' ' %.0s' '''' ' {1..97} | dd of=' filename ' bs=1 seek=19 conv=notrunc'];
system(cmds)
end