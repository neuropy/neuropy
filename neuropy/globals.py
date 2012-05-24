# mean spike rate thresh (Hz) that delineates normal vs "quiet" neurons. 0.1 Hz seems
# reasonable if you plot mean spike rate distributions for all the neurons in a given
# track. But, if you want a reasonable looking DJS histogram withouth a lot of missing
# netstates, you need to exclude more low firing rate cells, 0.5 works better:
QUIETMEANRATETHRESH = 0.1
