"""Global variables that can be modified by the user at the IPython command line.
Access programatically using:

get_ipython().user_ns
"""

"""Mean spike rate that delineates normal vs "quiet" neurons. 0.1 Hz seems reasonable if you
plot mean spike rate distributions for all the neurons in a given track. But, if you want a
reasonable looking DJS histogram withouth a lot of missing netstates, you need to exclude
more low firing rate cells, 0.5 works better"""
MINRATE = 0.1 # Hz

"""NeuronCode (Ising matrix) and network state parameters"""
CODEKIND = 'binary'
CODETRES = 20000 # us
CODEPHASE = 0 # deg
CODEWORDLEN = 10 # in bits
