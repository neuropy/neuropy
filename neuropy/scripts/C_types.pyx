def test():
    """Print out some diagnostic info about how various C types are defined
    on this particular machine"""
    print('char: %d bytes' % sizeof(char))
    print('unsigned char: %d bytes' % sizeof(unsigned char))
    print('short: %d bytes' % sizeof(short))
    print('unsigned short: %d bytes' % sizeof(unsigned short))
    print('int: %d bytes' % sizeof(int))
    print('unsigned int: %d bytes' % sizeof(unsigned int))
    print('long: %d bytes' % sizeof(long))
    print('unsigned long: %d bytes' % sizeof(unsigned long))
    print('float: %d bytes' % sizeof(float))
    print('double: %d bytes' % sizeof(double))
    print('')
    print('char vs unsigned char overflow:')
    cdef char cx = 127
    cx += 1
    print(cx)
    cdef unsigned char ucx = 127
    ucx += 1
    print(ucx)
    
