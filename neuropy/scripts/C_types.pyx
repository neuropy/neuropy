# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: profile=False

from libc.stdio cimport printf

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
    print('long long: %d bytes' % sizeof(long long))
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
    cdef int i = 5
    cdef int j = 2
    printf('5 / 2 = %d\n', (5 / 2))
    printf('5 / <double>2 = %.3f\n', (5 / <double>2))
    printf('5.0 / 2 = %.3f\n', (5.0 / 2))
    printf('<double>5 / 2 = %.3f\n', (<double>5 / 2))
