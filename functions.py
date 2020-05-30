import math as m 

#declare list of functions (all functions have parameteres a (coefficient) and dx (shift in x axis))
def linear(x, dx = 0, ax = 1, a = 1):
    '''
    y = a * (ax * x + dx)
    '''
    return a * (ax*x + dx)


def sin(x, dx = 0, ax = 1, a = 1):
    '''
    y = a * sin(ax * x + dx)
    '''
    return a * m.sin(linear(x,dx,ax))


def exp(x, dx = 0, ax = 1, a = 1):
    '''
    y = 0.1 * a * e ^ (ax * x + dx)
    '''
    return a * m.exp(linear(x,dx,ax)) 


def ln(x, dx = 0, ax = 1, a = 1):
    '''
    y = a * ln(ax * |x| + dx)
    '''
    tx = linear(x,dx,ax)
    if(tx == 0):
        return 0
    if(tx < 0):
        tx = abs(tx)
    return a * m.log(tx)

def sqrt(x, dx = 0, ax = 1, a = 1):
    '''
    y = a * sqrt(ax * |x| + dx)
    '''
    tx = linear(x,dx,ax)
    if(tx < 0):
        tx = abs(tx)
    return a * m.sqrt(tx)


def my_abs(x, dx = 0, ax = 1, a = 1):
    '''
    y = a * |ax * x + dx|
    '''
    return a * abs(linear(x,dx,ax))


#Если убрать деление на 10, то почему-то 
# все ГА ломаются
def sqr(x, dx = 0, ax = 1, a = 1):
    '''
    y = a * (ax * x + dx) ^ 2
    '''
    return a * linear(x/10,dx,ax)**2


def hyper(x, dx = 0, ax = 1, a = 1):
    '''
    y = a * (ax * x + dx) ^ 3
    '''
    return a * linear(x/10,dx,ax)**3

