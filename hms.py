#!/usr/bin/python -tt

import sys
from sage.all import *
import cProfile
import random


#This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    if '-profile' in sys.argv:
        cProfile.run('main()')
    else:
        main()



def main():
    pass


def datall_test(numsamples=10**4, maxheight=10**6, disc=5, with_w=True):
    '''
    Usually followed by
    
    var('a,b,x')
    model(x) = a*x+b
    find_fit(datall[lower:upper], model)
    scatter_plot(datall,markersize=1) + plot(m*x + k, (x,40,70))
    
    for appropriate choices of lower and upper, where m and k are the outputs of the find_fit.
    '''
    generic = False
    samples = get_random_I10s(disc, numsamples, maxheight, generic=generic, with_w=with_w)
    data = samples_to_counts(samples)
    datall = [(log(tmp[0])/log(10.0), log(RR(tmp[1]))/log(10.0)) for tmp in data if tmp[1] > 0]
    return datall


def samples_to_counts(samples, make_positive=True):
    '''
    Assumed all samples are positive. If not then add some abs().
    '''
    if make_positive:
        samples = [abs(tmp) for tmp in samples]
        #samples = [tmp for tmp in samples if tmp > 0]
        #samples = [-tmp for tmp in samples if tmp < 0]
    samples.sort()
    lower = floor(log(RR(samples[0]))/log(10.0))
    upper = ceil(log(RR(samples[-1]))/log(10.0))
    checkpoints = [10**(i/RR(2)) for i in range(2*lower, 2*upper+1)]
    next_checkpoint_counter = 0
    next_checkpoint = checkpoints[next_checkpoint_counter]
    data = []
    for (j,samp) in enumerate(samples):
        if samp > next_checkpoint:
            data.append((next_checkpoint, j))
            next_checkpoint_counter += 1
            next_checkpoint = checkpoints[next_checkpoint_counter]
    return data

#######################
### random sampling ###
#######################

def get_random_I10s(disc, numsamples, maxheight, generic=False, with_w=False):
    I10_list = []
    for _ in range(numsamples):
        (x,y,z) = get_random_point(maxheight, with_w=with_w)
        if generic:
            I10_list.append(get_I10_generic(x,y,z,disc))
        else:
            I10_list.append(get_I10(x,y,z,disc))
    return I10_list


def get_random_point(maxheight, sageint=True, with_w=True, with_RM5_cond=False):
    if not with_RM5_cond:
        (x,y,z) = get_random_point_naive(maxheight, sageint=sageint, with_w=with_w)
    else:
        (x,y,z) = get_random_point_RM5_cond(maxheight)
    return (x,y,z)


def get_random_point_naive(maxheight, sageint=True, with_w=False):
    x = random.randint(-maxheight, maxheight)
    y = random.randint(-maxheight, maxheight)
    z = 0
    while z == 0:
        z = random.randint(-maxheight, maxheight)
    if with_w:
        w = 0
        while w == 0:
            w = random.randint(-maxheight, maxheight)
        x = ZZ(x)/ZZ(z)
        y = ZZ(y)/ZZ(w)
        z = ZZ(1)
    if sageint and not with_w:
        x = ZZ(x)
        y = ZZ(y)
        z = ZZ(z)
    return (x,y,z)


def get_random_point_RM5_cond(maxheight, enforce_min=False):
    alpha = 0.5#random.random()*0.9 + 0.05
    beta = 0.5#random.random()*0.9 + 0.05
    if alpha < beta:
        tmp = beta
        beta = alpha
        alpha = tmp
    delta = 1 - alpha
    gamma = 1 - beta
    maxheight_x = int(maxheight**alpha)
    maxheight_y = int(maxheight**gamma) # !
    maxheight_z = int(maxheight**beta)
    maxheight_w = int(maxheight**delta)
    x = random.randint(-maxheight_x, maxheight_x)
    while enforce_min and (abs(x) < maxheight_x/2):
        x = random.randint(-maxheight_x, maxheight_x)
    y = random.randint(-maxheight_y, maxheight_y)
    while enforce_mind and (abs(y) < maxheight_y/2):
        y = random.randint(-maxheight_y, maxheight_y)
    z = random.randint(-maxheight_z, maxheight_z)
    while enforce_min and (abs(z) < maxheight_z/2):
        z = random.randint(-maxheight_z, maxheight_z)
    w = random.randint(-maxheight_w, maxheight_w)
    while enforce_min and (abs(w) < maxheight_w/2):
        w = random.randint(-maxheight_w, maxheight_w)
    x = ZZ(x)/ZZ(z)
    y = ZZ(y)/ZZ(w)
    z = ZZ(1)
    return (x,y,z)


####################################
### get Igusa-Clebsch invariants ###
####################################

def get_I10_generic(x,y,z,disc):
    if disc == 5:
        return RM5_I10_generic(x,y,z)
    if (disc == 8) or (disc == 2):
        return RM8_I10_generic(x,y,z)
    raise NotImplementedError


def get_I10(x,y,z,disc):
    if disc == 5:
        return RM5_IC_invs(x,y,z, rescale=True)[3]
    if (disc == 8) or (disc == 2):
        return RM8_IC_invs(x,y,z, rescale=True)[3]
    raise NotImplementedError


def RM5_IC_invs(x,y,z, rescale=True):
    I2 = RM5_I2(x,y,z)
    I4 = RM5_I4(x,y,z)
    I6 = RM5_I6(x,y,z)
    I10 = RM5_I10(x,y,z)
    IC_invs = (I2,I4,I6,I10)
    if rescale:
        IC_invs = get_rescaling(IC_invs)
    return IC_invs


def RM8_IC_invs(x,y,z, rescale=True):
    I2 = RM8_I2(x,y,z)
    I4 = RM8_I4(x,y,z)
    I6 = RM8_I6(x,y,z)
    I10 = RM8_I10(x,y,z)
    IC_invs = (I2,I4,I6,I10)
    if rescale:
        IC_invs = get_rescaling(IC_invs)
    return IC_invs


def get_rescaling(IC_invs):
    '''
    IC_invs is assumed to be a tuple (I2,I4,I6,I10) of rational numbers which is taken to be an element of weighted projective space with weights 1,2,3,5
    This returns the tuple (a*I2, a**2*I4, a**3*I6, a**5*I10) where a is the smallest positive integer such that all the entries of the returned tuple are integral
    '''
    (I2,I4,I6,I10) = IC_invs

    I2d = I2.denominator()
    I4d = I4.denominator()
    I6d = I6.denominator()
    I10d = I10.denominator()

    a2d = I2d
    a4d = 1
    for (fac, mult) in list(I4d.factor()):
        a4d *= fac**(ceil(ZZ(mult)/2))
    a6d = 1
    for (fac, mult) in list(I6d.factor()):
        a6d *= fac**(ceil(ZZ(mult)/3))
    a10d = 1
    for (fac, mult) in list(I10d.factor()):
        a10d *= fac**(ceil(ZZ(mult)/5))
    a = lcm([a2d,a4d,a6d,a10d])

    to_ret = (a*I2, a**2*I4, a**3*I6, a**5*I10)
    return to_ret


####################
### Generic I10s ###
####################

def RM5_I10_generic(x,y,z):
    to_ret = 8*(x**5 - 10*x**3*y**2 + 25*x*y**4 + 5*x**4*z - 50*x**2*y**2*z + 125*y**4*z - 5*x**3*z**2 + 25*x*y**2*z**2 - 45*x**2*z**3 + 225*y**2*z**3 + 108*z**5)**2
    return to_ret


def RM8_I10_generic(x,y,z):
    to_ret = -8 * (x+z)**6 * (x-z)**3 * z**3 * (32*y**4 + x**3*z - 72*y**2*z**2 + 27*z**4 - 16*x**2*y**2 + 9*x**2*z**2 - 56*x*y**2*z + 27*x*z**3)**2
    return to_ret


###########
### RM8 ###
###########

def RM8_I2(x,y,z):
    num = 96*y**4 - 5*x**3*z - 152*y**2*z**2 + 57*z**4 - (48*y**2 - 19*z**2)*x**2 - (104*y**2*z - 57*z**3)*x
    denom = 4*(2*y**2 - z**2)*(x + z)*z
    to_ret = num/denom
    return to_ret


def RM8_I4(x,y,z):
    num = 144*x**3*y**2 - 288*x*y**4 - 5*x**4*z + 488*x**2*y**2*z + 1312*y**4*z - 72*x**3*z**2 + 144*x*y**2*z**2 - 234*x**2*z**3 - 1800*y**2*z**3 + 567*z**5
    denom = 64*(2*y**2 - z**2)**2*z
    to_ret = num/denom
    return to_ret


def RM8_I6(x,y,z):
    num = 83968*y**8*z + 5*x**7*z**2 - 245248*y**6*z**3 + 248672*y**4*z**5 - 108360*y**2*z**7 + 17415*z**9 - (176*y**2*z - 93*z**3)*x**6 - (4608*y**4 - 2792*y**2*z**2 + 259*z**4)*x**5 - (18336*y**4*z - 17048*y**2*z**3 + 3955*z**5)*x**4 + (18432*y**6 - 43520*y**4*z**2 + 33328*y**2*z**4 - 8073*z**6)*x**3 - (4608*y**6*z - 18496*y**4*z**3 + 11936*y**2*z**5 - 1935*z**7)*x**2 - (18432*y**8 + 63488*y**6*z**2 - 143360*y**4*z**4 + 88344*y**2*z**6 - 17415*z**8)*x
    denom = 512*(2*y**2 - z**2)**3*(x + z)*z**2
    to_ret = num/denom
    return to_ret


def RM8_I10(x,y,z):
    num = -(32*y**4 + x**3*z - 72*y**2*z**2 + 27*z**4 - (16*y**2 - 9*z**2)*x**2 - (56*y**2*z - 27*z**3)*x)**2*(x + z)*(x - z)**3
    denom = 131072*(2*y**2 - z**2)**5*z**2
    to_ret = num/denom
    return to_ret


###########
### RM5 ###
###########

def RM5_I2(x,y,z):
    num = 4*x**2 - 20*y**2 - 6*z**2
    denom = 5*z**2
    to_ret = num/denom
    return to_ret


def RM5_I4(x,y,z):
    num = (x**2 - 5*y**2 - 9*z**2)**2
    denom = 100*z**4
    to_ret = num/denom
    return to_ret


def RM5_I6(x,y,z):
    num = 75*x**6 - 1125*x**4*y**2 + 5625*x**2*y**4 - 9375*y**6 + 72*x**5*z - 720*x**3*y**2*z + 1800*x*y**4*z - 1165*x**4*z**2 + 11650*x**2*y**2*z**2 - 29125*y**4*z**2 - 360*x**3*z**3 + 1800*x*y**2*z**3 + 5985*x**2*z**4 - 29925*y**2*z**4 - 6399*z**6
    denom = 25000*z**6
    to_ret = num/denom
    return to_ret


def RM5_I10(x,y,z):
    num = (x**5 - 10*x**3*y**2 + 25*x*y**4 + 5*x**4*z - 50*x**2*y**2*z + 125*y**4*z - 5*x**3*z**2 + 25*x*y**2*z**2 - 45*x**2*z**3 + 225*y**2*z**3 + 108*z**5)**2
    denom = 39062500*z**10
    to_ret = num/denom
    return to_ret        
