# Program to generate the image of the limit set of a Kleinian group (specifically a Schottky group) with two generators
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


# list to hold all reduced words of a given max length
charlst = []

# list to hold all resultant matrices of the reduced words in charlst
matlst = []

limpts = []


# Matrix assignment
a = np.matrix([[9+0j, 0+0j], [0+0j, 1+0j]])
b = np.matrix([[5+0j, 4+0j], [4+0j, 5+0j]])
A = np.linalg.inv(a)
B = np.linalg.inv(b)


'''
Overview of function: 
    Recursive function to generate all reduced words with chars in {a, A, b, B} of length 'max' and store them in lst.
    Written to be able to be called with an empty word, if desired.
Description of parameters: 
    'word' is the string to be added onto, empty or consisting of chars in {a, A, b, B}. 
    'letter' is the last letter of 'word', i.e. word[-1]. 
    'max' is our desired length for our final words to be stored in lst, sufficiently large to create accurate images.
    'mat' is the identity, or some composition of maps a, A, b, and B.
'''
def genwords(word, letter, max, mat):
    point = 1 + 1j
    if max > len(word):
        # case avoiding appending a's inverse, A
        if letter == 'a':
            newletter = 'a'
            newword = word + 'a'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, a)
            if truthmat[0]:
                pass
            #what should happen if !truthmat[0]?
            #that means that a. it isnt getting closer,
            #b. it is above some bound or
            #c. neither of those things but its not within our determined proximity
            #if it is a bad point we do not want to continue to calculate
            #if it is unknown we want to keep calculating
            #if it is a good point we want to graph(? or maybe keep calculating, find more words that begin like that?)
            #if we get that a 7 letter word is close enough to the boundary, do we care about the max-length
            #words beginning with those 7 letters? probably not right? just worried maybe some points will degenerate
            #and we'll miss out on points in the limit set. thinking if we know a shorter word is good
            #(close and getting closer) we can graph all(? or some) max-length words that start with the shorter word
            #I only bring this up because if we change the "close enough" bound we graph a different number of points
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'b'
            newword = word + 'b'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, b)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'B'
            newword = word + 'B'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, B)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

        # case avoiding appending b's inverse, B
        elif letter == 'b':
            newletter = 'a'
            newword = word + 'a'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, a)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'b'
            newword = word + 'b'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, b)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'A'
            newword = word + 'A'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, A)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

        # case avoiding appending A's inverse, a
        elif letter == 'A':
            newletter = 'B'
            newword = word + 'B'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, B)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'b'
            newword = word + 'b'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, b)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'A'
            newword = word + 'A'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, A)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

        # case avoiding appending B's inverse, b
        elif letter == 'B':
            newletter = 'a'
            newword = word + 'a'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, a)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'B'
            newword = word + 'B'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, B)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'A'
            newword = word + 'A'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, A)
            if truthmat[0]:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

        # case with empty (word == ''), (letter == ''), and (for the fctn to work) (mat = ID).
        else:
            genwords('a', 'a', max, mat)
            genwords('A', 'A', max, mat)
            genwords('b', 'b', max, mat)
            genwords('B', 'B', max, mat)
    else:
        charlst.append(word)

        # FIXME for testing
        # print charlst

        '''
        zarray = []
        for pt in limpts:
            zarray.append(0)
        plt.plot(limpts, zarray, 'bo', markersize=0.1)
        plt.show()
        '''



'''
Overview of function: 
    Given Mobius transformations in PSL(2,R) a and b (and inverses A and B) translates words in charlst to
    resultant matrices and stores the results in matlst. 
    Note: matrices are multiplied on the left, so if 'aBAA' is an element of charlst, the matrix AABa will be stored in 
    matlst.
'''
def multwords():
    for x in charlst:
        result = np.matrix([[1 + 0j, 0 + 0j], [0 + 0j, 1 + 0j]])
        for y in reversed(x):
            if y == 'a':
                result = np.matmul(a, result)
            if y == 'b':
                result = np.matmul(b, result)
            if y == 'A':
                result = np.matmul(A, result)
            if y == 'B':
                result = np.matmul(B, result)
        matlst.append(result)



'''
Overview of function: 
    Plots the limit points of the Kleinian group generated by a and b by applying elements of matlst to a point in the
    upper half plane. Some products, z, will be near the real axis, others will show signs of diverging. We consider 
    only those z near the real axis, and approximate the true limit point by z.real. 
'''


def plotlimpts():
    '''
    point = 1 + 1j
    limpts = []
    '''
    zarray = []
    #for mat in matlst:
        #x = (mat[0,0]*point + mat[0,1])/(mat[1,0]*point + mat[1,1])
        #if np.imag(x) <= 0.1:
         #   limpts.append(x)
    for pt in limpts:
        zarray.append(0)

    plt.plot(limpts, zarray, 'bo', markersize=0.5)
    plt.show()
    return


# Concerned that mat is not changed when this funciton returns in the recursion in genwords

def rescheck(mat, point, trans):
    y = (mat[0, 0] * point + mat[0, 1]) / (mat[1, 0] * point + mat[1, 1])
    mat = np.matmul(trans, mat)
    x = (mat[0, 0] * point + mat[0, 1]) / (mat[1, 0] * point + mat[1, 1])


    # also check if imag is decreasing once bound is checked. nested conditional
    if np.imag(x) <= 0.01:
        if np.imag(x-y) <= 0:
            limpts.append(np.real(x))
        # plt.plot(np.real(x), 0, 'bo', markersize=0.1)
            return [True, mat]
        else:
            return [False, mat]
    # what is a good upper bound to stop considering? or should we check decreasing condition vs hard bound?
    # tried to make it check if it gets closer each step
    # mentors said to throw out increasing condition
    # elif np.imag(x-y) >= 0:
        # return [False, mat]
    elif np.imag(x) >= 50:
        return [False, mat]
    else:
        return [False, mat]

'''
MAIN

'''
max = 20
matrix1 = np.matrix([[1+0j, 0+0j], [0+0j, 1+0j]])
genwords('', '', max, matrix1)
print ('num limit points graphed = ')
print (len(limpts))

print ('versus generating all: ')
print (4*3**(max-1))

plotlimpts()






