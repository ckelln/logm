# Program to generate the image of the limit set of a Kleinian group (specifically a Schottky group) with two generators
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


# list to hold all reduced words of a given max length
charlst = []

# list to hold all resultant matrices of the reduced words in charlst
matlst = []

# matrices to be overwritten later as Mobius transformations

a = np.matrix([[1+0j, 0+0j],[0+0j, 1+0j]])
b = np.matrix([[1+0j, 0+0j],[0+0j, 1+0j]])
A = np.matrix([[1+0j, 0+0j],[0+0j, 1+0j]])
B = np.matrix([[1+0j, 0+0j],[0+0j, 1+0j]])


'''
Overview of function: 
    Recursive function to generate all reduced words with chars in {a, A, b, B} of length 'max' and store them in lst.
    Written to be able to be called with an empty word, if desired.
Description of parameters: 
    'word' is the string to be added onto, empty or consisting of chars in {a, A, b, B}. 
    'letter' is the last letter of 'word', i.e. word[-1]. 
    'max' is our desired length for our final words to be stored in lst, sufficiently large to create accurate images.
'''
def genwords(word, letter, max):
    if max > len(word):
        # case avoiding appending a's inverse, A
        if letter == 'a':
            newletter = 'a'
            newword = word + 'a'
            genwords(newword, newletter, max)
            newletter = 'b'
            newword = word + 'b'
            genwords(newword, newletter, max)
            newletter = 'B'
            newword = word + 'B'
            genwords(newword, newletter, max)
        # case avoiding appending b's inverse, B
        elif letter == 'b':
            newletter = 'a'
            newword = word + 'a'
            genwords(newword, newletter, max)
            newletter = 'b'
            newword = word + 'b'
            genwords(newword, newletter, max)
            newletter = 'A'
            newword = word + 'A'
            genwords(newword, newletter, max)
         # case avoiding appending A's inverse, a
        elif letter == 'A':
            newletter = 'B'
            newword = word + 'B'
            genwords(newword, newletter, max)
            newletter = 'b'
            newword = word + 'b'
            genwords(newword, newletter, max)
            newletter = 'A'
            newword = word + 'A'
            genwords(newword, newletter, max)
        # case avoiding appending B's inverse, b
        elif letter == 'B':
            newletter = 'a'
            newword = word + 'a'
            genwords(newword, newletter, max)
            newletter = 'B'
            newword = word + 'B'
            genwords(newword, newletter, max)
            newletter = 'A'
            newword = word + 'A'
            genwords(newword, newletter, max)
        # case with empty (word == '') and (letter == '').
        else:
            genwords('a', 'a', max)
            genwords('A', 'A', max)
            genwords('b', 'b', max)
            genwords('B', 'B', max)
    else:
        charlst.append(word)


'''    
Overview of function: 
    Assigns Mobius transformations in PSL(2, R) to a, b, A, and B. 
'''
def pickmats():
    global a
    a = np.matrix([[9+0j, 0+0j], [0+0j, 1+0j]])
    global b
    b = np.matrix([[5+0j, 4+0j], [4+0j, 5+0j]])
    global A
    A = np.linalg.inv(a)
    global B
    B = np.linalg.inv(b)


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
    initpoint = 1 + 1j
    x3 = 20
    limpts = []
    for mat in matlst:
        # composition of elementary transforms to get action on the Reimann sphere extended to upper half space
        # endpoint is the result of extended action in the form of [x1 + ix2, x3]
        endpoint = comp(initpoint, mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1], x3)
        if endpoint[1] <= 0.1:
            limpts.append(endpoint[0])
    reals = np.real(limpts)
    imags = np.imag(limpts)
    plt.plot(reals, imags, 'bo', markersize=0.1)
    plt.show()
    return


'''
Overview of function: 
    Helper function to define composition of elem transforms to get action on the Reimann sphere extended 
    to upper half space. f2(s) \circ f1(t) \circ f3 \circ f1(r-t*(q/s) \circ f2(s) \circ point.
    Given [[q, r],[s,t]] in PSL(2,C), converts to extended Mobius transformation via composition of simple
    transformations.
    Returns a list of the form [point, x3] where point is of the form x1+ix2, a complex number. 
Description of parameters: 
    'point' is a complex number x1+ix2, geometrically thought of as (x1, x2) on the boundary of H^3. 
    'q, r, s, t' define the Mobius transformation, represented in PSL(2,C) as [[q, r],[s,t]].
    'x3' is the z-axis coordinate in H^3, completing the triple (x1, x2, x3).
'''
def comp(point, q, r, s, t, x3):
    var1 = q / (r*s - t*q)
    if s != 0:
        var2 = r - (t*q)/s
    else :
        # FIXME
        var2 = r
    first = f2(point, s, x3)
    second = f1(first[0], t, first[1])
    third = f3(second[0], second[1])
    fourth = f1(third[0], var1, third[1])
    return f2(fourth[0], var2, fourth[1])


'''
Overview of function: 
    Helper subfunction: addition of b and extension to upper half space
'''
def f1(point, r, x3):
    return [point + r, x3]


'''
Overview of function: 
    Helper subfunction: multiplication by q and extension to upper half space
'''
def f2(point, q, x3):
    return [q*point, np.absolute(q)*x3]


'''
Overview of function: 
    Helper subfunction: rotation and extension to upper half space
'''
def f3(point, x3):
    norm = np.absolute(point)**2 + x3**2
    return [np.conj(point) * norm, x3*norm]


'''
MAIN

'''

genwords('', '', 10)
pickmats()
multwords()
plotlimpts()

