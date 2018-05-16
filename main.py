# Program to generate the image of the limit set of a Kleinian group (specifically a Schottky group) with two generators
import numpy as np


# list to hold all reduced words of a given max length
charlst = []

# list to hold all resultant matrices of the reduced words in charlst
matlst = []

# matrices to be overwritten later as Mobius transformations
a = np.matrix([[0+0j, 0+0j],[0+0j, 0+0j]])
b = np.matrix([[0+0j, 0+0j],[0+0j, 0+0j]])
A = np.matrix([[0+0j, 0+0j],[0+0j, 0+0j]])
B = np.matrix([[0+0j, 0+0j],[0+0j, 0+0j]])



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
    a = np.matrix([[9+0j, 0+0j],[0+0j, 1+0j]])
    b = np.matrix([[5+0j, 4+0j],[4+0j, 5+0j]])
    A = np.linalg.inv(a)
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
    return



'''
MAIN


# length to set precision of plot
length = raw_input('Please pick the length for words to be generated')
genwords('','',length)
pickmats()
multwords()
plotlimpts()
'''
