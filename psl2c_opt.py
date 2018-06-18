# Program to generate the image of the limit set of a Kleinian group (specifically a Schottky group) with two generators
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import cmath as cm
import time

#values for zero and infinity
zero = 1e-6
infty = 1e6

# list to hold all reduced words of a given max length
charlst = []


# list to hold all points to be plotted
limpts = []


# Designation of matrix values so that tr(abAB) = -2 => limit set is a quasi-circle
y = .4
x = np.sqrt(1+y**2)
v = 2.0
u = np.sqrt(1+v**2)
k = 1/(y*v) + cm.sqrt((1/(y**2 * v**2))-1)

print k

a = np.matrix([[u+0j, 0+(k*v)*1j], [0-(v/k)*1j, u+0j]])
b = np.matrix([[x+0j, y+0j], [y+0j, x+0j]])
A = np.linalg.inv(a)
B = np.linalg.inv(b)
'''
If tr(abAB) = -2, circles are tangent to eachother, connected limit set.
abAB has to be a parabolic element
'''

'''
Overview of function: 
    Recursive function to generate all reduced words with chars in {a, A, b, B} of length 'max' and store them in lst.
    Written to be able to be called with an empty word, if desired.
    Also carries out the associated matrix multiplication and storing of limit points to be plotted, all inside of the 
    helper function 'rescheck'. 
Description of parameters: 
    'word' is the string to be added onto, empty or consisting of chars in {a, A, b, B}. 
    'letter' is the last letter of 'word', i.e. word[-1]. 
    'max' is our desired length for our final words to be stored in lst, sufficiently large to create accurate images.
    'mat' is the identity, or some composition of maps a, A, b, and B.
'''
def genwords(word, letter, max, mat):
    point = 1 + 1j
    x3 = 1
    if max > len(word):
        # case avoiding appending a's inverse, A
        if letter == 'a':
            newletter = 'a'
            newword = word + 'a'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, a, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'b'
            newword = word + 'b'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, b, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'B'
            newword = word + 'B'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, B, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

        # case avoiding appending b's inverse, B
        elif letter == 'b':
            newletter = 'a'
            newword = word + 'a'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, a, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'b'
            newword = word + 'b'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, b, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'A'
            newword = word + 'A'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, A, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

        # case avoiding appending A's inverse, a
        elif letter == 'A':
            newletter = 'B'
            newword = word + 'B'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, B, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'b'
            newword = word + 'b'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, b, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'A'
            newword = word + 'A'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, A, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

        # case avoiding appending B's inverse, b
        elif letter == 'B':
            newletter = 'a'
            newword = word + 'a'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, a, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'B'
            newword = word + 'B'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, B, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
                pass
            else:
                genwords(newword, newletter, max, truthmat[1])

            newletter = 'A'
            newword = word + 'A'

            # list of the form [bool, matrix]
            truthmat = rescheck(mat, point, A, x3)
            if truthmat[0] == 0 or truthmat[0] == 2:
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


'''
Overview of function: 
    Plots the limit points of the Kleinian group generated by a and b by applying elements of matlst to a point in the
    upper half plane. Some products, z, will be near the real axis, others will show signs of diverging. We consider 
    only those z near the real axis, and approximate the true limit point by z.real. 
'''
def plotlimpts():
    reals = np.real(limpts)
    imags = np.imag(limpts)
    plt.plot(reals, imags, 'bo', markersize=0.1)
   # plt.ylim(ymax=1.2)
   # plt.ylim(ymin=-1.2)
    plt.axes().set_aspect('equal', 'datalim')
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
        first = f2(point, s, x3)
        second = f1(first[0], t, first[1])
        third = f3(second[0], second[1])
        fourth = f1(third[0], var1, third[1])
        return f2(fourth[0], var2, fourth[1])
    else :
        first = f2(point, q, x3)
        second = f1(first[0], r, first[1])
        return f2(second[0], 1/t, second[1])
    


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
    return [np.conj(point) / norm, x3 / norm]


'''
Overview of function:
    Helper function: find the fixed points of a Mobius transformation
'''
def fixpt(mat):
    det = np.linalg.det(mat)
    sqrtdet = np.sqrt(det)
    trace = np.trace(mat)
    if abs(mat[1, 0]) >= zero:
        if abs(trace-2*sqrtdet) <= zero or abs(trace+2*sqrtdet) <= zero:
            sqrtdiscr = np.sqrt(trace**2-4*det) 
	    p1 = (mat[0, 0]-mat[1, 1] + sqrtdiscr)/(2*mat[1, 0])
	    p2 = p1
	    return [p1 ,p2]
	else: 
	    sqrtdiscr = np.sqrt(trace**2-4*det) 
	    p1 = (mat[0, 0]-mat[1, 1] + sqrtdiscr)/(2*mat[1, 0])
	    p2 = p1 - sqrtdiscr/(mat[1, 0])
	    return [p1 ,p2]
    elif abs(mat[0, 0]-mat[1, 1]) >= zero:
	if mat[0, 0] < mat[1, 1]:
	    p1 = mat[0, 1]/(mat[1, 1]-mat[0, 0])
	    p2 = infty
	    return [p1,p2]
	else:
	    p1 = infty
	    p2 = mat[0, 1]/(mat[1, 1]-mat[0, 0])
	    return [p1,p2]
    else:
	p1 = infty
	p2 = infty
	return [p1,p2]


'''
Overview of function: 
     Checks the precision of the Mobius transforms in order that excessive matrix multiplication and word generation
     be avoided. Returns a number, used to indicate the proper course of action, and the matrix that is the result of 
     mat * trans. 
Description of parameters: 
    'mat' is the matrix resulting from multiplaction of Mobius transformations hitherto. 
    'point' and 'x3' together make a point (x1,x2,x3) [with point = (x1, x2)] in the upper half space on which the 
    transformations are acting. 
    'trans' is the next transformation to be multiplied. 
'''
def rescheck(mat, point, trans, x3):
    y = comp(point,  mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1], x3)
    mat = np.matmul(mat, trans)
    x = comp(point,  mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1], x3)

    if x[1] <= 0.0005:
        if (x[1]-y[1]) < 0:
            limpts.append(x[0])
            return [0, mat]
        else:
            return [1, mat]
    elif x[1] >= 5:
        return [2, mat]
    else:
        return [1, mat]


'''
MAIN

'''

t0=time.time()
max = 100
matrix1 = np.matrix([[1+0j, 0+0j], [0+0j, 1+0j]])
genwords('', '', max, matrix1)
t1=time.time()
total=t1-t0
print ('num limit points graphed = ')
print (len(limpts))

print ('versus generating all: ')
print (4*3**(max-1))
print ('total time =')
print (total)
print ('fixed point of the commutator word abAB:')
print (fixpt(np.matmul(np.matmul(a,b),np.matmul(A,B))))
print ('fixed point of the commutator word bABa:')
print (fixpt(np.matmul(np.matmul(b,A),np.matmul(B,a))))
print ('fixed point of the commutator word ABab:')
print (fixpt(np.matmul(np.matmul(A,B),np.matmul(a,b))))
print ('fixed point of the commutator word BabA:')
print (fixpt(np.matmul(np.matmul(B,a),np.matmul(b,A))))


plotlimpts()
