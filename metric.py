#tidy this up
from sympy import  *


# This line just sets up pretty printing for later.

init_printing(use_latex='mathjax')

# These are the most commently used coordinate symbols. This should be extended.

t, x, y, z, c, G, M = symbols('t x y z c G M')
r, theta, phi = symbols('r theta phi')
lamda = Symbol("lamda") #note the spelling!

# Here is where the heavy lifting is done. Write a metric class which contains the names of the coordinates, 
# a convenience variable storing the number of coordinates, and an array with the metric elements. 

# There are convenience functions to return the dual, derivatives of elements, connection coefficients, and so on.

# It might make sense to tabulate and store the connections once on instantiation.

class Metric():
    def __init__(self,coordinates):
        #initially set up with 1's on diagonal
        self.coordinates = coordinates
        self.dimension = len(coordinates)
        self.array = MutableDenseNDimArray([0] * self.dimension * self.dimension, [self.dimension, self.dimension])
        for i in range(self.dimension):
            self.array[i,i]=1
            
    def rank(self):
        return self.array.rank()
    
    def setDiagonal(self, diag):
        if len(diag) == self.dimension:
            for i in range(self.dimension):
                self.array[i,i] = diag[i]
        else:
            print ('Incorrect number of elements')
        
    def dual(self):
        return self.array.tomatrix().inv() 
    
    def d(self, i, j, k):
        #for convenience - derivative wrt k of i,j'th element
        return diff(self.array[i,j], self.coordinates[k])
    
    def connection(self, m, i, j):
        #indices to match MathWorld
        sum  = 0
        dual = self.dual()
        for k in range(self.dimension):
            sum += dual[k, m] * (self.d(i,k,j)+ self.d(j,k,i) - self.d(i,j,k))
        return simplify(sum / 2)
    
    def riemann(self, i, j, k, l):
        sum = 0
        sum += diff (self.connection(i,j,l), self.coordinates[k])
        sum -= diff (self.connection(i,j,k), self.coordinates[l])
        for m in range(self.dimension):
            sum += self.connection(m,j,l)*self.connection(i,m,k)
            sum -= self.connection(m,j,k)*self.connection(i,m,l)
        return trigsimp(expand_trig(simplify(sum)))
    
    def ricci(self):
        ricciArray = MutableDenseNDimArray([0] * self.dimension * self.dimension, [self.dimension, self.dimension])
        for i in range(2):
            for j in range(2):
                sum = 0
                for m in range(self.dimension):
                    # note unusual contraction used in RGC 
                    sum += self.riemann(m, i, j, m)
                ricciArray[i, j] = sum
        return ricciArray
    
    def ricciScalar(self):
        ricciArray = self.ricci()
        sum = 0
        dual = self.dual()
        for m in range(self.dimension):
            for n in range(self.dimension):
                sum += dual[m, n] * ricciArray[m, n]
        # usual CAS kick to properly simplify trig expressions
        return trigsimp(expand_trig(simplify(sum)))
        
    def geodesicEquation(self, i):
        myEquation = 0
        if i in range(self.dimension):
            myEquation += Derivative(self.coordinates[i], lamda, lamda)
            for j in range(self.dimension):
                for k in range(self.dimension):
                    myEquation += self.connection (i, j, k) * Derivative (self.coordinates[j], lamda) * Derivative (self.coordinates[k], lamda)
        # put into '= 0' form for clarity
        return Eq(myEquation,0)
    
    def gaussianCurvature (self):
        # lower the first index on the magic Riemann component
        r_0101 = 0
        for  i in range(self.dimension):
            r_0101 += self.array[i,0]*self.riemann(i,1,0,1)
        # divide by the determinant of the metric
        g = self.array.tomatrix().det()
        return r_0101 / g