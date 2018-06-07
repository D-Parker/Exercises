# python3
# Simplex method for optimization - Faster version

from sys import stdin
import copy
import itertools
#from decimal import *
#getcontext().prec = 20


#n, m = list(map(int, stdin.readline().split()))
#A = []
#for i in range(n):
#  A += [list(map(int, stdin.readline().split()))]
#b = list(map(int, stdin.readline().split()))
#c = list(map(int, stdin.readline().split()))

#n = 3
#m = 2
#
#A = [[-1,-1],[1,0],[0,1]]  # bounded
#b = [-1,2,2]
#c = [-1,2]




def CreateMatrix(A,b,n,m):

    w = [ ([0]*(m)) for i in range(n+1)] # w is coefficient matrix A
    
    y = [0]*(n+m+2)
    
    
    for i in range(len(b)):
        y[i+1] = b[i]
    
    y[-1] = VBN
    
    for i in range(len(A)):
        for j in range(len(A[0])):
            w[i+1][j] = A[i][j]
            
           
    diag = [ ([0]*(m)) for i in range(m+1)]
    for i in range(m):
        for j in range(m):
            if i ==j:
                diag[i][j] = -1
                
                #diagonal is negative <=0 because all x- must be positive
                
    for i in range(m,m+1):
        for j in range(len(diag[0])):
            diag[i][j] = 1
    
    w = w + diag
    
    return w,y  # Matrix and vector are now prepared

# Generate list of all n+m permutation squares of size m
# Solve each m, returning b
# Test that b satisfies all of the inequalities
# Break if an inequality is violated
# Calculate the objective value

def ComputeBest(w,y,c):
    
    test = []

    for i in range(1,len(w)):
        test.append(i)

    perm_count = set(itertools.combinations(test, m))

    perm_list=[0]

    for i in perm_count:
        perm_list.append(list(i))
        
    sol_max_id = -1
    sol_max = -float('inf')
    plug_max =-float('inf')

    for i in range(1,len(perm_list)):
        
        w1 = copy.deepcopy(w)
        y1 = copy.deepcopy(y)
        
        test_a = []
        test_b = []
        
        iteration = i
        
        for j in perm_list[i]:
            test_a.append(w1[j])
            test_b.append(y1[j])
                        
        equation = ReadEquation(test_a,test_b)
        solution = SolveEquation(equation)
        
        #now check solution against all inequalities
        
        for k in range(1,len(test)+1):
                  
            stop = 0
            sol_sum = 0
                                 
            for l in range(len(solution)):
                sol_sum = sol_sum+ solution[l]*w[k][l]
                
            sol_test =  y[k]
                        
            if  (sol_sum <= y[k] + 1e-3)==False: # if violates a constraint. Comparison #1
                stop = 1
                break #no solution
        if stop ==1:
            continue
           
        #solution has one objective value        
        sol_obj = 0
        for i in range(len(c)):
            sol_obj = sol_obj+ solution[i]*c[i]
            
        eq = perm_list[iteration]
        c_test = False
        for i in eq:
            if i == test[-1]:
                c_test=True
                
        if c_test==False:
            if sol_obj > sol_max:  #EPS:
                sol_max_id = iteration 
                sol_max = sol_obj
                max_solution = solution
            
        if c_test==True:
            if sol_obj > plug_max:  #EPS:
                plug_max_id = iteration 
                plug_max = sol_obj
                plug_solution = solution
        
        #find max objective value and associated equations    

    if sol_max_id == -1:
        return -1,0 #no solution
    
    if plug_max > sol_max:
        return -2,0
    
    return 0, max_solution

            
class Equation:
    def __init__(self, a, b):
        self.a = a
        self.b = b


class Position:
    def __init__(self, column, row):
        self.column = column
        self.row = row


def ReadEquation(a,b):

    return Equation(a, b)


def SelectPivotElement(a, used_rows, used_columns):
    # Selects the first free element.
    pivot_element = Position(0, 0)

    while used_columns[pivot_element.column]==True:
        if pivot_element.column >= len(a[0])-1:
            return -1
        pivot_element.column += 1
      
    while used_rows[pivot_element.row]==True or a[pivot_element.row][pivot_element.column]==0:
        if pivot_element.row >= len(a)-1:
            return -1
        pivot_element.row += 1
        
    return pivot_element


def SwapLines(a, b, used_rows, pivot_element):
    a[pivot_element.column], a[pivot_element.row] = a[pivot_element.row], a[pivot_element.column]
    b[pivot_element.column], b[pivot_element.row] = b[pivot_element.row], b[pivot_element.column]
    used_rows[pivot_element.column], used_rows[pivot_element.row] = used_rows[pivot_element.row], used_rows[pivot_element.column]
    pivot_element.row = pivot_element.column;


def ProcessPivotElement(a, b, pivot_element):
    p = a[pivot_element.row][pivot_element.column] 
    
    for i in range(len(a[pivot_element.row])):
        a[pivot_element.row][i] = a[pivot_element.row][i]/p
       
    b[pivot_element.row]= b[pivot_element.row]/p
    
    #test section
    a[pivot_element.row][pivot_element.column] = 1
    #
        
    for i in range(len(a)):
        if i != pivot_element.row:
            z = a[i][pivot_element.column]
            for j in range(len(a[0])):
                a[i][j] = a[i][j] - z *a[pivot_element.row][j]
                
                #test#
                if j == pivot_element.column:
                    a[i][j] = 0
                #test#
            
            b[i] = b[i] - z*b[pivot_element.row]

    pass


def MarkPivotElementUsed(pivot_element, used_rows, used_columns):    
    used_rows[pivot_element.row] = True
    used_columns[pivot_element.column] = True


def SolveEquation(equation):
    a = equation.a
    b = equation.b
    size = len(a)
    used_columns = [False] * size
    used_rows = [False] * size
    
    for step in range(size):
        pivot_element = SelectPivotElement(a, used_rows, used_columns)
        if pivot_element == -1:
            return b
        SwapLines(a, b, used_rows, pivot_element)
        ProcessPivotElement(a, b, pivot_element)
        MarkPivotElementUsed(pivot_element, used_rows, used_columns)
        
    for i in range(len(b)):
        if abs(b[i]-0) < 1e-6:
            b[i] == 0
        
        if abs(b[i]-1) < 1e-6:
            b[i] == 1
      
    return b


if __name__ == '__main__':
    print('Enter number of constraints and variables:\nExample:\n\
          3 2 \n')
    n,m = list(map(int, stdin.readline().split())) # number of constraints
   
    A = []
    print('Enter equation coefficients:\n\nExample:\n\
          -1 -1\n\
           1  0 \n\
           0  1 \n\
          corresponds to: \n\
          -x1 -  x2 \n\
           x1 + 0x2 \n\
          0x1 +  x2 \n          ')
    
    for i in range(n):
      A += [list(map(int, stdin.readline().split()))]
      
    print('For Ax=b, enter b coefficients:\n\
       \nExample: -1 2 2 \n')
    b = list(map(int, stdin.readline().split()))
    
    print('Enter objective function coefficients:\n\n\
        Example:\n\
        -1 2 \n\
        corresponds to: \n\
        -x1 + 2x2 \n' )
    c = list(map(int, stdin.readline().split())) # objective function coefficients
    
    #EPS = 1e-15
    VBN = 10e9
    EPS = 1e-6
    PRECISION = 20
    
    (w,y) = CreateMatrix(A,b,n,m)
    (anst,ansx) = ComputeBest(w,y,c)
    
    if anst == -1:
      print("No solution")
    if anst == 0:  
      print("Bounded solution")
      print(' '.join(list(map(lambda x : '%.18f' % x, ansx))))
    if anst == -2:
      print("Infinity")