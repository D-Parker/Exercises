# python3
# Simplex method for optimization - Slow version
# Implementation of CLRS 3rd edition Chapter 29


from sys import stdin
import copy


def create_vectors(A,b,c,n,m):
    
    w= [[0]*(n+m+1) for i in range(n+m+1)]

    for i in range(len(A)):
        for j in range(len(A[0])):    
            w[m+i+1][j+1] = A[i][j]
        
    A = w

    b_ = [0] *(n+m+1)
    
    for i in range(len(b)):
        b_[m+1+i] = b[i]
   
    b = b_

    c = [0] + c + [0]*(n+m-len(c))
    
    return A,b,c


def simplex(A,b,c):
    
    A,b,c = create_vectors(A,b,c,n,m)

    (N,B,A,b,c,v,l,e)=init_simplex(A,b,c)
    
    if e == -1:
        return N,B,A,b,c,v,l,e #no solution
    
    d = [None]*(n+m+1)
    x = True
    while x is True:
        x = False
        #choose entering variable with smallest index
        for j in N:
            if c[j] > 0:
                x = True
                e = j
                break
            
        if x == False:
            break
                    
        for i in B:
            if A[i][e] >0:
                d[i] = b[i]/A[i][e]
            else:
                d[i] = float('Inf')
                                      
        #choose leaving variable with smallest index     
        min_d = float('Inf')
       
        for l in B:       
            min_d = min(min_d, d[l])
    
        for l in B:
            if d[l] == min_d:
                l_min = l
                break
        l = l_min
      
        if d[l] == float('Inf'):
    #            print('yuck')
            e = -2
            return N,B,A,b,c,v,l,e
    #            return N,B,A,b,c,v,l,e #unbounded
        else:
            (N,B,A,b,c,v,l,e) = pivot(N,B,A,b,c,v,l,e)
            
    return N,B,A,b,c,v,l,e
            

def init_simplex(A,b,c):
    
    o_c = copy.deepcopy(c)
    
    k = -1
    
    min_b = float('inf')
    for i in range(m+1,m+n+1):
        min_b = min(b[i],min_b)
        if min_b == b[i]:
            k = i

    N = set()
    for i in range(1,m+1):
        N.add(i)

    B=set()
    for i in range(m+1,n+m+1):
        B.add(i)
    
    l =  k
    v = 0
    
    if b[k] >= 0:
        return N,B,A,b,c,v,l,0
    
    #convert to slack form

    N.add(0)


    for i in range(m+1,n+m+1):
        A[i][0] = -1       
    
    c= [0] *(n+m+1)
    c[0] = -1
    
    l =  k
    v = 0
    
    e = 0
    
    (N,B,A,b,c,v,l,e) = pivot(N,B,A,b,c,v,l,0)
    
    #basic solution is now feasible

    x = True
    d = [None]*(n+m+1)
    while x is True:
        #choose entering variable with smallest index  
        x = False
        for j in N:
            if c[j] > 0:
                x = True
                e = j
                break
        if x == False:
            break
        
        for i in B:
            if A[i][e] >0:
                d[i] = b[i]/(A[i][e])
            else:
                d[i] = float('Inf')

        min_d = float('Inf')
        for l in B:       
            min_d = min(min_d, d[l])
            #choose leaving variable with smallest index  
        for l in B:
            if d[l] == min_d:
                l_min = l
                break
        
        l = l_min

        (N,B,A,b,c,v,l,e) = pivot(N,B,A,b,c,v,l,e)
    
    #optimal solution found at this point, if it exists
    
    feasible = True
    
    if abs(0-v) > .001:
        feasible = False
    
    if feasible == False:
        e = -1
        return N,B,A,b,c,v,l,e        
    
    if 0 in B:
        pivot(N,B,A,b,c,v,0,max(N))
    N.remove(0)
        
        #manipulate into slack form that has a feasible basic solution
        
        #zero out everything with x0 term
        
    for i in range(len(A)):
        A[i][0]=0
        
    for i in range(len(A[0])):
        A[0][i]=0
        
    for i in range(len(c)):
        c[i] = 0
        
    for i in range(len(c)):            
        if i not in B:
            c[i] = o_c[i]
        
    for i in range(len(c)):            
        if i in B:
           for j in range(len(A[i])):
               c[j]= c[j] + o_c[i]*(-A[i][j])
                
    return N,B,A,b,c,v,l,e
    
    
def pivot(N,B,A,b,c,v,l,e):

    b_ =[0] *(n+m+1)
    c_ =[0] *(n+m+1)
    
    a_ = [ ([0]*(n+m+1)) for i in range(n+m+1)]
        
    b_[e] = (b[l])/A[l][e]
    for j in (N-{e}):
        a_[e][j] = (A[l][j])/(A[l][e]) 
    a_[e][l]=1/(A[l][e]) 
    
    #compute coefficients of the remaining constraints  
    for i in (B-{l}):
        b_[i] = b[i] - (A[i][e])*b_[e]
        for j in (N-{e}):
            a_[i][j] = A[i][j] -(A[i][e])*a_[e][j]
        a_[i][l] = -(A[i][e])*a_[e][l]
        
    #compute objective function
    v_ = v + c[e]*b_[e]
    for j in (N-{e}):
        c_[j] = c[j] - (c[e])*a_[e][j]
    c_[l] = -(c[e])*a_[e][l]
    #Compute new sets of basic and nonbasic variables
    
    if e in N:
        N.remove(e)
    N.add(l)
    if l in B:
        B.remove(l)
    B.add(e)
    
    return (N,B,a_,b_,c_,v_,l,e)


def solve_problem(A, b, c):  
  
    (N,B,A,b,c,v,l,e) = simplex(A,b,c)
    
    value = 0
    
    if e == -1:
        value = -1 #no solution
        
    if e == -2:
        value = 1 # infinity
    
    solution = [0] * m
    
    for i in range(1,m+1):
        solution[i-1] = b[i]
    
    return value, solution  
  

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
        
           
    (anst, ansx) = solve_problem(A, b, c)
    
    if anst == -1:
      print("\nNo solution")
    if anst == 0:  
      print("\nBounded solution")
      print(' '.join(list(map(lambda x : '%.18f' % x, ansx))))
    if anst == 1:
      print("Infinity")