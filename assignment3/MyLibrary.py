"""
Author: Swaroop Ramakant Avarsekar (2011172)
"""

def lud(mat):
    for i in range(len(mat)):
        for j in range(len(mat)):
            if i>0 and i<=j: # transforming upper triangular matrix
                prod=0
                for k in range(i):
                    prod+=mat[i][k]*mat[k][j]
                mat[i][j]=mat[i][j]-prod #formula
            if i>j: #transforming lower triangular matrix
                prod=0
                for k in range(j):
                    prod+=mat[i][k]*mat[k][j]
                mat[i][j]=(mat[i][j]-prod)/mat[j][j]
    return mat

def fsub(mat,b): #forward substition
    y=[0 for i in range(len(mat))]
    y[0]=b[0][0]
    for i in range(len(mat)): #formula
        prod=0
        for j in range(i):
            prod+=mat[i][j]*y[j]
        y[i]=b[i][0]-prod
    return y

def gauss_seidel(a,b):
    guess=[[0],[0],[0],[0],[0],[0]] # Let this be the guess solution
    prev=[] # let prev list be filled with zeros
    for row in range(len(a)):
        somelist=[]
        for column in range(len(a[0])):
            somelist.append(0)
        prev.append(somelist)
        
    count=0
    new=guess
    for somenum in range(10000): # setting up iteration limit          
        for i in range(len(a)):
            prev[i][0]=new[i][0] # assigning prev to new
    
        for i in range(len(a)): # Formula 
            temp=b[i][0]
            for j in range(i):
                temp=temp-a[i][j]*new[j][0]  
            for j in range(i+1,len(a)):
                temp=temp-a[i][j]*prev[j][0]
            new[i][0]=temp/a[i][i]
        count+=1
        
        d=0 #precision check
        for num in range(len(prev)):
            d+=abs(new[num][0]-prev[num][0])
        if d<0.000001:
            break
    print("Iterations=",count)
    print('Solution=',new)

def cholesky(mat):
    import math
    if mat==transpose(mat): # checking for symmetric matrix
        for i in range(len(mat)):
            add=0
            #print(add)
            for j in range(i):
                add+=(mat[i][j])**2
            mat[i][i]=math.sqrt(mat[i][i]-add)

            for j in range(i+1,len(mat)):
                tot=0
                for k in range(i):
                    tot+=mat[i][k]*mat[k][j]
                mat[j][i]=(mat[j][i]-tot)/mat[i][i]
                mat[i][j] = mat[j][i]

        for i in range(len(mat)): # triangular part zero
            for j in range(len(mat[0])):
                if j>i:
                    mat[i][j]=0
        return mat
    
    else:
        print('Asymmetric matrix,cannot be Cholesky decomposed')

def chol_fsub(mat,coeff):
    y=[]
    for i in range(0,len(mat)):#i=1
        prod=0
        for j in range(i):#for j in range(0,0)
            prod+=(mat[i][j]*y[j])
        y.append((coeff[i][0]-prod)/mat[i][i])
    return y

def chol_bsub(mat1,y):
    x=[0]*len(mat1)
    #print(x)
    for i in range(len(mat1)-1,-1,-1): #decrement
        prod=0
        for j in range(i+1,len(mat1)):
            prod+=mat1[i][j]*x[j]
        #print(i)
        x[i]=((y[i]-prod)/mat1[i][i])
    return x

def transpose(mat):
    mat1 = [[mat[j][i] for j in range(len(mat))] for i in range(len(mat[0]))]
    return mat1

def patrix(m):
    for i in range(len(m)):
        print(m[i])
'''****************'''
def prng(x,c,n): 
    n=int(n) # n is number of elements
    numlist=[] # x is the seed
    for i in range(n):
        numlist.append(c*x*(1-x))
        x=numlist[i] 
    return numlist

def lcg(a,c,m,x,n): 
    n=int(n)
    numlist=[]
    for i in range(n):
        numlist.append(((a*x+c)%m)) # storing rand no.s in list
        x=numlist[i] #changing the input x
    return numlist

def randwalk(a,c,m,x1,x2,n):
    l1=lcg(a,c,m,x1,n) # Lists containing random numbers are l1 and l2
    l2=lcg(a,c,m,x2,n) #x1 and x2 are seeds
    x=[0] #x and y will be the coordinates of the random walk plot
    y=[0]
    d2=0
    for i in range(len(l1)):
        l1[i]=2*(l1[i]/m)-1 # converting rand no.s b/w -1 to 1
        l2[i]=2*(l2[i]/m)-1
        x.append(x[i]+l1[i]) # walking
        y.append(y[i]+l2[i])
        d2+=((x[i+1]-x[i])**2+(y[i+1]-y[i])**2) # summing the squares of the distances
    print('RMS=',(d2/n)**0.5)
    print('net displacement=',((x[-1]-x[0])**2+(y[-1]-y[0])**2)**0.5)
    import matplotlib.pyplot as plt
    plt.plot(x,y) 
    plt.scatter(x[0],y[0],c='b') # initial point(blue)
    plt.scatter(x[-1],y[-1],c='r') # final point(red)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axhline(y=0,c='k',linewidth=1)# line y=0
    plt.axvline(x=0,c='k',linewidth=1)# line x=0
'''************************'''
def oddsum(N):
    N=int(N)
    total=0
    for num in range(2*N): #listing all the first N odd nos 
        if num%2!=0:
            total+=num  # adding the sum
    print('Sum of first', N,'odd terms=',total)
    
def factorial(N):
    N=int(N)
    endnum=1
    for i in range(1,N+1): #zero is not being considered as it will multiply by zero
        endnum=endnum*i  #increment of endnum till the number N
    print(N,'! =',endnum)
    
def apsum(a,d,N): # a is first term, d is common difference, N is no. of terms
    t=0
    N=int(N)
    for i in range(N):
        #print(a)
        t+=a #t is total sum
        a+=d
    print('Sum of first',N,'terms of AP=',t)
    
def gpsum(a,r,N): # a is first term, r is common ratio, N is no. of terms
    t=0
    N=int(N)
    for i in range(N):
        #print(a)
        t+=a #t is total sum
        a=a*r
    print('Sum of first',N,'terms of GP=',t)

def hpsum(a,d,N): # a is first term, d is common difference, N is no. of terms
    t=0
    N=int(N)
    for i in range(N):
        #print((1/a))
        t=t+(a) 
        a=1/(1/(a)+d) 
    print('Sum of first',N,'terms of HP=',t)

def series1(n):
    n=int(n)
    summation=0
    sums,num=[],[] 
    for n in range(1,n+1): # n=0 is not considered
        #summation+=((-1)**(n+1))/(2**n)
        #summation+=float("{:.4f}".format(((-1)**(n+1))/(2**n)))
        summation +=((-1)**(n+1))/(2**n)
        sums.append(summation) #appending sum in list for different n
        num.append(n)
    #print(float("{:.4f}".format(summation)))
    print('Sum of series',round(summation,4))
    import matplotlib.pyplot as plt
    plt.plot(num,sums,markersize=5)
    plt.xlabel('n')
    plt.ylabel('Sum')

def dotprod(c,d):
    dot_p=0
    for i in range(len(c)):
        dot_p+=c[i][0]*d[i][0] # adding after multiplication of the corresponding elements of two matrices
    print(dot_p)
    
def matrixmult(a,b):
    ans=[]
    for m in range (len(a)): # creating a matrix with all entries to be in 0, according to dimension of the result
        row=[]
        for n in range (len(b[0])):
            row.append(0)
        ans.append(row)
    #print(ans)
    for i in range(len(a[0])):   # traversing through column of a
        for j in range(len(b[0])):   # traversing through column of b
            for k in range(len(b)):   # traversing through row of b
                ans[i][j]+=float(a[i][k]*b[k][j])
    for i in range(len(ans)):
        print(ans[i])
    
'''******************'''    
class myComplex:
    def __init__(self,real,imag):
        self.real=real
        self.imag=imag
        
    def add(self,acom):
        return myComplex(self.real+acom.real, self.imag+acom.imag)
    
    def prod(self,acom): #acom is another complex no.
        return myComplex((self.real*acom.real-self.imag*acom.imag),(self.real*acom.imag+self.imag*acom.real)) # adds real & imag part of acom to previously defined complex no.
    
    def mod(self):
        import math
        print(math.sqrt(self.real**2+self.imag**2))
    
    def show(self):
        print(self.real,f'{self.imag:+}','i')
        
'''******************''' 


