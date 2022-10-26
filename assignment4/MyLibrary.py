"""
Author: Swaroop Ramakant Avarsekar (2011172)
"""

def polycoeff(x,y,k): #x-data, y-data, k-degree of polynomial
    def sumik(x,k):
        add=0
        for i in range(len(x)):
            add+=(x[i])**k
        return add

    def sumxyik(x,y,k):
        add=0
        for i in range(len(x)):
            add+=((x[i])**k)*y[i]
        return add
    X,Y=[],[]
    for i in range(k+1):  # Filling X,Y matrix with 0s accordingly
        rowx=[]
        rowy=[0]
        X.append(rowx)
        Y.append(rowy)
        for j in range(k+1):
            rowx.append(0)
    for i in range(len(X)):  # Filling X matrix with required elements
        for j in range(len(X)):
            if i==j:
                X[i][j]=sumik(x,2*i)

            if i!=j and i>j:
                X[i][j]=sumik(x,i+j)
            X[j][i]=X[i][j]
    for i in range(len(Y)): # Filling Y matrix with required elements
        Y[i][0]=sumxyik(x,y,i)
    return X,Y


def leastsq(x,y,w): #x-xdata, ydata, w-weights
    chisq=0
    sy,s,sx,sxx,sxy,syy=0,0,0,0,0,0
    for i in range(len(x)):
        sy+=y[i]/((w[i])**2)
        s+=1/((w[i])**2)
        sx+=x[i]/((w[i])**2)
        sxy+=x[i]*y[i]/((w[i])**2)
        sxx+=x[i]*x[i]/((w[i])**2)
        syy+=y[i]*y[i]/((w[i])**2)
    D=s*sxx-(sx)**2
    a1=(sxx*sy-sx*sxy)/D
    a2=(sxy*s-sx*sy)/D
    erra1,erra2=sxx/D,s/D
    r2=((sxy)**2)/(sxx*syy)
    return a1,a2,r2#,erra1,erra2


def deflation(cf,div): #(x-div)
    for i in range(1,len(cf)):
        cf[i]+=cf[i-1]*div
    cf.pop() # removing last element
    return cf

def der(cflist):
    n=len(cflist)-1
    arr = []
    for i in range(n+1):
        arr.append((n-i)*cflist[i])
    arr.pop() #popping last element for der
    return arr

def doubder(cflist): 
    n=len(cflist)-1
    arr = []
    for i in range(n+1):
        arr.append((n-i)*(n-i-1)*cflist[i])
    arr.pop()  #popping last 2 elements for doub der
    arr.pop()
    return arr

def in_func(cflist,x): # substituting x in function
    add=0
    for i in range(len(cflist)):
        add+=cflist[i]*(pow(x,(len(cflist)-i-1)))
    return add

def laguerre(cflist,b):
    import math
    if len(cflist) == 1: #breaks recursion
        return

    if abs(in_func(cflist,b))<0.0001: # if b is the solution
        print(b)
        cf=deflation(cflist,b)
        return laguerre(cf, b)
    else:
        n = len(cflist)-1
        G=in_func(der(cflist),b)/in_func(cflist,b)
        H=G**2-(in_func(doubder(cflist),b)/in_func(cflist,b))
        ap=G+math.sqrt(abs((n-1)*(n*H-G**2))) #with plus sign
        am=G-math.sqrt(abs((n-1)*(n*H-G**2))) #with minus sign
        if abs(ap)<abs(am):
            a=n/am
        else:
            a=n/ap
        b-=a
        return laguerre(cflist, b)

def bracketing(f,a,b): #a<b
    #global count
    count=0
    if f(a)*f(b)<0:
        #print('a')
        return a,b,count
    
    if f(a)*f(b)>0: #else  #f(a) and f(b) have same sign
        
        if abs(f(a))<abs(f(b)):
            a=a-0.5*(b-a)   # shifting a to left
            count+=1
            return bracketing(f,a,b)
        if abs(f(a))>abs(f(b)): 
            b=b+0.5*(b-a)  # shifting b to right
            count+=1
            return bracketing(f,a,b)
        
    if count>11:
        print('choose different intervals')

def bisection(f,a,b,E,D,ctr):
    a,b,count=bracketing(f,a,b)
    print(round(a,7),round(b,7),ctr)
    if abs(b-a)<E and f(a)<D: #precision check
        return round(a,7),round(b,7),ctr
    else:
        c=(a+b)/2 # bisecting the interval
        if f(c)*f(a)<0:
            ctr+=1
            return bisection(f,a,c,E,D,ctr) # shifting b
        if f(c)*f(b)<0:
            ctr+=1
            return bisection(f,c,b,E,D,ctr) #shifting a

def newraph(f,fd,x,E,D): #x-guess,E-epsilon,D-delta
    ctr=0
    a,b=x-f(x)/fd(x),x
    while abs(a-b)>E and abs(f(b))>D: #precision check
        b=a
        a=b-f(b)/fd(b)
        ctr+=1
        print(a,b)
    return a,b,ctr

def regula(f,a,b,E,D):
    c,ctr=b-((b-a)*f(b))/(f(b)-f(a)),0
    #print(round(c,7),ctr)
    if f(a)*f(b)<0:
        if abs(b-a)>E or (f(a)>D and f(b)>D):  #precision check
            while abs(a-b)>E:
                if f(a)*f(c)<0:
                    b=c
                if f(b)*f(c)<0:
                    a=c
                if abs(f(c))<D:
                    return round(b,7),round(c,7),ctr        #upto 6 decimal places      
                c=b-((b-a)*f(b))/(f(b)-f(a))
                ctr+=1    
                print(round(b,7),round(c,7),ctr)
    else:   # if the given bracket doesn't satisfy
        print(round(b,7),round(c,7),ctr) 
        a,b,count=bracketing(f,a,b)
        return regula(f,a,b,E,D)

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

def jacobi(new,a,b):
    for i in range(len(a)):
        temp=b[i][0]
        for j in range(len(a[0])):
            if i!=j:
                p=a[i][j]*new[j][0]
                temp=temp-p
        new[i][0]=temp/a[i][i]
    return new

def gauss_jacobi(c,d):
    x=[]
    for i in range(len(c)):
        x.append([0])
    check1=ddom_already(c,d)
    if check1==1: #already diagonally dominant 
        a,b=c,d
    else:   #make diagonally dominant
        a,b=ddom(c,d)
    prev,new=[[0]]*len(a),jacobi(x,a,b)
    count=0
    d=[]
    for i in range(len(a)):
        d.append(0)
    while True:
        new=jacobi(new,a,b)
        count+=1    
        for i in range(len(a)):
            d[i]=abs(new[i][0]-prev[i][0])
        if min(d)<0.000001:
            break
        else:
            count+=1
            for i in range(len(new)):
                prev[i][0]=new[i][0]
    return count,new


def ddom(a,b):
    sumlist=[]
    for row in range (0,len(a)):
        add=0
        for col in range (0, len(a[row])):
            add+=abs(a[row][col])
        sumlist.append(add)
    check=0
    for i in range (0,len(a)):
        if abs(a[i][i])>= (sumlist[i]/2): 
            pass
        else:
            for j in range (i+1,len(a)):
                if abs(a[j][i]) >= (sumlist[j]/2): 
                    a[j],a[i]=a[i],a[j] 
                    b[j],b[i]=b[i],b[j]
                    check = 1
    if check==1:  # if check is 1 implies diagonally dominant
        return a,b
    else:
        return ("diagonal dominant not possible")
    

def ddom_already(a,b):
    sumlist=[]
    for row in range (0,len(a)):
        add=0
        for col in range (0, len(a[row])):
            add+=abs(a[row][col])
        sumlist.append(add)
    check=0
    for i in range (0,len(a)):
        if abs(a[i][i])>=(sumlist[i]-a[i][i]): #already diagonally dominant
            check=1
        else:
            check=0
    return check
        
def gauss_seidel(a,b):
    '''check1=ddom_already(c,d)
    if check1==1: #already diagonally dominant 
        a,b=c,d
    else:   #make diagonally dominant
        a,b=ddom(c,d)'''
    guess=[] # Let this be the guess solution
    prev=[] # let prev list be filled with zeros
    for row in range(len(a)):
        guess.append([0])
        prev.append([0])
        
    count=0
    new=guess
    for somenum in range(100): # setting up iteration limit          
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
    print("iterations=",count)
    print('solution=',new)

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