"""
Author: Swaroop Ramakant Avarsekar (2011172)
"""

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

def series(n):
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