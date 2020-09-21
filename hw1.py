#Gauss-Legendre quadrature subroutine
def quadrature_sum(a0,b0,f):
    x = []
    w = []
    S = []
    for i in range(0,20):
        x.append(0.5*(b0-a0)*M[i][0]+0.5*(b0+a0))
        w.append(0.5*(b0-a0)*M[i][1])
        S.append(w[i]*f(x[i]))
    return sum(S)
#Adaptive integrator subroutine
def adaptive_integrator(a0,b0,f,epsilon):
    #Start adpative integrator.
    left_p = []; right_p = []; s = []
    #Initialisation
    left_p.append(a0); right_p.append(b0); s.append(quadrature_sum(a0,b0,f))
    j = 0; J = 0    #uppercase J used to record number of subtintervals used (see line 38).
    I = 0
    for i in range(0,10000):
        c = 0.5*(left_p[j]+right_p[j])
        s1 = quadrature_sum(left_p[j],c,f)
        s2 = quadrature_sum(c,right_p[j],f)
        if abs(s1+s2-s[j])>epsilon:
            #Increase range of updated lists
            left_p.append(0);   right_p.append(0);  s.append(0)
            #Now, proceed with algo.
            left_p[j+1] = left_p[j]
            right_p[j+1] = 0.5*(left_p[j]+right_p[j])
            s[j+1] = s1
            left_p[j] = 0.5*(left_p[j]+right_p[j])
            s[j] = s2
            j += 1; J += 1
        else:
            I = I + s1 + s2
            j -= 1
            if (j==-1):
                return I,i,2*(J+1) #the integral has been evaluated


#Nodes and weights for M = 20 (from Mathematica)
M = [[-0.9931285991851, 0.0176140071391521], [-0.9639719272779, \
  0.0406014298003869], [-0.9122344282513, \
  0.0626720483341091], [-0.8391169718222, \
  0.083276741576705], [-0.7463319064602, \
  0.101930119817240], [-0.6360536807265, \
  0.118194531961518], [-0.5108670019508, \
  0.131688638449177], [-0.3737060887154, \
  0.142096109318382], [-0.2277858511416, \
  0.149172986472604], [-0.076526521133, \
  0.152753387130726], [0.076526521133, \
  0.152753387130726], [0.227785851142, \
  0.149172986472604], [0.373706088715, \
  0.142096109318382], [0.510867001951, \
  0.131688638449177], [0.636053680727, \
  0.118194531961518], [0.746331906460, \
  0.101930119817240], [0.839116971822, \
  0.083276741576705], [0.912234428251, \
  0.0626720483341091], [0.963971927278, \
  0.0406014298003869], [0.993128599185, 0.0176140071391521]]

import math as m
#special constants
pi = m.pi 
e = m.exp(1)

#Main program
inputtest = 'y'
while (inputtest =='y'):
    try:
        a0 = eval(input("Enter left endpoint."))
        b0 = eval(input("Enter right endpoint."))
        inputfun = input("Enter function of x.  Please preface any math function with m.  For example, for cos(x) use m.cos(x).")
        epsilon = float(input("Enter desired precision for integral."))
        f = lambda x: eval(inputfun)
        I,i,J = adaptive_integrator(a0,b0,f,epsilon)
        break
    except:
        print('Mistake in entering input')
        inputtest = input('Try again (y/n)?')

if (inputtest != 'y'):
    exit()
else:
    print('The integral is:')
    print(I)
    print('computed with')
    print(i)
    print('iterations')
    print('and')
    print(J)
    print('subintervals')
