#import important library
import math
import matplotlib.pyplot as plt


phi = 3
E_F = 0
F = 1E7
x = []
y = []
for i in range(200):

    a = (1.54143E-6*F*F)/phi
    b = (6.83089E7 * math.pow(phi,1.5)/F)
    c = (1.02463E8 * math.sqrt(phi)*E_F)/F
    y0 = (3.79469E-4 * math.sqrt(F)/phi)
    #print("y0={}\n".format(y0))

    p = ((((3*1.179)*(math.log(8/ math.sqrt(1.179))/math.log(2.71828)))/8)+(1.179/16))*(y0*y0)
    q = ((3*1.179)*y0*y0*(math.log(y0)/math.log(2.71828))/8)
    r = ((1.179*1.179)*(1-(math.log(8/math.sqrt(1.179))/math.log(2.71828)))*pow(y0,4))/32
    s = ((1.179*1.179)*pow(y0,4)*(math.log(y0)/math.log(2.71828))/32)

    v_y0 = 1-p+q+r+s;

    dv_dy = (-2)*y0*((((3*1.179)*(math.log(8/math.sqrt(1.179))/math.log(2.71828)))/8)+(1.179/16))+((6*1.179)*y0*(math.log(y0)/math.log(2.71828))/8)+((3*y0*1.179)/8)+(((1.179*1.179)*(1-(math.log(8/math.sqrt(1.179))/math.log(2.71828)))*pow(y0,3))/8)+(1.179*1.179*pow(y0,3)*(math.log(y0)/math.log(2.71828))/8)+(1.179*1.179*pow(y0,3)/32)
    t_y0 = (v_y0-((2*y0*dv_dy)/3))

    tt_y0 = (t_y0*t_y0)

    u_y0 = (1+(c*t_y0))* math.exp(-(c*t_y0))

    #print("a={}\nb={}\nc={}\nv(y0)={}\ndv_dy={}\n".format(a,b,c,v_y0,dv_dy))

    J = ((a*(1-u_y0)*math.exp(-(b*v_y0)))/tt_y0)

    #print("t(y0)= {}\nt^2(yo)= {} \nu(y0)={}\n\ncurrent density: {} A/cm^2\n".format(t_y0,tt_y0,u_y0,J))

    #store x and y value here
    x.append(E_F)
    y.append(J)
    
    #increment
    E_F = E_F + 0.2

print(x)
print(y)

plt.figure()
plt.plot(x,y)
plt.title('The graph of Current Density:\n')
plt.xlabel('E_F')
plt.ylabel('J')
# plt.show()
plt.savefig("mygraph.png")