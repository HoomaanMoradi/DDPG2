from math import cos, sin, pi
import numpy as np
import matplotlib.pyplot as plt



class doublePendulum():

    def __init__ (self, g=9.81, l1=0.042, l2=0.042, m1=0.055, m2=0.055, ts=0.05):

        self.grav=g
        self.len1=l1
        self.len2=l2
        self.mass1=m1
        self.mass2=m2
        self.stepsize=ts
        self.theta1=pi
        self.omega1=0
        self.theta2=pi
        self.omega2=0

    def reset(self):
        self.theta1=pi
        self.omega1=0
        self.theta2=pi
        self.omega2=0
        st_res = np.zeros((4,), dtype=np.float32)
        st_res[0] = self.theta1
        st_res[1] = self.omega1
        st_res[2] = self.theta2
        st_res[3] = self.omega2
        return st_res

   



    def step(self, torque1, torque2):

    
        dif = self.theta1 - self.theta2

        numerator1   = -(self.grav*((2*self.mass1)+self.mass2)*sin(self.theta1)) - (self.mass2*self.grav*sin(self.theta1 - (2*self.theta2))) - (2*sin(dif)*self.mass2*(((self.omega2**2)*self.len2) + ((self.omega1**2)*self.len1*cos(dif))))+(torque1)
        denominator1 = self.len1*((2*self.mass1) + self.mass2 - (self.mass2*(cos((2*self.theta1) - (2*self.theta2)))))

        omega1prime = numerator1/denominator1

        numerator2   = 2*sin(dif)*(((self.omega1**2)*self.len1*(self.mass1+self.mass2)) + (self.grav*(self.mass1+self.mass2)*cos(self.theta1)) + ((self.omega2**2)*self.len2*self.mass2*cos(dif)))+(torque2)
        denominator2 = self.len2*((2*self.mass1) + self.mass2 - (self.mass2*cos((2*self.theta1) - (2*self.theta2))))

        omega2prime = numerator2/denominator2

        self.omega1 += self.stepsize * omega1prime
        self.omega2 += self.stepsize * omega2prime
        self.theta1 += self.stepsize * self.omega1
        self.theta2 += self.stepsize * self.omega2

        if self.theta1 > 2*pi:
            ct = int(self.theta1 / 2 / pi)
            self.theta1 -= ct * 2 * pi
        elif self.theta1 < -2 * pi:
            ct = int((-self.theta1)/2/pi)
            self.theta1 += ct * 2 * pi
        if self.omega1 > 10.0:
            self.omega1 = 10.0
        elif self.omega1 < -10.0:
            self.omega1 = -10.0


        if self.theta2 > 2*pi:
            ct = int(self.theta2 / 2 / pi)
            self.theta2 -= ct * 2 * pi
        elif self.theta2 < -2 * pi:
            ct = int((-self.theta2)/2/pi)
            self.theta2 += ct * 2 * pi
        if self.omega2 > 10.0:
            self.omega2 = 10.0
        elif self.omega2 < -10.0:
            self.omega2 = -10.0

        st_res = np.zeros((4,), dtype=np.float32)
        st_res[0] = self.theta1
        st_res[1] = self.omega1
        st_res[2] = self.theta2
        st_res[3] = self.omega2
        rwd = -5*self.theta1*self.theta1-5*self.theta2*self.theta2-1*self.omega1*self.omega1-1*self.omega2*self.omega2


        return st_res, float(rwd)



def main():
    env = doublePendulum()
    st = env.reset()
    ts_rec = []
    omega1points = []
    omega2points = []
    theta1points = []
    theta2points = []
    for i in range(100):
        ts_rec.append(i)
        theta1points.append(st[0])
        omega1points.append(st[1])
        theta2points.append(st[2])
        omega2points.append(st[3])
        st, _ = env.step(-0.013, -0.013)

     
if __name__ == '__main__':
    main()
    
