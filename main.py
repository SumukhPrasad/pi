import math
from time import time


def sqrtOfBigNumber(n, one):
     """
     Return the square root of n as a fixed point number with the one
     passed in.  It uses a second order Newton-Raphson convergence.  This
     doubles the number of significant figures on each iteration.
     """
     floating_point_precision = 10**16
     n_float = float((n * floating_point_precision) // one) / floating_point_precision
     x = (int(floating_point_precision * math.sqrt(n_float)) * one) // floating_point_precision
     n_one = n * one
     while 1:
         x_old = x
         x = (x + n_one // x) // 2
         if x == x_old:
             break
     return x


def piChudnovsky(one=1000000):
     """
     Calculate pi using Chudnovsky's series
     This calculates it in fixed point, using the value for one passed in
     """
     i = 1
     a_i = one
     a_sum = one
     b_sum = 0
     C = 640320
     C3_OVER_24 = C**3 // 24
     while 1:
          a_i *= -(6*i-5)*(2*i-1)*(6*i-1)
          a_i //= i*i*i*C3_OVER_24
          a_sum += a_i
          b_sum += i * a_i
          i += 1
          if a_i == 0:
               break
     total = 13591409*a_sum + 545140134*b_sum
     pi = (426880*sqrtOfBigNumber(10005*one, one)*one) // total
     return pi 

if __name__ == "__main__":
     for i in range(1,7):
          digits = 10**i
          one = 10**digits

          start =time()
          pi = piChudnovsky(one)
          f = open("pi.txt", "w")
          f.write(str(pi))
          f.close()
          print("chudnovsky: digits",digits,"time",time()-start)