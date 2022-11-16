# add to initial conditions file
err = 5 # the percentage error for tolerence
# import it into main file like others

c = [-1,1] # creating a list number between -1 and 1

import random
def random_number_generator():
  f = random.choice(c)
  fx = f*random.random()
  return fx

tol = np.zeros(N)
for n in range(N):
  tol[n] = (100+err*(random_number_genrator()))/100
 
tol[0] = 0

# change the calculation of acceleration into the one below. Having tolerence
a_i[n] = (k/m) * (x_i[n-1]*tol[n-1] + x_i[n+1]*tol[n+1] - 2*x_i[n]*tol[n]) * (1 + alp*(x_i[n+1]*tol[n+1] - x_i[n-1]*tol[n-1]))

#implimenting this calculations into the main file will produce a result with tolerence
