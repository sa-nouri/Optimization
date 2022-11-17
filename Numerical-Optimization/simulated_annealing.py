#@title Import libraries

import numpy as np
import matplotlib.pyplot as plt

from random import randint, random 
from math import exp, log

#@title Simulated Annealing

class SimAnl():
    ''' Simulated Annealing'''

    def __init__(self, func, x0, temparature,  step_max=2000,
                 cooling_coe=0.9, bounds=[]):
        
        self.t = temparature
        self.main_temp = temparature
        self.step_max = step_max
        self.hist = []
        self.cost_func = func
        self.bounds = bounds[:]
        self.x0 = x0
        self.current_state = self.x0
        self.current_energy = func(self.x0)
        self.best_state = self.current_state
        self.best_energy = self.current_energy
        self.get_neighbor = self.move_continuous
        self.update_t = self.cooling_func
        self.cooling = cooling_coe

        # begin optimizing
        self.step, self.accept = 1, 0

        while self.step < self.step_max:
            proposed_neighbor = self.get_neighbor()
            E_n = self.cost_func(proposed_neighbor)

            dE = E_n - self.current_energy

            # determine if we should accept the current neighbor
            if random() < self.safe_exp(float(-dE) / float(self.t)):
                self.current_energy = E_n
                self.current_state = proposed_neighbor[:]
                self.accept += 1

            # check if the current neighbor is best solution so far
            if E_n < self.best_energy:
                self.best_energy = E_n
                self.best_state = proposed_neighbor[:]

            # persist some info for later
            self.hist.append([
                self.step,
                self.t,
                self.current_energy,
                self.best_energy])

            # update some stuff
            self.t = self.update_t(self.step)
            self.step += 1

        # generate some final stats
        self.acceptance_rate = self.accept / self.step


    def move_continuous(self):
        
        neighbor = [item + ((random() - 0.2)) for item in self.current_state]

        # clip to upper and lower bounds
        if self.bounds:
            for i in range(len(neighbor)):
                x_min, x_max = self.bounds[i]
                neighbor[i] = min(max(neighbor[i], x_min), x_max)

        return neighbor
    
    def cooling_func(self, step):
        return self.main_temp /  (1 + self.cooling * step)

    def safe_exp(self, x):
        try: return exp(x)
        except: return 0

    def results(self):
        print('+------------------------ RESULTS -------------------------+\n')

        print(f'  initial temp: {self.main_temp}')
        print(f'    final temp: {self.t:0.6f}')
        print(f'     max steps: {self.step_max}')
        print(f'    final step: {self.step}\n')
        print(f'  final energy: {self.best_energy:0.6f}\n')

        print('+-------------------------- END ---------------------------+')
        return self.hist


def func(x):
    return (x[0]**3) - 60*(x[0]**2) + 900*(x[0]) + 100

#@title Process

x0 = [19]
T0 = 100 #@param{}
opt = SimAnl(func, x0, step_max=20000, temparature=T0, cooling_coe=0.9, bounds=[])

results = np.array(opt.results())

steps = results[:, 0]
temperatures = results[:, 1]
cur_energies = results[:, 2]
best_energies = results[:, 3]

# Plot figures

plt.figure()
# plt.plot(temperatures, 'b')
plt.plot(cur_energies, 'b')
plt.xlabel('steps')
# plt.ylabel('temparatures')
plt.ylabel('Energies')
plt.title('Initial temperature: 100')
plt.grid()