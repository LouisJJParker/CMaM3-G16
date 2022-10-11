import math
from re import I
from tkinter import CENTER
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

class graph(object):
    def __init__(self):
        self.L = 1 # string length
        self.n_b = 17 # Number of oscillators
        self.rho = 400 # mass per unit length
        self.alpha = 0.25 #non-linear coefficient
        self.k_young = 0.1 # Young's modulus
        self.init_type = 1 # initial conditions

        self.dt = 0.1 # time step
        self.final_time = 12000 # total simulation time

    def init(self):
        return self.patches

    def create_oscillators(self):
        self.oscillator = []
        x_coords = np.linspace(0,self.L,self.n_b)
        for i in x_coords:
            self.oscillator.append(Oscillator(i, self.n_b))

    def animate(self, i):
        self.simulate()

        for j in range(len(self.oscillator)):
            self.patches[j].center = [self.oscillator[j].x,self.oscillator[j].y]
        return self.patches

    def simulate(self):
        self.find_acceleration()
        for i in self.oscillator:
            i.accelerate()
        for i in self.oscillator:
            i.move()

    def run(self):
        #Creating plot elements

        fig = plt.figure()
        ax = plt.axes()

        self.patches = []

        for i in range(self.n_b):
            self.patches.append(plt.Circle([self.oscillator[i].x,self.oscillator[i].y],0.01))

        for i in range(0, len(self.patches)):
            ax.add_patch(self.patches[i])

        #ax.axis('scaled')
        ax.set_xlim(0,1)
        ax.set_ylim(-1,1)

        anim = FuncAnimation(fig, self.animate, init_func = self.init, 
                            frames = self.final_time, repeat = False,
                            interval = 0.0001, blit = True)

        plt.show()

    def find_acceleration(self):
        for i in range(1,self.n_b-1):
            self.oscillator[i].a = ((self.k_young/self.rho) * 
            (self.oscillator[i-1].y+self.oscillator[i+1].y-2*self.oscillator[i].y) *
            (1+self.alpha*(self.oscillator[i+1].y-self.oscillator[i-1].y)))


class Oscillator(object):
    def __init__(self,x,i):
        self.i = i #index
        self.x = x #x position
        self.y = np.sin(2*np.pi*x) #displacement
        self.pos = np.array([self.x,self.y])
        self.a = 0 #acceleration
        self.v = 0 #velocity

    def accelerate(self):
        self.v += self.a*2

    def move(self):
        self.y += self.v*2
    

Graph=graph()
Graph.create_oscillators()
Graph.run()