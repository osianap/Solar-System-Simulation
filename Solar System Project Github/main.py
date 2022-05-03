# -*- coding: utf-8 -*-
"""
Created on Tue May  3 20:28:01 2022

@author: 2010o
"""
import simulation.proj_main as simulation


def main():
    
    sim = simulation.Solar_system_simulation()
    sim.run()
    analysis = simulation.Solar_system_analysis(sim)
    analysis.run()
    
    
main()
