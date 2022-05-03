# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 17:32:18 2022

@author: Osian ap Sion s1811930

A Solar System Simulation and Analysis Using Newtonian Physics with the Beeman Numerical Integration Algorithm for the First Four Planetary Bodies with an Artificial Satellite
"""
#importing necessary files and Agg backend for matplotlib export to mp4 files
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
matplotlib.use('Agg')

#class defining an object within the solar system
class Solar_system_object():
    #Initialising variables. Three acceleration variables are needed for Beeman algorithm
    def __init__(self, name, size, colour, mass, position, velocity):
        self.name = name
        self.size = size
        self.colour = colour
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.acceleration_minus = np.array([0,0])
        self.acceleration = self.acceleration_minus.copy()
        self.acceleration_plus = self.acceleration_minus.copy()
    #copy function in case a duplicate object is needed
    def copy(self):
        return Solar_system_object(self.name, self.size, self.colour, self.mass, self.position.copy(), self.velocity.copy())
        
#class defining the main simulation
class Solar_system_simulation():
    #function to help with processing strings from the initial values file 
    def clean_string(self, string, character):
        output_string = ''
        for char in string:
            if char == character:
                continue
            else: output_string += char
            
        return output_string
    #function for converting initial values from text file to objects and initialising the variables for the instance of this class
    def get_values(self):
        
        filename = 'data/Initial Values.txt'#str(input("File to open: "))
        filein = open(filename, "r")
        lines = filein.readlines()
        
        self.object_list = []
        
        objects = False
        simulation = False
        
        for line in lines:
            
            current_line = self.clean_string(line, ' ')
            
            if current_line.startswith('*'):
                
                if (current_line[current_line.count('*'):len(current_line) - 1] == 'OBJECTS'):
                    objects = True
                    
                elif current_line[current_line.count('*'):len(current_line) - 1] == 'SIMULATION':
                      simulation = True
                      objects = False
            #Process the part of the text file that defines the simulation objects          
            elif objects == True:
                
                parameter_list = current_line.split(',')
                new_parameter_list = []
                
                for parameter in parameter_list:
                    
                    if parameter.isalpha() == False:
                        if parameter.find('e') != -1:
                            parameter = parameter.split('e')
                            parameter = float(parameter[0])*(10**float(parameter[1]))
                        else: 
                            parameter = float(parameter)
                            
                    new_parameter_list.append(parameter)
                #append defined objects to list for later use             
                self.object_list.append(Solar_system_object(new_parameter_list[0], float(new_parameter_list[1]), new_parameter_list[2], float(new_parameter_list[3]), np.array([float(new_parameter_list[4]), float(new_parameter_list[5])]), np.array([float(new_parameter_list[6]), float(new_parameter_list[7])])))
            
            #Process the part of the text file that defines simulation variables
            elif simulation == True:
                
                parameter_list = current_line.split(',')
                new_parameter_list = []
                
                for parameter in parameter_list:
                    
                    if parameter.isalpha() == False:
                        if parameter.find('e') != -1:
                            parameter = parameter.split('e')
                            parameter = float(parameter[0])*(10**float(parameter[1]))
                        else: 
                            parameter = float(parameter)
                    
                    new_parameter_list.append(parameter)
                    
                self.size = int(new_parameter_list[0])
                self.density = int(new_parameter_list[1])
                self.gravitational_constant = float(new_parameter_list[2])
                self.precision = float(new_parameter_list[3])
                self.duration_input = float(new_parameter_list[4])
                
        filein.close()
    
    def __init__(self):
        
        self.get_values()
        self.duration = int(15768*self.duration_input*self.precision)
        self.time_step = int(2000/self.precision)
        self.time = 0.0
        
        
        #Initalise variables for storing and accessing data
        self.log = []
        self.object_variables = []
        self.object_index = 0
        self.size_index = 1
        self.colour_index = 2
        self.mass_index = 3
        self.position_index = 4
        self.velocity_index = 5
        self.acceleration_index = 6
    
    #calculate the acceleration due to other objects, acceleration is not conserved in this model
    #Calculate the change in acceleration, velocity and position for each time step
    def step_acceleration(self, obj):
        
        obj.acceleration_minus = obj.acceleration
        obj.acceleration = obj.acceleration_plus
        obj.acceleration_plus = 0
        for other in self.object_list:
            if obj != other:
                position_sum_of_squares = float(np.sum(np.power(other.position - obj.position, 2)))
                obj.acceleration_plus += (self.gravitational_constant*(np.array((obj.mass*other.mass)/position_sum_of_squares))*((other.position - obj.position)/np.sqrt(position_sum_of_squares)))/obj.mass
            
    #update position, velocity and acceleration for one time step 
    def step_velocity(self, obj):
        
        #update acceleration
        
        obj.velocity = obj.velocity + (1/6)*(2*obj.acceleration_plus + 5*obj.acceleration - obj.acceleration_minus)*self.time_step
    
    def step_position(self, obj):
        obj.position = obj.position + obj.velocity*self.time_step +(1/6)*(4*obj.acceleration - obj.acceleration_minus)*(self.time_step**2)
   
    
    def time_step_forward(self):
        #iterate through all objects and run their inividual variable update functions
        
        for obj in self.object_list:
            self.step_acceleration(obj)
            self.step_velocity(obj)
        #positions updated seperately to improve simulation accuracy
        for obj in self.object_list:
            self.step_position(obj)
            
        self.time += self.time_step
    #begin the simulation
    def run(self):
        
        for t in range(self.duration):
            self.log_variables()
            self.time_step_forward()
       
    #store object details for a single time-step    
    def log_variables(self):
        self.object_variables = []
        for obj in range(len(self.object_list)):
            self.object_variables.append([])
            self.object_variables[obj].append(self.object_list[obj].name)
            self.object_variables[obj].append(self.object_list[obj].size)
            self.object_variables[obj].append(str(self.object_list[obj].colour))
            self.object_variables[obj].append(self.object_list[obj].mass)
            self.object_variables[obj].append(self.object_list[obj].position.copy())
            self.object_variables[obj].append(self.object_list[obj].velocity.copy())
            self.object_variables[obj].append(self.object_list[obj].acceleration.copy())
        self.log.append(self.object_variables)
    
#Class that is used for the analysis of the simulation, makes extensive use of the log variable in the simulation class             
class Solar_system_analysis():
    
    def __init__(self, Solar_system_simulation):
        #initialise variables for the data analysis
        self.sim = Solar_system_simulation
        self.kinetic_energy_list = np.array([])
        self.potential_energy_list = np.array([])
        self.satellite_mars_distances_list = np.array([])
        self.satellite_earth_distances_list = np.array([])
        
        #Initalise variables for the animation
        self.current_time = 0
        self.xpos = np.linspace(-0.5*self.sim.size, 0.5*self.sim.size, self.sim.density)
        self.ypos = self.xpos
        self.sim_steps_per_frame = 100 #Not necessary to animate every simulation step, this makes the animation process considerably quicker
        self.niter = int((self.sim.duration/self.sim_steps_per_frame) - 1) #-1 to stop niter going out of bounds 
    #kinetic energy calculation for one object    
    def kinetic_energy(self, current_obj_trace):
        
        kinetic_energy = 0
        velocity_vec = current_obj_trace[self.sim.velocity_index]
        velocity_mag_sq = velocity_vec[0]**2+velocity_vec[1]**2
        kinetic_energy += (1/2)*(current_obj_trace[self.sim.mass_index])*(velocity_mag_sq)
        
        return kinetic_energy
        
    #potential energy calculation for one object
    def potential_energy(self, current_obj_trace, time):
        
        potential_energy = 0
        
        for other in self.sim.log[time]:
            
            if current_obj_trace[self.sim.object_index] != other[self.sim.object_index]:
                distance = np.sqrt(float(np.sum(np.power(other[self.sim.position_index] - current_obj_trace[self.sim.position_index], 2))))
                potential_energy += -(1/2)*self.sim.gravitational_constant*current_obj_trace[self.sim.mass_index]*other[self.sim.mass_index]/distance
        
        return potential_energy
    
    #calculate kinetic and potential energy for all simulation objects and print tese to file
    def calc_energies(self):
        fileout = open("Output\Total Energy.txt","w")
        fileout.write("The Total Energy of the System, in Joules: " + "\n\n")
        for t in range(self.sim.duration):
            total_kinetic_energy = 0
            total_potential_energy = 0
            
            for obj_trace in self.sim.log[t]:
                
                total_kinetic_energy += self.kinetic_energy(obj_trace)
                total_potential_energy += self.potential_energy(obj_trace, t)
                
            self.kinetic_energy_list = np.append(self.kinetic_energy_list, total_kinetic_energy)
            self.potential_energy_list = np.append(self.potential_energy_list, total_potential_energy)
            
            total_energy = total_kinetic_energy + total_potential_energy
            time = t*self.sim.time_step/(60*60*24*365.25)
            fileout.write(str(total_energy) +'J, ' + str(time) + ' Years ' + "\n")
        
        fileout.close()
    
    #calculate the orbital periods of each object and write to file
    #This works by timing the gap between the x coordinate changing from negative to positive
    def orbital_periods(self):
        
        fileout = open("Output\Orbital Periods.txt","w")
        fileout.write("The orbital periods of each object in Earth years: " + "\n\n")
        
        period_list = np.array([])
        for obj_trace in range(len(self.sim.log[0])):
            start_t = 0
            x_sign = 0
            x_sign_previous = 1
            delay = False
            begin_timer = False
            period_sum = 0
            orbits = 0
            
            for t in range(self.sim.duration):
                current_obj_trace = self.sim.log[t][obj_trace]
                position = current_obj_trace[self.sim.position_index]
                
                
                
                if x_sign != x_sign_previous and np.sign(x_sign) > 0:
                    if begin_timer == True:
                        period_sum += (t - start_t)*self.sim.time_step
                        start_t = t
                        orbits += 1
                    elif begin_timer == False:
                        if delay == True:
                            begin_timer = True
                        delay = True
                        
                        start_t = t
                      
                x_sign_previous = x_sign
                x_sign = np.sign(position[0])
            
            if orbits > 0:
                period = (period_sum/orbits)/(60*60*24*365.25)
            else:
                period = 'na'
            
            period_list = np.append(period_list, period)
                
        print(period_list)
        for period in range(len(period_list)):
            fileout.write(str(self.sim.log[0][period][self.sim.object_index])+ ": " + str(period_list[period]) + " Years\n")
        fileout.close()
    
    #calculate the distance of the satellite from mars and earth. writes the shortest mars distance to file
    def satellite_distances(self):
        
        satellite_index = 5
        mars_index = 4
        earth_index = 3
        self.satellite_mars_distances_list = []
        self.satellite_earth_distances_list = []
        
        for t in range(self.sim.duration):
            
            satellite_mars_distance = np.linalg.norm(self.sim.log[t][mars_index][self.sim.position_index] - self.sim.log[t][satellite_index][self.sim.position_index])
            satellite_earth_distance = np.linalg.norm(self.sim.log[t][earth_index][self.sim.position_index] - self.sim.log[t][satellite_index][self.sim.position_index])
            
            self.satellite_mars_distances_list.append(satellite_mars_distance)
            self.satellite_earth_distances_list.append(satellite_earth_distance)
        
        fileout = open("Output\Other Data.txt","w")
        fileout.write("The distance the satellite gets to mars: " + str(min(self.satellite_mars_distances_list)/1000) + " km\n")
        fileout.write("In a time of " + str(self.satellite_mars_distances_list.index(min(self.satellite_mars_distances_list))*self.sim.time_step/(60*60*24*365.25)) + " Years")
        
        fileout.close()
    #Graph all of the collected data from the previous functions and export to file
    #axes and figures are cleared before plotting to overwriite previous graphs exported to file 
    def create_graphs(self):
        self.x = np.linspace(0, self.sim.duration*self.sim.time_step/(60*60*24*365.25), self.sim.duration) #define x axis for the data, since they all run the length of the simulation
        s1 = self.potential_energy_list                                                                    #define y axes
        s2 = self.kinetic_energy_list
        s3 = self.potential_energy_list + self.kinetic_energy_list
        s4 = self.satellite_mars_distances_list
        s5 = self.satellite_earth_distances_list
        
        plt.figure(1)
        plt.clf()
        plt.cla()
        plt.ylabel("Potential Energy (J)")
        plt.xlabel("Time (years)")
        plt.title("Potential Energy Against Time")
        plt.plot(self.x, s1)
        plt.savefig('Output\Potential Energy Over Time.png', dpi = 300)
        
        plt.figure(2)
        plt.clf()
        plt.cla()
        plt.ylabel("Kinetic Energy (J)")
        plt.xlabel("Time (years)")
        plt.title("Kinetic Energy Against Time")
        plt.plot(self.x, s2)
        plt.savefig('Output\Kinetic Energy Over Time.png', dpi = 300)
        
        plt.figure(3)
        plt.clf()
        plt.cla()
        plt.ylabel("Total Energy (J)")
        plt.xlabel("Time (years)")
        plt.title("Total Energy Against Time")
        plt.plot(self.x, s3)
        plt.savefig('Output\Total Energy Over Time.png', dpi = 300)
        
        plt.figure(4)
        plt.clf()
        plt.cla()
        plt.ylabel("Distance (m)")
        plt.xlabel("Time (years)")
        plt.title("Distance From Satellite to Mars Against Time")
        plt.plot(self.x, s4)
        plt.savefig('Output\Distance From Satellite to Mars Over Time.png', dpi = 300)
        
        plt.figure(5)
        plt.clf()
        plt.cla()
        plt.ylabel("Distance (m)")
        plt.xlabel("Time (years)")
        plt.title("Distance From Satellite to Earth Against Time")
        plt.plot(self.x, s5)
        plt.savefig('Output\Distance From Satellite to Earth Over Time.png', dpi = 300)
        
  
    #animate the simulation
    def animate(self, i):
        
        for obj_trace in range(len(self.sim.log[self.current_time])):

            current_obj_trace = self.sim.log[self.current_time][obj_trace]
            self.patches[obj_trace].center = (current_obj_trace[self.sim.position_index][0], current_obj_trace[self.sim.position_index][1])
            
        self.current_time += self.sim_steps_per_frame
        
        return self.patches
    
    #run the animation
    def run_animation(self):
        fig = plt.figure(10)  #make sure the figure is not referring to an already defined figure
        plt.clf()
        plt.cla()
        ax = plt.axes()
        
        
        # Set up formatting for the animation file
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='Osian ap Sion'), bitrate=3600)
        
        #create a list of animation objects to append to
        self.patches = []
        
        #append the animation objects 
        for obj_trace in self.sim.log[self.current_time]:
            self.patches.append(plt.Circle((obj_trace[self.sim.position_index][0], obj_trace[self.sim.position_index][1]), radius = obj_trace[self.sim.size_index], color = str(obj_trace[self.sim.colour_index]), animated = True))

        for i in range(len(self.patches)):
            ax.add_patch(self.patches[i])
            
        ax.axis('off')
        ax.set_xlim(self.xpos[0], self.xpos[-1])
        ax.set_ylim(self.ypos[0], self.ypos[-1])
        plt.title("Solar System Animation")
        
        #animate and save the animation to file
        self.anim = animation.FuncAnimation(fig, self.animate, frames = self.niter, repeat = False, interval = 5, blit = True)
        self.anim.save('Output\Solar System Animation.mp4', writer=writer)

    def run(self):
        self.calc_energies()
        self.orbital_periods()
        self.satellite_distances()
        self.create_graphs()
        self.run_animation()
            
        