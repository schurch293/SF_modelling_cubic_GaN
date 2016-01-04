import pylab as lab
import math

"""
SCHRODINGER EQUATION SOLVER IN CONFINEMENT POTENTIAL - varying stacking fault separation

This code solves the Schrodinger equation using Numerov's method, defined by the equation relating adjacent wavefunction values:

    y_n = (2*y_n-1(1-(5h^2/12)*K_n-1)-y_n-2(1+(h^2/12)*K_n-2))/(1+(h^2/12)*K_n)

where h is separation between adjacent points in the wavefunction. The PDE is in the form:

    y'' + Ky =0

meaning that, for the shrodinger equation,

    psi''+(2*m/(hbar**2 * c**2))*(E-V)*psi = 0

    K_n = (2*m/(hbar**2 * c**2))*(E_n-V)
    
where parameters are in units of nm,eV throughout unless stated otherwise.

The dimensionless form of the equation is used, where:

    phi'' + b(e-V/E_0)phi = 0   

    phi = phi(z)

    z=x/x_0
    e=E/E_0

    B=2*m*x_0**2*E_0/(hbar**2 * c**2)

A solution is reached by calculating the wavefunction for an initial energy, using the Numerov equation. The energy is then modified
iteratively until the wavefunction at the edge of the potential converges to zero.

This code is specialised to stacking faults in cubic GaN, calculating the potential drop across the SF
for a range of SF separations, then using numerov's method to determine the ground state confinement energy
for each separation.
"""

#parameter definitions
DeltaSF = 0.7775*(10**(-9))   
Eoff = 0.07                         #band offset between Hex/Cubic GaN
#Epot = 0.303                #included calculation for this parameter
hbar = 6.582119514E-16
me = 510998.9
m=2*me                              #effective mass of carrier 
#Steps = 1000
MaxSwitch = 20
c=299792458                         #speed of light in free space

#calculation of potential drop across SF
Phex = -0.034   #added value, polarisation of Hex Gan, C/m^2
epcub = 9.7     #added value, dielectric constant of cubic GaN
ephex = 10.28   #added value, dielectric constant of Hex GaN
eps = 8.8541878*(10**(-12)) #added value, vacuum permittivity, F/m


#loop to vary separation of SFs

Max_sep = 20*DeltaSF         #Max separation of SFs considered
No_sep = 1             #number of separations to consider

L=lab.linspace(340*DeltaSF/3,340*DeltaSF/3,No_sep)      #separation of SFs to solve equation
ground_states = [0]*No_sep                  #array to store ground state energies
SF_drops = [0]*No_sep
s=0                                         #SF separation loop tracking parameter

PSIout = []

while s<No_sep:

    DeltaX=L[s]                 #SF separation

    MaxZ = (2*DeltaX+DeltaSF)/DeltaSF   #dimensionless potential range
    Epot = -Phex/(eps*((ephex/DeltaSF)+(epcub/(2*DeltaX))))      #output in eV
    stepsize = 0.01E-9
    Steps = (2*DeltaX+DeltaSF)/stepsize
    StepSize =MaxZ/Steps

    SF_drops[s] = Epot
    
    #normalisation parameters
    x_0 = DeltaSF
    E_0 = Epot

    B=(2*m*(x_0**2)*E_0)/((hbar**2)*(c**2))


    #potential definition for SF system
    def V(z):
        x=x_0*z
        if x <DeltaX:                               #outside SF
            return ((-Epot)/(DeltaX))*x+(Epot+Eoff)
        elif x<DeltaX+DeltaSF:                      #inside SF
            return (Epot/DeltaSF)*(x-DeltaX)                    #changed 0.33 to Epot
        else:                                       #outside SF
            return ((-Epot)/(DeltaX))*(x-DeltaX-DeltaSF)+(Epot+Eoff)
    

    #wavefunction prefactor definition
    def K(X,E):
        return B*(E-(V(X)/E_0))



    xRange = lab.linspace(0,MaxZ,Steps)


    #normalised potential calculation
    potential = []                                  
    for Val in xRange:
        potential.append(V(Val)/E_0)


    #initialisation of loop parameters
    converge = False
    Energy = 0.1                   #initial energy
    delta = 0.01                    #initial energy change after iteration
    direction = 1                   #parameter to track phase of wavefunction at edge of potential after iteration
    switches = 0                    #number of wavefunction phase changes
    n_tries = 0
    epsilons = []                   #array for energy values after each iteration
    switchStep = []
    while converge == False:
        #if n_tries % 10 == 0:
        #    print( n_tries)
        y=0.0
        step = 2
        x=step*StepSize
        k_minus_2 = K(x-2*StepSize,Energy) # k_0
        k_minus_1 = K(x-StepSize, Energy) # k_1
        a = 0.1                         #value for 2nd initial wavefunction
        y_minus_2 = 0 # y_0
        y_minus_1 = a # y_1
        x_out = []
        y_out = []

        #loop to calculate dimensionless wavefunction across the potential via numerov
        while step<len(xRange)-2:
            #print(x,step,StepSize)
            step +=1
            x+= StepSize

            #calculation of numerov's equation, y is the dimensionless wavefunction
            k = K(x,Energy)
            b=(StepSize**2)/12
            y = ( 2*(1-5*b*k_minus_1) * y_minus_1 - (1+b*k_minus_2) * y_minus_2 ) / (1 + b * k)

            x_out.append(x)             #output of wavefunction
            y_out.append(y)

            y_minus_2 = y_minus_1       #preparation of numerov parameters for next loop
            y_minus_1 = y
            k_minus_2 = k_minus_1
            k_minus_1 = k

        #conditions for phase change of wavefunction at edge of potential
        if y_out[-1] >0:
            #no change in phase
            if direction >0:
              Energy = Energy +delta
          
            #change in phase
            else:
              switches +=1
              #direction *=-1
              print(str(switches))              #print to track progress of code
              delta = delta/2                   #reduce energy step for convergence in next iteration
              Energy = Energy+delta             #modify energy value
              epsilons.append(Energy)
              switchStep.append(n_tries)

              #plot of wavefunction for every 10 switches
              #if switches % 10==0:
               # maxy=max(y_out)
                #yplot = [y/maxy for y in y_out]
                #lab.plot(x_out, yplot, label="$\epsilon = "+repr(epsilon)+"$")
                #lab.legend(loc=1)
                #switches = 0
        elif y_out[-1] < 0:
            #no change in phase
            if direction <0:
                Energy = Energy +delta
            #change in phase
            else:
              switches +=1
              #direction *=-1
              print(str(switches))              #print to track progress of code
              delta = delta/2                   #reduce energy step for convergence in next iteration
              Energy = Energy-delta             #modify energy value
              epsilons.append(Energy)
              switchStep.append(n_tries)
              #if switches % 10==0:
              #  maxy=-min(y_out)
              #  yplot = [y/maxy for y in y_out]
              #  lab.plot(x_out, yplot, label="$E/E_0 = "+repr(Energy)+"$")
              #  lab.legend(loc=1)
                #switches = 0
            
        #convergence if wavefunction is zero at end of potential
        else:
            #print(epsilon,delta,y_out[0],y_out[-1],n_tries)
            converge = True
            ground_states[s] = Energy
        #force end of loop if reach max switches
        if (switches >=MaxSwitch):
            print(Energy,delta,y_out[0],y_out[-1],n_tries)
            converge = True
            ground_states[s] = Energy*E_0
        else:
            #print(Energy,delta,y_out[0],y_out[-1],n_tries)
            n_tries+=1
    s += 1
    #plot output figures
    PSIout.append(y_out)

'''
    #plot of dimensionless wavefunction solution
    lab.figure(2+5*s)
    norm = max(y_out)
    lab.plot(x_out, y_out, label="E/E_0 = "+repr(Energy)+"$") #legend label is normalised energy
    lab.xlabel("x")
    lab.ylabel("y")
    lab.title("Output")
    lab.legend(loc=1)
'''
    #plot of dimensionless potential, wavefunction and energy
norm = max(y_out)
lab.figure(3+5*s)
lab.plot(xRange,potential,color = 'b', label="Potential")
normalised = [(y/(3*norm))+Energy for y in y_out]       #normalised wavefunction to plot at energy level
lab.plot(x_out,normalised,color='r',label="E = "+repr(Energy*E_0)+" eV")      #legend contains energy in eV, changed units to ev from $
lab.axhline(Energy,color='r')
lab.xlabel("z(x/x_0)")      #changed x to z in label
lab.ylabel("Energy (E/E_0)")
lab.legend(loc=1)
'''
    #plot of energy vs iterations
    lab.figure(4+5*s)
    lab.plot(switchStep,epsilons)
    lab.xlabel("Iterations")
    lab.ylabel("Energy/E_0")

    #plot of magnitude of percentage deviation with iteration
    lab.figure(5+5*s)
    variation = [((((ep-epsilons[-1])/epsilons[-1])**2)**(0.5)) for ep in epsilons]
    lab.semilogy(switchStep,variation)
    lab.xlabel("Iterations")
    lab.ylabel("% deviation from result")
    lab.show()
'''

#output graph of separation variation
lab.figure(100)
lab.scatter(L*10**9,ground_states)
lab.xlabel("SF separation (nm)")
lab.ylabel("Ground State Energy (eV)")

##for n in range(0,30):
##    lab.figure(n)
##    lab.plot(x_out,PSIout[n])
##    lab.title(L[n])
    
lab.show()

#output separation variation to file
file = open('output.txt.','w')
file.write('x/m ,\t\t\t Potential drop/eV ,\t\t Ground State/eV ,\n')
for i in range(0,(No_sep)):
    print >>file,L[i],',\t',SF_drops[i],',  \t\t',ground_states[i], ','
file.close()
    

