import numpy as np
import matplotlib.pyplot as plt
from IPython.display import clear_output
import os
import contextlib
np.seterr(all='raise')

class NoSignChange(Exception):
    def __init__(self, message='No sign change in specified bounds'):
        # Call the base class constructor with the parameters it needs
        super(NoSignChange, self).__init__(message)

class BoundError(Exception):
    def __init__(self, message='Bounds on v0 are divergent'):
        # Call the base class constructor with the parameters it needs
        super(BoundError, self).__init__(message)

class solver():

    def __init__(self, gamma=1, a=10):
        """
        Parameters:
            gamma: numerical constant, from equation of state
            a: numerical constant, equal to GM/cs^2 r0
        """
        self.gamma = gamma
        self.a = a
        self.W0 = np.log(1)
        
    def ndf1(self, psi, coords):
        '''Dimensionless f1 function, takes t, (x, y, z) and returns 1-exp[2y-z]'''
        try:
            #return 1-np.exp(2*coords[1])
            return 1-np.exp(2*coords[1]-coords[2])
        except OverflowError or UnderflowError:
            print("Error in ndf1 at ", psi, ", ", coords)
            return -1

    def ndf2(self, psi, coords):
        '''Dimensionless f2 function, takes t, (x, y, z) and returns (GM/r0*cs^2)exp[-x-z]-2'''
        try:
            #return np.exp(-coords[0])*a-2
            return self.a*np.exp(-coords[0]-coords[2])-2
        except OverflowError or UnderflowError:
            print("Error in ndf2 at ", psi, ", ", coords)
            return -1

    def dx(self, t, state):
        '''Infinitesimal x step'''
        return self.ndf1(t, state)
        #return state[0]
    def du(self, t, state):
        '''Infinitesimal u step'''
        return self.ndf2(t, state)
        #return state[1]
    def dw(self, t, state):
        '''Infinitesimal w step'''
        return -(self.gamma-1)*(2*self.ndf1(t, state)+self.ndf2(t, state))
        #return state[2]

    def adx(self, t, state):
        '''Absolute value of dx'''
        return abs(self.ndf1(t, state))
    def adu(self, t, state):
        '''Absolute value of du'''
        return abs(self.ndf2(t, state))
    def adw(self, t, state):
        '''Absolute value of dw'''
        return -(self.gamma-1)*(2*abs(self.ndf1(t, state))+abs(self.ndf2(t, state)))

    def coupledRKStep(self, funs, state, dt):
        """RK4 step for array of coupled functions that evolve a state, for step size dt
        Returns array (t, new state)
        funs and state[1] should have the same length

        Parameters:
            funs: array of functions [dx, du, dw] to evolve a state
            state: numpy array [t, [x, u, w]]
            dt: generalized time step size

        Returns:
            numpy array [t, [x, u, w]] incremented by one RK4 step
        """
        assert isinstance(funs, list), "Expected array input"
        for i in funs:
            assert callable(i), "Expected functions in list"
        assert isinstance(state, (np.ndarray, list)), "Expected array input"
        assert np.shape(state) == (2, ), "Expected 1D array with 2 elements"
        assert isinstance(state[1], (np.ndarray, list)), "Expected array in second element"
        assert np.shape(state[1]) == (len(funs), ), "State and function arrays are different lengths"
        assert isinstance(dt, (float, int, np.float64)), "Expected float input"


        karr = np.array([np.zeros(4) for i in range(len(funs))])
        for j in range(4):
            for i in range(len(funs)):
                if j == 0:
                    karr[i][j] = dt*funs[i](state[0], state[1])
                else:
                    karr[i][j] = dt*funs[i](state[0]+dt/2, state[1]+np.array([karr[k][j-1]/2 for k in range(len(karr))]))
        for i in range(len(funs)):
            karr[i][1] = 2*karr[i][1]
            karr[i][2] = 2*karr[i][2]
        step = np.array([state[1][i]+np.sum(karr[i])/6 for i in range(len(karr))])
        return np.array([state[0]+dt, step])


    def percentChange(self, curr, step):
        """Takes a current state and a new state (vectors) and calculates the percent change

        Parameters:
            curr: numpy array (x, u)
            step: numpy array (x', u')

        Returns:
            Percent change from curr to step"""

        assert isinstance(curr, (np.ndarray, list)), "Expected array input"
        assert isinstance(step, (np.ndarray, list)), "Expected array input"
        assert len(curr) == len(step), "Array lengths are not the same"

        return 100*abs(np.linalg.norm((step-curr)/np.linalg.norm(curr)))

    def adaptRK(self, currState, pc, t, dt, funs):
        """Iteratively adapts dt to ensure the change in (x, u) is between .1 and 1 percent
        This new dt can be either larger or smaller, independent of the initial percent change
        Takes a state (x, u, w), a percent change, t, dt, and the set of functions (dx, du, dw) being used to evolve the system

        Parameters:
            currState: numpy array (x, u, w)
            pc: percent change between currState and the next RK step with t, dt
            t: generalized time of currState
            dt: generalized time step size
            funs: array of functions [dx, du, dw] to evolve a state

        Returns:
            dt, adjusted so that the percent change is between .1 and 1 percent"""

        assert isinstance(funs, (np.ndarray, list)), "Expected array input"
        for i in funs:
            assert callable(i), "Expected functions in list"
        assert isinstance(currState, (np.ndarray, list)), "Expected array input"
        assert np.shape(currState) == (len(funs), ), "Expected 1D array input"
        assert isinstance(dt, (float, int, np.float64)), "Expected float input"
        assert isinstance(t, (float, int, np.float64)), "Expected float input, got {0}".format(type(t))
        assert isinstance(pc, (float, int, np.float64)), "Expected float input"

        ddt = 1.5
        i = 0
        itermax = 10000
        if pc > 1:
            pc2 = 1e10 #initialize dummy percent changes, used to track movement in % change while finding new dt
            prevpc2 = 1e10
            while pc2 > 1 and i < itermax:
                dt = dt*ddt
                step2 = self.coupledRKStep(funs, [t, currState], dt) #calculate hypothetical next step using new dt
                pc2 = self.percentChange(currState, step2[1])
                if pc2 > prevpc2: #if we're moving in the wrong direction, invert our change in dt
                    ddt = 1/ddt
                prevpc2 = pc2
                i = i+1
            if i == itermax: print("Max iteration count exceeded in adaptation")
            return dt #once we've found a working dt, take a step using it
        elif pc < .1: 
            pc2 = 1e-10 #initialize dummy percent changes, used to track movement in % change while finding new dt
            prevpc2 = 1e-10
            while pc2 < .1 and i < itermax:
                dt = dt*ddt
                step2 = self.coupledRKStep(funs, [t, currState], dt) #calculate hypothetical next step using new dt
                pc2 = self.percentChange(currState, step2[1])
                if pc2 < prevpc2: #if we're moving in the wrong direction, invert our change in dt
                    ddt = 1/ddt
                prevpc2 = pc2
                i = i+1
            if i == itermax: print("Max iteration count exceeded in adaptation")
            return dt #once we've found a working dt, take a step using it

    def generateFunc(self, state0, itermax=10000, AV=True, xrange=10, urange=5):
        """Generates a trace of wind behavior with initial conditions x0 and u0 (dimensionless) using the RK4 method with adaptation in dt
        Takes x0, u0, max iteration count and returns a 2D array tracking t, x, and u

        Parameters:
            state0: initial state of system (x0, u0, w0)
            itermax: maximum iteration count for loop
            AV: use the absolute value of f1 and f2 (boolean)
            xrange: maximum x value to display on plot (displays (1, xrange))
            urange: maximum u value to display on plot (displays (0, urange))

        Returns:
            numpy array, [0] contains t values, [1] contains x values, [2] contains u values, [3] contains w values for the wind curve"""

        assert isinstance(state0, (list, np.ndarray)), "Expected array input"
        assert len(state0) == 3, "Expected state array of length 3"
        assert isinstance(itermax, int), "Expected integer input"
        assert itermax > 0, "Expected positive itermax"
        assert isinstance(AV, bool), "Expected boolean input"
        assert isinstance(xrange, (float, int, np.float64)), "Expected numerical input"
        assert xrange > 1, "Expected xrange > 1"
        assert isinstance(urange, (float, int, np.float64)), "Expected numerical input"
        assert urange > 0, "Expected positive urange"

        xsol = np.array([state0[0]])
        usol = np.array([state0[1]])
        wsol = np.array([state0[2]])
        tarray = np.array([0])
        t = 0
        dt = .01
        i = 0
        currState = np.array([xsol[-1], usol[-1], wsol[-1]])

        if AV:
            funs = [self.adx, self.adu, self.adw]
        else:
            funs = [self.dx, self.du, self.dw]

        #Main loop for adaptive RK solver
        #Exit conditions are based on values for exp(x)
        #Using zero points for f1 and f2 only works if you change ndf1 and ndf2 to return absolute values, and then you don't see the full curve
        #Setting a max iteration count doesn't always work well - the solution curves may be cut
        while np.exp(currState[0]) > 1e-6 and np.exp(currState[0]) < xrange and i < itermax:

            #Load the current position of the system to determine if adaptation is necessary
            currState = np.array([xsol[-1], usol[-1], wsol[-1]])

            #Calculate the next integration step using the RK function defined above
            step = self.coupledRKStep(funs, [t, currState], dt)

            #Calculate percent change from current state to predicted next state
            
            pc = self.percentChange(currState, step[1])

            #If the percent change is too large or too small, change dt to compensate
            if pc > 1 or pc < .1:
                dt = self.adaptRK(currState, pc, t, dt, funs)
                step = self.coupledRKStep(funs, [t, currState], dt)

            xsol = np.append(xsol, step[1][0]) #update solution curves with next step
            usol = np.append(usol, step[1][1])
            wsol = np.append(wsol, step[1][2])
            t = t+dt
            i = i+1
            tarray = np.append(tarray, t)
        return np.array((tarray, xsol, usol, wsol))

    def makePlot(self, u0, AV=True, xrange=10, urange=5):
        """Generates a plot of one wind curve
        Takes u0 and generates a curve
        No return, prints a plot of the wind curve

        Parameters:
            u0: initial value of u
            AV: use absolute value of f1 and f2 functions (boolean)
            xrange: maximum x value to display on plot (displays (1, xrange))
            urange: maximum u value to display on plot (displays (0, urange))

        Displays a plot of one wind curve"""

        assert isinstance(u0, (float, int, np.float64)), "Expected float input"
        assert u0 > 0, "Expected positive input"
        assert isinstance(AV, bool), "Expected boolean input"
        assert isinstance(xrange, (float, int, np.float64)), "Expected numerical input"
        assert xrange > 1, "Expected xrange > 1"
        assert isinstance(urange, (float, int, np.float64)), "Expected numerical input"
        assert urange > 0, "Expected positive urange"

        plt.figure(1)
        plt.xlim(1, xrange)
        plt.ylim(0, urange)

        func = self.generateFunc(np.array([0, np.log(u0), self.W0]), 10000, AV, xrange, urange)
        plt.scatter(np.exp(func[1]), np.exp(func[2]-func[3]/2), s=.5);
        plt.title("Velocity vs Radius (Dimensionless), v0 = %1.9f cs"%u0)
        plt.xlabel("r/r0")
        plt.ylabel("v/cs")
        

    def makePlots(self, vmin, vmax, dv, AV=True, xrange=10, urange=5, showAll=False):
        """Generates a plot of wind curves
        Takes vmin, vmax, dv and generates a curve for u0 (x0 = 0), then increments v0 by dv and generates a new curve, repeating until v0 = umax
        Expects vmin, vmax, and dv scaled by 1/cs
        No return, prints a plot of the different wind curves

        Parameters:
            umin: starting u value
            umax: maximum u value
            du: increment of u
            AV: Use absolute value of f1 and f2 functions (boolean)
            xrange: maximum x value to display on plot (displays (1, xrange))
            urange: maximum u value to display on plot (displays (0, urange))
            showAll: boolean, True will plot additional graphs of velocity and temperature

        Displays a plot of several different wind curves"""

        assert isinstance(vmin, (float, int, np.float64)), "Expected float input"
        assert vmin > 0, "Expected positive input"
        assert isinstance(vmax, (float, int, np.float64)), "Expected float input"
        assert vmax > 0, "Expected positive input"
        assert isinstance(dv, (float, int, np.float64)), "Expected float input"
        assert dv > 0, "Expected positive input"
        assert vmax > vmin, "Maximum less than minimum"
        assert isinstance(AV, bool), "Expected boolean input"
        assert isinstance(xrange, (float, int, np.float64)), "Expected numerical input"
        assert xrange > 1, "Expected xrange > 1"
        assert isinstance(urange, (float, int, np.float64)), "Expected numerical input"
        assert urange > 0, "Expected positive urange"

        plt.figure(1)
        plt.xlim(1, xrange)
        plt.ylim(0, urange)


        for i in np.arange(vmin, vmax, dv):
            func = self.generateFunc(np.array([0, np.log(i), self.W0]), 10000, AV, xrange, urange)
            if i == vmin:
                data = [func]
            else:
                data.append(func)

        plt.figure(1)
        plt.title("Velocity vs Radius (Dimensionless)")
        for i in range(len(data)):
            plt.scatter(np.exp(data[i][1]), np.exp(data[i][2]-data[i][3]/2), s=.5, label='v0/cs = %g'%np.exp(data[i][2][0]));
        plt.xlabel("r/r0")
        plt.ylabel("v/cs")
        plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
        if showAll:
            plt.figure(2)
            plt.title("Temperature vs Radius (Dimensionless)")
            for i in range(len(data)):
                plt.scatter(np.exp(data[i][1]), np.exp(data[i][3]), s = .5, label = 'v0/cs = %g'%np.exp(data[i][2][0]));
            plt.xlabel("r/r0")
            plt.ylabel("v/cs")
            plt.legend(loc = "upper left", bbox_to_anchor = (1, 1))
            plt.figure(3)
            plt.title("Temperature vs Velocity (Dimensionless)")
            for i in range(len(data)):
                plt.scatter(np.exp(data[i][2]), np.exp(data[i][3]), s = .5, label = 'v0/cs = %g'%np.exp(data[i][2][0]));
            plt.xlabel("r/r0")
            plt.ylabel("v/cs")
            plt.legend(loc = "upper left", bbox_to_anchor = (1, 1))
        

    def findZeros(self, v0):
        """Finds the t values where f1 and f2 reach zero, and returns the difference
        Takes starting velocity v0, expected to be scaled by 1/cs
        Returns tu-tx

        Parameters:
            v0: initial value of v/cs"""

        assert isinstance(v0, (float, int, np.float64)), "Expected numerical input"
        assert v0 > 0, "Expected positive input"

        u0 = np.log(v0)
        xsol = np.array([0])
        usol = np.array([u0])
        wsol = np.array([self.W0])
        t = 0
        dt = .01
        currState = np.array([xsol[-1], usol[-1], wsol[-1]])
        xfound = False
        ufound = False
        tu = 0
        tx = 0

        #Main loop, uses RK solver and iterates until f1 or f2 changes sign, and returns the t value where that takes place    
        while not ufound and not xfound:

            #Load the current position of the system to determine if adaptation is necessary
            currState = np.array([xsol[-1], usol[-1], wsol[-1]])

            #Calculate the next integration step using the RK function defined above
            step = self.coupledRKStep([self.adx, self.adu, self.adw], np.array([t, currState]), dt)

            #Calculate percent change from current state to predicted next state
            
            pc = self.percentChange(currState, step[1])

            #If the percent change is too large or too small, change dt to compensate
            if pc > 1 or pc < .1:
                dt = self.adaptRK(currState, pc, t, dt, [self.adx, self.adu, self.adw])
                step = self.coupledRKStep([self.adx, self.adu, self.adw], np.array([t, currState]), dt)

            #if ndf2 changes sign, we have found its turnover point and can exit the loop
            if np.sign(self.ndf2(t, step[1])) != np.sign(self.ndf2(t, currState)):
                xfound = True
                tx = t

            #if ndf1 changes sign, we have found its turnover point and can exit the loop
            if np.sign(self.ndf1(t, step[1])) != np.sign(self.ndf1(t, currState)):
                ufound = True
                tu = t

            #update solution curves with next step
            xsol = np.append(xsol, step[1][0])
            usol = np.append(usol, step[1][1])
            wsol = np.append(wsol, step[1][2])
            t = t+dt

        #return difference between zeros in ndf1 and ndf2
        return (tu-tx)

    def findVboundary(self, guess, increment=1e-4, maxprecision=1e-10, itermax=10000):
        """Locates the boundary value of v0 at which f1 and f2 pass through zero at the same time using a bisection method

        Parameters:
            guess: starting value of v0
            increment: initial step size for incrementing v0
            maxprecision: sets limit on how small dv can be before exiting the loop
            itermax: maximum iteration count

        Returns:
            Boundary value for v0"""

        assert isinstance(guess, (float, int, np.float64)), "Expected numerical input"
        assert guess > 0, "Expected positive input"
        assert isinstance(increment, (float, int, np.float64)), "Expected numerical input"
        assert isinstance(maxprecision, (float, int, np.float64)), "Expected numerical input"
        assert maxprecision > 0, "Expected positive input"
        assert isinstance(itermax, int), "Expected integer input"
        assert itermax > 0, "Expected positive input"

        v0 = guess
        dv = increment
        i = 0

        #Main loop, while dv is larger than the max precision we want to continue refining our search
        while abs(dv) > maxprecision and i < itermax:

            #When we find a sign change in the zero for a given v0, that implies that the solution has changed character between a breeze and a non-physical one
            #i.e. whichever of f1 and f2 changed sign first, that has now reversed
            #when that happens, we shrink dv until we avoid the sign change, in order to approximate where exactly that occurs in v0
            while (v0+dv) <= 0:
                dv = dv/2
                
            self.findZeros(v0+dv)
            
            while np.sign(self.findZeros(v0+dv)) != np.sign(self.findZeros(v0)):
                dv = dv/2
                while (v0+dv) <= 0:
                    dv = dv/2
            v0 = v0+dv
            i = i+1
        if i >= itermax:
            print("Max iteration count exceeded")
            raise Exception("MaxIterExceeded")
        return v0+dv
    
    def findV0(self, lowerguess, upperguess, dv, maxprecision=1e-10, itermax=10000, xrange=10, urange=5, show=True, showPlot=False):
        """Finds boundary values for v0 and estimates an exact value for it

        Parameters:
            lowerguess: estimate of lower bound on v0/cs
            upperguess: estimate of upper bound on v0/cs
            dv: initial step size for incrementing v0/cs
            maxprecision: sets limit on how small dv can be before exiting the loop
            itermax: maximum iteration count
            xrange: maximum x value to display on plot (displays (1, xrange))
            urange: maximum u value to display on plot (displays (0, urange))
            show: show plot of v0 once it is found (boolean)

        Returns:
            Average of upper and lower bounds on v0

        Prints bounds on v0, estimated value, and estimated error
        Displays plot of the wind curves for the boundary values of v0 if show = True"""

        assert isinstance(lowerguess, (float, int, np.float64)), "Expected numerical input"
        assert lowerguess > 0, "Expected positive input, got %"%lowerguess
        assert isinstance(upperguess, (float, int, np.float64)), "Expected numerical input"
        assert upperguess > 0, "Expected positive input, got %f"%upperguess
        assert isinstance(dv, (float, int, np.float64)), "Expected numerical input"
        assert dv > 0, "Expected positive input, got %"%dv
        assert isinstance(maxprecision, (float, int, np.float64)), "Expected numerical input"
        assert maxprecision > 0, "Expected positive input"
        assert isinstance(itermax, int), "Expected integer input"
        assert itermax > 0, "Expected positive input"
        assert isinstance(xrange, (float, int, np.float64)), "Expected numerical input"
        assert xrange > 1, "Expected xrange > 1"
        assert isinstance(urange, (float, int, np.float64)), "Expected numerical input"
        assert urange > 0, "Expected positive urange"
        
        if np.sign(self.findZeros(upperguess)) == np.sign(self.findZeros(lowerguess)):
            print("No sign change in specified range")
            raise NoSignChange
        dv=(upperguess-lowerguess)/2
        upper = self.findVboundary(upperguess, -dv, maxprecision, itermax)
        lower = self.findVboundary(lowerguess, dv, maxprecision, itermax)
        
        if show:
            print("Lower bound on v0: ", lower)
            print("Upper bound on v0: ", upper)
            if abs(lower-upper) > 1e-2:
                print("divergent bounds, exiting")
                raise BoundError
            print("Estimated v0: ", (lower+upper)/2)
            print("Estimated error: ", abs((upper-lower)/2))
        if showPlot:
            self.makePlot(lower, False, xrange, urange)
            self.makePlot(upper, False, xrange, urange)

        return (lower+upper)/2
    
    def gammaSearch(self, a=10, g0=None, dg=.0025, glim=5/3, lower=1e-10, upper=.9, itermax=100):
        """Searches through gamma values for a given a and returns a table of gamma values and associated critical velocities
        
        Parameters:
            a: GM/cs^2 r0; constant value
            g0: starting gamma value, defaults to the gamma the class was initialized with
            dg: gamma increment (positive or negative)
            glim: gamma value at which the search will stop
            lower: estimate of lower bound on v0/cs for first gamma value
            upper: estimate of upper bound on v0/cs for first gamma value
            itermax: maximum iterations for searching a given gamma value for v0
            
        Returns:
            List of ordered pairs [gamma, v0]
            Prints scatter plot of v0 vs gamma
        """
        assert isinstance(a, (float, int, np.float64)), "Expected numerical input"
        assert a > 0, "Expected positive input"
        assert g0 is None or isinstance(g0, (float, int, np.float64)), "Expected numerical input"
        assert g0 is None or g0 > 0, "Expected positive input"
        assert isinstance(dg, (float, int, np.float64)), "Expected numerical input"
        assert isinstance(glim, (float, int, np.float64)), "Expected numerical input"
        assert isinstance(lower, (float, int, np.float64)), "Expected numerical input"
        assert lower > 0, "Expected positive input"
        assert isinstance(upper, (float, int, np.float64)), "Expected numerical input"
        assert upper > lower, "Nonsensical bounds"
        assert upper < 1, "Upper bound should be less than 1"
        
        gdata = np.array([])

        if a != None:
            self.a = a
        if g0 != None:
            self.gamma = float(g0)
        else:
            g0 = float(self.gamma)

        i = 0
        while (i < itermax) and ((self.gamma <= glim and np.sign(dg) == 1) or (self.gamma >= glim and np.sign(dg) == -1)):
            if i == 0: print("Searching gamma  = ", self.gamma)
            if self.gamma == g0:
                try:
                    with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
                        gdata = [self.gamma, self.findV0(lower, upper, (upper-lower)/2, show = False)]
                    self.gamma = self.gamma+dg
                    i = 0
                    clear_output()
                except NoSignChange:
                    if i == 0:
                        print("No sign change, decrementing bounds")
                    upper = float(lower)
                    lower = lower/2
                    i = i+1
            else:
                try:
                    with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
                        gdata = np.vstack((gdata, [self.gamma, self.findV0(lower, upper, (upper-lower)/2, show=False)]))
                    self.gamma = self.gamma+dg
                    i = 0
                    clear_output()
                except NoSignChange:
                    if i == 0:
                        print("No sign change, decrementing bounds")
                    upper = float(lower)
                    lower = lower/2
                    i = i+1
        if i >= itermax: 
            print("Max iteration count exceeded at gamma  = ", self.gamma)
            print("No sign change above v0/cs  = ", lower)
            
        if len(gdata) > 0:    
            plt.figure(1)
            plt.title("Critical velocity vs. Gamma")
            plt.xlabel("Gamma")
            plt.ylabel("v0/cs")
            plt.scatter(gdata[:, 0], gdata[:, 1])
        else:
            print("No data collected")
        return gdata