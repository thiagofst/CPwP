class rungekutta(object):

    def __init__(self, *functions):

        """
        Runge-Kutta 4th order method to solve ODEs """

        self.f = functions
        self.t = 0


    def solve(self, y, h, n):

        """
        Solve an ODE system """

        t = []
        res = [] # Empty arreys for results
        for i in y:
            res.append([])

        while self.t <= n and h != 0: # While t is smaller than n and h is not equal to 0
            t.append(self.t)
            y = self._solve(y, self.t, h)
            for c, i in enumerate(y):
                res[c].append(i)

            self.t += h # t = t + n

            if self.t + h > n:
                h = n - self.t

        return t, res


    def _solve(self, y, t, h):

        functions = self.f

        k1 = []
        for f in functions:
            k1.append(h * f(t, *y))

        k2 = []
        for f in functions:
            k2.append(h * f(t + .5*h, *[y[i] + .5*h*k1[i] for i in range(0, len(y))]))

        k3 = []
        for f in functions:
            k3.append(h * f(t + .5*h, *[y[i] + .5*h*k2[i] for i in range(0, len(y))]))

        k4 = []
        for f in functions:
            k4.append(h * f(t + h, *[y[i] + h*k3[i] for i in range(0, len(y))]))

        return [y[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0 for i in range(0, len(y))]
