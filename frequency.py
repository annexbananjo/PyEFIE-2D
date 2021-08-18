import numpy as np
from math import pi
from object import Object


class Frequency(Object):
    """ Class frequency
    """

    def __init__(self, idObj=-1, freq=1e9):
        """ Init the frequency object
        """
        Object.__init__(self, idObj)
        if freq <= 0.0:
            raise ValueError("Frequency must be greater than zero")
        self.freq = freq
        self.omega = 2.0*pi*freq

    def __str__(self):
        """print info about the freq
        """
        s0 = 'Frequency obj with idObj: {}'.format(self.getId())
        s2 = '(freq, omega): {:.3f}; {:.3f}'.format(self.freq, self.omega)
        return s0 + "\n" + s2


if __name__ == '__main__':
    freq = Frequency()
    print(freq)
