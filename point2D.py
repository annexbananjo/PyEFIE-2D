import numpy as np
from numpy.lib.scimath import sqrt


class Point2D():
    """ Class for handle a simple point in 2D
    """

    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y

    def __lt__(self, other):
        if (self.x < other.x) and (self.y < other.y):
            return True
        else:
            return False

    # Print the class point
    def __str__(self):
        s1 = '(x, y) = {}, {}'.format(self.x, self.y)
        return s1


def compute_distance(p1, p2):
    """ Compute the distance between two points
    """
    dx = p1.x - p2.x
    dy = p1.y - p2.y
    return sqrt(dx**2 + dy**2)


if __name__ == '__main__':
    p = Point2D()
    print(p)
    p1 = Point2D(1.3, 454.9)
    print(p1)
    # Compute distance beween two points
    print(compute_distance(p, p1))

    p2 = Point2D(0.0, 1.0)
    p3 = Point2D(-1.0, 0.0)
    print(compute_distance(p2, p3))

    p4 = Point2D(1.0, 0.0)
    p5 = Point2D(-1.0, 0.0)
    print(compute_distance(p4, p5))

    print(p4 < p5)
    print(p5 < p4)
    print(p4 < p1)
