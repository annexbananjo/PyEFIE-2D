import numpy as np
from object import Object
from point2D import Point2D, compute_distance


class Rectangle(Object):
    """ Create a rectangle object
    """

    def __init__(self, idObj=-1, pMin=Point2D(), pMax=Point2D()):

        self.pMin = pMin
        self.pMax = pMax
        self.type = "Rectangle"
        Object.__init__(self, idObj)

    def isInternal(self, aPoint):
        """Check if a point is within this rectangle
        """
        if (aPoint.x >= self.pMin.x and aPoint.x <= self.pMax.x) \
                and (aPoint.y >= self.pMin.y and aPoint.y <= self.pMax.y):
            return True
        else:
            return False


# Print the class Rectangle


    def __str__(self):
        s0 = 'Object is: '
        s1 = self.type
        s11 = 'idObj: {}'.format(self.getId())
        s2 = 'pMin: {}'.format(self.pMin)
        s3 = 'pMax: {}'.format(self.pMax)
        return s0 + s1 + "\n" + s11 + "\n" + s2 + "\n" + s3


class Circle(Object):
    """ Create a circle object
    """

    def __init__(self, idObj=-1, radius=0.0, center=Point2D()):
        self.radius = radius
        self.center = center
        self.type = "Circle"
        Object.__init__(self, idObj)

    def isInternal(self, aPoint):
        """Check if a point is within this circle
        """
        if compute_distance(self.center, aPoint) <= self.radius:
            return True
        else:
            return False

# Print the class Rectangle
    def __str__(self):
        s0 = 'Object is: '
        s1 = self.type
        s11 = 'idObj: {}'.format(self.getId())
        s2 = 'center: {}'.format(self.center)
        s3 = 'radius: {}'.format(self.radius)
        return s0 + s1 + "\n" + s11 + "\n" + s2 + "\n" + s3


if __name__ == '__main__':
    rec1 = Rectangle()

    pMin = Point2D(1, 2)
    pMax = Point2D(10, 20)
    rec2 = Rectangle(pMin=pMin, pMax=pMax)
    print(rec2.isInternal(pMax))

    circ1 = Circle()

    circ2 = Circle(radius=10.0, center=Point2D(0.0, 0.0))
    print(circ2.isInternal(pMax))
    print(circ2.isInternal(Point2D()))

    # append a couple of mix objs
    print("\n")
    vPrObj = []
    vPrObj.append(rec1)
    vPrObj.append(rec2)
    vPrObj.append(circ1)
    vPrObj.append(circ2)
    for obj in vPrObj:
        print(obj)
        print("\n")
