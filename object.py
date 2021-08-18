# A base class object
class Object():
    """ A class Obj
    """

    def __init__(self, idObj=-1):
        self.idObj = idObj

    def getId(self):
        return self.idObj

    def isInternal(self, aPoint):
        """ Abstract method
        """
        raise NotImplementedError()
