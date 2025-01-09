class Lattice:
    """
    A simple Singleton class.
    """

    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Lattice, cls).__new__(cls)
        return cls._instance

    def __init__(self, sample):
        # Internal dictionary to store the views, naming convention : att1, att2, ... attn where attn is the comparison attribute
        self.shelf = dict()
