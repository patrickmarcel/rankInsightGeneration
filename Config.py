import configparser
import json

class Config:
    """
    A simple Singleton class.
    """

    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Config, cls).__new__(cls)
        return cls._instance

    def __init__(self, fpath, USER):  # Constructor with a parameter
        if not hasattr(self, 'initialized'):
            config = configparser.ConfigParser()
            self.fpath = fpath
            config.read(self.fpath)

            # Database connection parameters
            self.dbname = config[USER]['dbname']
            self.user = config[USER]['user']
            self.password = config[USER]['password']
            self.host = config[USER]['host']
            self.port = int(config[USER]['port'])

            # Cube info
            self.table = config["Common"]['table']
            self.measures = json.loads(config.get("Common", "measures"))
            self.groupbyAtt = json.loads(config.get("Common", "groupbyAtt"))
            self.sel = config["Common"]['sel']
            self.meas = config["Common"]['meas']
            self.measBase = config["Common"]['measBase']
            self.function = config["Common"]['function']
            self.prefs = json.loads(config.get("Common", "preferred"))

            self.initialized = True