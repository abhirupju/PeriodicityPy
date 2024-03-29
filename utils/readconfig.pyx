import ConfigParser

configs = {}
Config = None
def init(filename):
    global Config
    Config = ConfigParser.ConfigParser()
    Config.read(filename)
    """
    configs = {}
    for s in Config.sections():
        configs[s] = ConfigSectionMap(Config, s)
    """
    
def ConfigSectionMap(Config, section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

def getnParticle():
    Config.getint("Particles", "number")
