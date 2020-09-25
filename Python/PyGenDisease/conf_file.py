
class ConfFile:
    def __init__(self):
        self.key_vals = {}


    def InitFromFile(self, path):
        # Read user defined key-vals from config file.
        with open(path, 'r') as f:
            self.InitKeyVals([l for l in f])


    def InitKeyVals(self, lines):
        for l in lines:
            l = l.strip()
            if len(l) == 0 or l.startswith('#'):
                continue                    # Skip comments and empty lines

            key_val = [sp.strip() for sp in l.split(':')]
            self.key_vals[key_val[0]] = key_val[1]

        # Add hard-coded key-vals
        self.key_vals['NUM_CHROMOSOMES'] = '2'


    def Print(self, log=None):
        for k in self.key_vals:
            key_val_str = '{} -> {}'.format(k, self.key_vals[k])
            print(key_val_str) if log == None else log(key_val_str)



if __name__ == '__main__':
    print('This is a module :)')

