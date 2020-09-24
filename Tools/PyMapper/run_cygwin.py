import os
import subprocess
import sys
import time

import logger



CYGWIN_PATH = 'c:\\cygwin64\\bin'



def AddCygwinPath():
    env_path = os.environ['PATH']
    paths = [ p.lower() for p in env_path.split(';') ]
    if CYGWIN_PATH not in paths:
        logger.Log('Add cygwin path . . .')
        new_env = CYGWIN_PATH + ';' + env_path
        os.environ['PATH'] = new_env
    else:
        logger.Log('cygwin path exists!')



def RunCommand(cmd):
    logger.Log(' START : ' + time.ctime(time.time()))
    cmd = [str(c) for c in cmd]
    logger.Log(' '.join(cmd))
    ret = subprocess.call(cmd)
    logger.Log('ret : ' + str(ret))
    logger.Log(' END   : ' + time.ctime(time.time()))
    return ret



def RemoveCygwinPath():
    env_path = os.environ['PATH']
    paths = [ p.lower() for p in env_path.split(';') ]
    if CYGWIN_PATH in paths:
        logger.Log('remove cygwin path . . .')
        paths.remove(CYGWIN_PATH)
        new_env = ''
        for p in paths:
            new_env += p + ';'
        os.environ['PATH'] = new_env



def main():
    AddCygwinPath()                 # Add cygwin path.
    RunCommand(sys.argv[1:])        # Run executable.
    RemoveCygwinPath()              # Clean up cygwin path.


if __name__ == '__main__':
    main()

