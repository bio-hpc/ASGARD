#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import subprocess
import sys


class Execute(object):

    def __init__(self, cfg):
        self.cfg = cfg

    def run (self, command):
        """

            If command fails then exit()
        """
        try:
            if self.cfg.p_debug:
                print( BColors.WARNING + "Debug: "+command + BColors.ENDC)
                #out = subprocess.check_output(command+" 2>>"+ self.cfg.log, shell=True in command)

                out = subprocess.check_output(command, shell=' ' in command)

            else:
                #out = subprocess.check_output(command+" 2>>"+ self.cfg.log, shell=True in command)
                out = subprocess.check_output(command, shell=True)


        except subprocess.CalledProcessError as e:
            ret = e.returncode
            if ret in (1, 2):
                print ("failed")
                print(command)
                exit()

            elif ret in (3, 4, 5):
                print ("the command failed very much")
                print(command)
                exit()


        if sys.version_info[0] == 3:
            return str(out, 'utf-8')
        else:
            return str(out)

    def run_with_fail (self, command):
        try:
            if self.cfg.p_debug:
                print( BColors.WARNING + "Debug: "+command + BColors.ENDC)
                #out = subprocess.check_output(command+" 2>>"+ self.cfg.log, shell=True in command)
                out = subprocess.check_output(command, shell=' ' in command)
                print ("out: ",out)
            else:
                #out = subprocess.check_output(command+" 2>>"+ self.cfg.log, shell=True in command)
                out = subprocess.check_output(command, shell=' ' in command)
        except subprocess.CalledProcessError as e:
            ret = e.returncode
            if ret in (1, 2):
                return("failed")
            elif ret in (3, 4, 5):
                return ("the command failed very much")


        if sys.version_info[0] == 3:
            return str(out, 'utf-8')
        else:
            return str(out)
        

class BColors(object):
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def __init__(self):
        pass
