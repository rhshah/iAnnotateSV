'''
Created on 20/02/2015
@author: Ronak H Shah

'''
import pytest
import filecmp
import os
from subprocess import Popen
import shlex
def TestIO():
    this_dir, this_filename = os.path.split(__file__)
    new_dir = os.path.dirname(this_dir)
    inputFile = os.path.join(new_dir, "data", "test", "testData.txt")
    outFile = os.path.join(this_dir, "testResult.txt")
    cmpFile = os.path.join(new_dir, "data", "test", "testResult.txt")
    scriptFile = os.path.join(new_dir, "iAnnotateSV.py")
    cmd = "python " + scriptFile + " -r hg19 -i " + inputFile + " -o " + outFile + " -d 3000"
    args = shlex.split(cmd)
    if(os.path.isfile(outFile)):
        os.remove(outFile)
    else:
        proc = Popen(args)
        proc.wait()
        retcode = proc.returncode
        if(retcode >= 0):
            code = 1
        else:
            assert 0
        assert filecmp.cmp(outFile, cmpFile)
TestIO()