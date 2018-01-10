'''
Created on 20/02/2015
@author: Ronak H Shah

'''
import filecmp
import os
import sys
from subprocess import Popen
import shlex
import nose
import logging


def setup_module():
    this_dir, this_filename = os.path.split(__file__)
    new_dir = os.path.dirname(this_dir)
    inputFile = os.path.join(new_dir, "data", "test", "testData.txt")
    ctFile = os.path.join(new_dir, "data", "canonicalInfo",
                          "canonical_transcripts_cv6.txt")
    outFilePrefix = "testResult"
    outFileTxt = os.path.join(new_dir, "testResult_Annotated.txt")
    outFileXlsx = os.path.join(new_dir, "testResult_Annotated.xlsx")
    outFileJson = os.path.join(new_dir, "testResult_Annotated.json")
    outFileFtxt = os.path.join(new_dir, "testResult_functional.txt")
    cmpFile = os.path.join(new_dir, "data", "test", "testResult.txt")
    scriptFile = os.path.join(new_dir, "iAnnotateSV.py")
    cmd = "python " + scriptFile + " -r hg19 -i " + inputFile + " -ofp " + \
        outFilePrefix + " -o " + new_dir + " -c " + ctFile + " -v"
    logging.info("test_iAnnotateSV:%s", cmd)
    args = shlex.split(cmd)
    if(os.path.isfile(outFileTxt)):
        os.remove(outFileTxt)
        os.remove(outFileXlsx)
        os.remove(outFileJson)
        os.remove(outFileFtxt)
    try:
        proc = Popen(args)
        proc.wait()
        retcode = proc.returncode
        if(retcode >= 0):
            pass
    except:
        e = sys.exc_info()[0]
        logging.info(
            "Running of python command: %s \n has failed. The exception produced is %s Thus we will exit", cmd, e)
        sys.exit(1)


def teardown_module():
    this_dir, this_filename = os.path.split(__file__)
    new_dir = os.path.dirname(this_dir)
    outFileTxt = os.path.join(new_dir, "testResult_Annotated.txt")
    outFileXlsx = os.path.join(new_dir, "testResult_Annotated.xlsx")
    outFileJson = os.path.join(new_dir, "testResult_Annotated.json")
    outFileFtxt = os.path.join(new_dir, "testResult_functional.txt")
    if(os.path.isfile(outFileTxt)):
        os.remove(outFileTxt)
        os.remove(outFileXlsx)
        os.remove(outFileJson)
        os.remove(outFileFtxt)


def test_text_fileSimilarity():
    this_dir, this_filename = os.path.split(__file__)
    new_dir = os.path.dirname(this_dir)
    outFileTxt = os.path.join(new_dir, "testResult_Annotated.txt")
    cmpFileTxt = os.path.join(new_dir, "test", "testResult_Annotated.txt")
    nose.tools.ok_(filecmp.cmp(outFileTxt, cmpFileTxt),
                   msg="The current result text file and the original result text file for iAnnotate are not the same")


if __name__ == '__main__':
    nose.main()
