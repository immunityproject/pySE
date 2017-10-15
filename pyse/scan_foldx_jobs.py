import os
from argparse import ArgumentParser
from os import path, mkdir, listdir, remove, chmod, getcwd
from random import random
from stat import S_IXUSR, S_IWUSR, S_IRUSR, S_IRGRP, S_IXGRP, S_IROTH, S_IXOTH
from shutil import copyfile, move, rmtree
from time import sleep
from socket import gethostname
import subprocess
import sys
import re
import logging


logging.basicConfig(filename='scan_foldx_jobs.log', level=logging.DEBUG)


parser = ArgumentParser()
parser.add_argument('--job', dest='job', default=None,
                    help='only this one job directory will be processed, regardless of lockfile')
parser.add_argument('--jobs-dir', dest='jobs_dir', default='/epitopedata/jobs',
                    help='parent directory of jobs directories')
parser.add_argument('--working-dir', dest='working_dir', default='foldx_job',
                    help='where to copy job files before running FoldX')
parser.add_argument('--foldx-exe', dest='foldx_exe',
                    default='FoldX.linux32', help='the name of the FoldX executable')
parser.add_argument('--foldx-dir', dest='foldx_dir', default='/epitopedata/ImmunityProject/foldx/foldx3b4hack_1',
                    help='where the FoldX distribution files can be found')
parser.add_argument('--exit', dest='exit_when_done', action='store_true',
                    help='if set, then the process exits after all jobs are completed')
args = None


def copyConfigurationFrom(dirName):
    config_files = ['list.txt', 'individual_list.txt']
    logging.debug('copying config files from ' + dirName)
    for fileName in config_files:
        sourceFileName = path.join(dirName, fileName)
        destFileName = path.join(args.working_dir, fileName)
        copyfile(sourceFileName, destFileName)


def retry(func, retries=3):
    for i in range(retries):
        try:
            func()
            return
        except:
            sleep(10)
    raise sys.exc_info()[1]


def moveResultsTo(dirName):
    for fileName in listdir(args.working_dir):
        if is_output_file(fileName):
            sourceFileName = path.join(args.working_dir, fileName)
            destFileName = path.join(dirName, fileName)
            retry(lambda: move(sourceFileName, destFileName))
            logging.debug('Moved ' + sourceFileName + ' to ' + destFileName)


def is_output_file(fileName):
    return fileName.lower().endswith('.fxout') or re.match(r'.*_\d+(?:_\d+)?\.pdb', fileName)


def deleteResults(dirName):
    for fileName in listdir(dirName):
        if is_output_file(fileName):
            retry(lambda: remove(path.join(dirName, fileName)))


def runFoldx(dirName):
    exe = path.join(getcwd(), args.working_dir, args.foldx_exe)
    logging.debug('Running ' + exe + ' -runfile runfile.txt')
    proc = subprocess.Popen(
        [args.foldx_exe, '-runfile', 'runfile.txt'], cwd=args.working_dir, executable=exe)
    result = proc.wait()
    foldxOutputFilename = path.join(dirName, 'foldx-output.txt')
    try:
        with open(foldxOutputFilename, 'w') as errFile:
            errFile.write("stdout:\n")
            errFile.write(proc.stdout)
            errFile.write("\n\nstderr:\n")
            errFile.write(proc.stderr)
            mkdir(path.join(args.working_dir, 'error-' + str(int(random(1000)))))
    except Exception, e:
        logging.error('Could not write ' + foldxOutputFilename, e)
    if result != 0:
        logging.error('FoldX returned code ' +
                      proc.returncode + ' in ' + dirName)
        deleteResults(dirName)
    else:
        moveResultsTo(dirName)
        logging.debug('Successfully completed FoldX run in ' + dirName)


def processDir(dirName):
    copyConfigurationFrom(dirName)
    runFoldx(dirName)


foundWork = False


def scan_jobs_dir():
    global foundWork
    while True:
        foundWork = False
        for root, dirs, files in os.walk(args.jobs_dir, followlinks=True):
            check_job_dir(root, files)
        if not foundWork:
            logging.info('No available jobs')
            if args.exit_when_done:
                exit(0)
        sleep(600)  # sleep before looping up to the top job again


def check_job_dir(dir_name, filenames):
    global foundWork
    logging.debug('Checking ' + dir_name)
    if 'individual_list.txt' not in filenames:
        return
    lockfile_name = 'lockfile'
    if lockfile_name in filenames:
        logging.debug(dir_name + ' is already locked')
    else:
        lockfile_path = path.join(dir_name, lockfile_name)
        try:
            with open(lockfile_path, 'w') as lock_file:
                lock_file.write(gethostname() + '\n')
            logging.debug('Created lockfile ' + lockfile_path)
        except:
            logging.warning('Could not create lockfile ' + lockfile_path)
        else:
            try:
                logging.debug('Processing ' + dir_name)
                processDir(dir_name)
                logging.debug('Completed ' + dir_name)
                completed_filename = path.join(dir_name, 'completed')
                try:
                    with open(completed_filename, 'w') as lock_file:
                        lock_file.write(gethostname() + '\n')
                    logging.debug('Created completed file ' +
                                  completed_filename)
                except Exception as e:
                    logging.warning(
                        'Could not write completed file ' + completed_filename)
            except Exception as e:
                logging.debug('Error in ' + dir_name)
                logging.exception(e)
                deleteResults(dir_name)
                remove(lockfile_path)
        foundWork = True


def setupWorkingDir():
    rmtree(args.working_dir)
    mkdir(args.working_dir)
    for fileName in listdir(args.foldx_dir):
        destFileName = path.join(args.working_dir, fileName)
        if not path.exists(destFileName):
            sourceFileName = path.join(args.foldx_dir, fileName)
            retry(lambda: copyfile(sourceFileName, destFileName))
    chmod(path.join(args.working_dir, args.foldx_exe), S_IXUSR |
          S_IWUSR | S_IRUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)


def main():
    global args
    args = parser.parse_args()
    setupWorkingDir()
    if args.job is not None:
        processDir(args.job)
    else:
        scan_jobs_dir()


if __name__ == '__main__':
    main()
