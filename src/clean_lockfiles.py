import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--jobs-dir', dest='jobs_dir', default='/epitopedata/jobs',
                    help='parent directory of jobs directories')
parser.add_argument('--debug', action='store_true',
                    help='turn on debug printouts')
args = parser.parse_args()


for root, dirs, files in os.walk(args.jobs_dir, followlinks=True):
    if 'lockfile' in files and not 'completed' in files:
        lockfilePath = os.path.join(root, 'lockfile')
        os.remove(lockfilePath)
        if args.debug:
            print 'removed ' + str(lockfilePath)
