import sys
from contextlib import contextmanager
from functools import wraps
import os
import logging
from hashlib import sha256
from time import time, sleep
import numpy as np



if sys.version_info >= (3, 3):
    # Python 3.3 and up have a native `replace` method
    from os import replace
elif sys.platform.startswith("win"):
    def replace(src, dst):
        # TODO: on Windows, this will raise if the file is in use,
        # which is possible. We'll need to make this more robust over
        # time.
        try:
            os.remove(dst)
        except OSError:
            pass
        os.rename(src, dst)
else:
    # POSIX rename() is always atomic
    from os import rename as replace


@contextmanager
def atomic_write(filepath, binary=False, fsync=False):
    """ Writeable file object that atomically updates a file (using a temporary file). In some cases (namely Python < 3.3 on Windows), this could result in an existing file being temporarily unlinked.
    :param filepath: the file path to be opened
    :param binary: whether to open the file in a binary mode instead of textual
    :param fsync: whether to force write the file to disk
    """

    tmppath = filepath + '~'
    while os.path.isfile(tmppath):
        tmppath += '~'
    try:
        with open(tmppath, 'wb' if binary else 'w') as file:
            yield file
            if fsync:
                file.flush()
                os.fsync(file.fileno())
        replace(tmppath, filepath)
    finally:
        try:
            os.remove(tmppath)
        except (IOError, OSError):
            pass


def timeit(func, name=""):
    """Summary

    Args:
            func (TYPE): Description
    """
    @wraps(func)
    def wrapped(*args, **kwargs):
        start_time = time()
        func_name = name if name else func.__name__

        func_ret = func(*args, **kwargs)
        end_time = time()
        exec_time = end_time - start_time

        # if it's an instance method, logging.debug self.iterations
        if len(args) >= 1:
            logging.debug("function: {1}, {2} permutations, time elapsed: {0}".format(
                exec_time, func_name, args[0].iterations))
        else:
            logging.debug(
                "function: {1}, time elapsed: {0}".format(exec_time, func_name))

        logging.debug("----------")
        return func_ret

    return wrapped


@contextmanager
def open_file_or_string(filein, *args, **kwargs):
    try:
        with open(filein, *args, **kwargs) as f:
            yield f
    except FileNotFoundError:
        logging.warn(
            "Did not find file {filein}, resorting to treating {filein} as a string instead".format(filein=filein))
        string = filein
        yield [string]

def possible_genotypes(pos1, pos2):
    """
    takes in three-character string of form A/T
    returns list of possible combinations of first char and third char
                    (AA, AT, TT)
    """
    possible = [pos1 * 2, ''.join(sorted(pos1 + pos2)), pos2 * 2]
    return possible


def randomizer(first, second):
    """
    takes in tuple of 2 lists
    randomly shuffles elements across lists
    returns tuple of randomly shuffled lists
    """
    big = first[:]
    big.extend(second)
    np.random.shuffle(big)
    first_len = len(first)
    return (big[:first_len], big[first_len:])



@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    logging.debug(f"Currently in directory: {prevdir}")
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def make_working_dir(p_file, override_folder):
    # make working directory, if already exist, raise error
    working_dir = "{}_all_files".format(p_file)
    try:
        os.makedirs(working_dir)
    except OSError as e:
        # if e.errno == errno.EEXIST:
        if override_folder:
            remove_folder(working_dir)
            make_working_dir(p_file, override_folder)
        else:
            raise OSError(
                "Directory '{}' already exist, please choose another name or delete the existing one".format(working_dir))
        # else:
        #     raise OSError(
        #         "Ecountered problems different than directory exist. Please report error to the developers")

    return working_dir


def remove_folder(path):
    import shutil

    # Try to remove tree; if failed show an error using try...except on screen
    try:
        shutil.rmtree(path)
    except OSError as e:
        logging.debug("Error: %s - %s." % (e.filename, e.strerror))


def make_hash(self, *args):
    h = sha256()
    for arg in args:
        try:
            arg = arg.encode()
        except:
            pass
        h.update(arg)
    return h.hexdigest()  # hex()
