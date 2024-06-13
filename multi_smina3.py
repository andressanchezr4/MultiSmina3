# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 19:29:38 2024

@author: andres.sanchez
"""

import sys
import os
import subprocess as sp
import tempfile
import multiprocessing as mp
import itertools
import time
import gzip
import argparse

try:
    import progressbar
    PROGRESSBAR = True
except ImportError:
    print("'progressbar' module not found! Not using a progressbar")
    PROGRESSBAR = False
    

# int controlling amount of output written out, higher is more.
# default = 2
# 0 -> None
DEBUG = 2

# boolean for whether to clean up the tmp files
CLEANUP = True

# list of temp files to clean up
temp_files = []

# global variable for counting how many molecules are in the input file
TOTAL_COMPOUNDS = 0

# global variable for counting how many molecules have been processed so far
TOTAL_DONE = 0

def not_dollars(line):
    ''' 
    Check to see if the input line is the end of the molecule.
    
    Arguments:
    line -- line from sdf file
    
    Returns:
    boolean
    '''
    return "$$$$" != line.strip('\n')

def parse_molecules(file):
    '''
    A generator function that reads in sdf files 
    one molecule at a time.
    Can also read in a compressed .sdf.gz file
    
    Arguments:
    file -- sdf(.gz) file name
    
    Returns:
    block -- list of strings of lines of the molecule
    '''
    if DEBUG > 2:
        print('file:', file)
    if os.path.splitext(file)[1] == '.sdf': 
        with open(file, 'r') as lines:
            while True:
                block = list(itertools.takewhile(not_dollars, lines))
                if not block:
                    break
                yield block + ['$$$$\n']
                
    elif os.path.splitext(file)[1] == '.gz': 
        with gzip.open(file, 'rt') as lines:
            while True:
                block = list(itertools.takewhile(not_dollars, lines))
                if not block:
                    break
                yield block + ['$$$$\n']

def splitligs(ligfile, num):
    '''
    Split the input sdf file into the specified
    number of temp files. They are spread evenly
    among the partitions with the first ligand in
    the input going into the first partition and
    the second ligand going into the second partition
    and so on.
    
    Arguments:
    ligfile -- ligand file (.sdf)
    num -- number of file to split into
    
    Returns:
    tmpfiles -- list of temporary files
    '''
    global TOTAL_COMPOUNDS
    
    # Initialize the temporary files
    tmpfiles = []
    for i in range(num):
        # Create a temp file 
        tmpfile = tempfile.NamedTemporaryFile(mode='w+b', suffix='.sdf.gz', dir='/tmp/', delete=False)
        gzipfile = gzip.GzipFile(tmpfile.name, mode='wb', fileobj=tmpfile)
        tmpfiles.append(gzipfile)
        
    # Add the temp files to the global list of temporary files
    temp_files.extend([tmp.name for tmp in tmpfiles])
    
    # Index for which temp file to write to
    index = 0

    # For each compound in the sdf file, write to each tmp file
    for lig in parse_molecules(ligfile):
        TOTAL_COMPOUNDS += 1
        tmpfiles[index].write(''.join(lig).encode('utf-8'))  # Write to the temp file
        index += 1
        
        # Reset the counter
        if index >= num:
            index = 0
    
    # Throw a trailing endline and close the files
    for tmpfile in tmpfiles:
        tmpfile.write('\n'.encode('utf-8'))
        tmpfile.close()
    
    return tmpfiles

def reassemble_ligs(outfile, minfiles):
    '''
    De-interleave the results from smina.
    Takes the first compound from the first file and so on.
    
    Arguments:
    outfile -- output file name
    minfiles -- smina processed names
    
    Returns:
    outfile -- output file name
    '''
    # Open the minimized files and get the parse_molecules generator for each
    generators = [parse_molecules(minfile) for minfile in minfiles]
    
    if os.path.splitext(outfile)[1] == '.gz':
        fout = gzip.open(outfile, 'wt')  # Open output file in text mode for writing
    else:
        fout = open(outfile, 'w', encoding='utf-8')  # Open output file in UTF-8 encoding
    
    index = 0
    
    try:
        while True:
            for genny in generators:
                if DEBUG > 2:
                    print('reassemble_ligs:', index)
                # Get the lig from the next generator
                lig = next(genny)
                fout.write(''.join(lig))
                
    # When the first generator runs out of items, stop
    except StopIteration:
        if DEBUG > 3:
            print('Split files:', ' '.join([file for file in minfiles]))
        if DEBUG > 0: 
            print('Results reassembled:', outfile)
    
    fout.close()
    
    return outfile

def count_processed(stdout):
    '''
    Count the number of molecules that are processed in the
    bit of stdout pulled from the output of smina/smina
    #, and update the TOTAL_DONE variable.
    and return the value
    
    Arguments:
    stdout -- stdout from smina/smina output
    
    Returns:
    num_processed -- number processed in the input
    '''
    if DEBUG > 4:
        print('\ncount_processed:', stdout)
    num_processed = stdout.count('Affinity:')
    return num_processed
    
def wrap_run_smina_init(q):
    ''' Set the queue object into the function '''
    wrap_run_smina.q = q
    
def wrap_run_smina(argsque):
    '''
    Wrapper for run_smina that splits the arguments.
    
    Arguments:
    argsque -- tuple of smina arguments and the multiprocessing queue instance
    
    Returns:
    out -- output of run_smina(args, que)
    '''
    # Split the arguments
    args, que = argsque
    
    # Set the queue as an object of run_smina
    run_smina.q = wrap_run_smina.q
    
    # Pass to run_smina
    out = run_smina(args, que)
    
    return out

def run_smina(args, que):
    '''
    Run smina with the specified list of arguments.
    
    Arguments:
    args -- list of smina arguments
    que -- multiprocessing queue
    
    Returns:
    out -- dump of smina output
    '''
    if DEBUG > 3:
        print('run_smina: args:', args, 'queue:', que)
    
    # Try statement is to allow for Control-C to kill the child processes
    try:
        if DEBUG > 2:
            print('run_smina:', ' '.join(args))
            
        p = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
        
        # Dump of smina stdout
        out = ''
        while True:
            # Read in a line from the standard out
            output = p.stdout.readline()
            if DEBUG > 4:
                print('output:', output)
            
            # Add it to what has been seen so far
            out += output
            
            if DEBUG > 3:
                print('checking to break')
            
            # If there's nothing, you are done
            if output == '':
                if DEBUG > 3:
                    print('BREAKING')
                break
            
            # Count how many were processed
            num_processed = count_processed(output)
            if DEBUG > 3:
                print('number processed:', num_processed)
            # If some were processed, pass it back to the parent process
            if num_processed > 0:
                if DEBUG > 4:
                    print('PUTTING', num_processed, 'into queue')
                run_smina.q.put(num_processed)
        
        if DEBUG > 4:
            print('Leaving run_smina while loop')
             
    except KeyboardInterrupt:
        pass
    
    return out

def mult_run_smina(tmpfiles, outfile, args):
    '''
    Run each split file on smina using
    the specified arguments.
    
    Arguments:
    tmpfiles -- list of split ligand temporary files
    outfile -- output file name
    args -- list of smina arguments excluding the --ligand specification
    
    Returns:
    outfile -- output file name
    '''
    global TOTAL_DONE
    
    if DEBUG > 0:
        print('Running', sys.argv[4])
    if DEBUG > 3:
        print('splitfiles')
        for splitfile in tmpfiles:
            print(splitfile)
    
    maxproc = len(tmpfiles)
    
    # Initialize another set of temp files to capture the output of smina
    minfiles = []
    for i in range(maxproc):
        minfiles.append(tempfile.NamedTemporaryFile(mode='w+b', suffix='.sdf.gz', dir='/tmp/', delete=False))
    
    # Add the temp files to the global list of temp files
    temp_files.extend([tmp.name for tmp in minfiles])
    
    if PROGRESSBAR:
        # Set up the progress bar
        pbar = progressbar.ProgressBar(widgets=[progressbar.Percentage(), progressbar.Bar()], maxval=TOTAL_COMPOUNDS).start()
        
    # Make a queue object for communication between the parent and the children
    que = mp.Manager().Queue()
    
    # Initialize the multiprocessing pool
    pool = mp.Pool(None, wrap_run_smina_init, [que])
    
    # Add the tmp file name to the smina command and the queue
    smina_args = [(args + ['--ligand', tmpfile.name, '--out', minfile.name], que) for tmpfile, minfile in zip(tmpfiles, minfiles)]
    
    if DEBUG > 3:
        print('smina_args:')
        for line in smina_args:
            print(line)
    
    # Get the results. Wrapped in try-except to allow Control-C to kill the child processes
    try:
        # Setup the pool map
        p_map = pool.map_async(wrap_run_smina, smina_args)
        
        # While the processes are not done
        while not p_map.ready():
            if DEBUG > 3:
                print('mult_run_smina: TOTAL_DONE:', TOTAL_DONE)
            if not que.empty():
                if DEBUG > 3:
                    print('getting from queue')
                TOTAL_DONE += que.get(True)
            else:
                if DEBUG > 3:
                    print('sleeping')
                time.sleep(0.25)  # Check every quarter second
                
            if PROGRESSBAR:
                # Update the progress bar
                pbar.update(TOTAL_DONE)
            else:
                # Just print how many are done
                sys.stdout.write('%s\r' % str(TOTAL_DONE))
                sys.stdout.flush()
    
        if PROGRESSBAR:
            # Close the progress bar
            pbar.finish()
        else:
            print('')  # Write a blank line
        
        # Get the results and block
        results = p_map.get(True)
        
        # Close and join the pool 
        pool.close()
        pool.join()
        
    except KeyboardInterrupt:
        print('parent received control-c')
        return

    if DEBUG > 3:
        for result in results:
            for line in result.split('\n'):
                print(line)
    
    # Close the files
    [file.close() for file in tmpfiles]
    [file.close() for file in minfiles]
    
    # Count how many compounds are in the input and output files
    if DEBUG > 2:
        print('input: ', '\n'.join([file.name + ' ' + str(open(file.name, 'rb').read().count(b'$$$$')) for file in tmpfiles]))      
        print('output:', '\n'.join([file.name + ' ' + str(open(file.name, 'rb').read().count(b'$$$$')) for file in minfiles]))
    
    # Reassemble the ligands from the tmp files
    outfile = reassemble_ligs(outfile, [file.name for file in minfiles])

    # Delete the tmp files
    if CLEANUP:
        for minfile in minfiles:
            if os.path.exists(minfile.name):
                os.remove(minfile.name)
        for tmpfile in tmpfiles:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)

    return outfile

def clean_up(file_list):
    '''
    Delete all the temp files in the input list
    
    Arguments:
    file_list -- list of files to delete
    
    Returns:
    None
    '''
    for filen in file_list:
        if os.path.exists(filen):
            os.remove(filen)

def parse_args():
    '''
    Parse command line arguments from sys.argv
    
    Arguments:
    None -- reads from sys.argv
    
    Returns:
    cpus -- number of cpus to use 
    ligfile -- input ligand file name
    outfile -- output file name
    smina_commands -- commands for smina
    
    '''
    parser = argparse.ArgumentParser(description='Run smina (or vina) using multiple threads.\nPass arguments as you would to smina.\nIf you specify the smina executable, it must be before all other args except -o or -l.')
    parser.add_argument('-l', '--ligand', dest='ligfile', required=True, help='input ligand file (.sdf or .sdf.gz)')
    parser.add_argument('-o', '--out', dest='outfile', required=True, help='output ligand file (.sdf or .sdf.gz)')
    parser.add_argument('-c', '--cpu', dest='cpus', required=False, type=int, default=mp.cpu_count(), help='Number of CPUs to use. Default: detect number in system')
    
    args, unknown = parser.parse_known_args()
    
    ligfile = args.ligfile
    outfile = args.outfile
    cpus = int(args.cpus)
    smina_commands = ['/home/andres/Downloads/smina.static'] + unknown  # Assume smina is the command if not specified
    
    if DEBUG > 3:
        print('command line args:')
        print('sys.argv:', sys.argv)
        print('ligfile:', ligfile)
        print('outfile:', outfile)
        print('cpus:', cpus)
        print('smina_commands:', smina_commands)
        
    return cpus, ligfile, outfile, smina_commands

def main():
    '''
    Main script
    '''
    CPUS, ligfile, outfile, smina_commands = parse_args()
    
    if DEBUG > 1:
        print('Splitting input into', CPUS, 'parts')
    splits = splitligs(ligfile, CPUS)
    
    mult_run_smina(splits, outfile, smina_commands)
    
    if DEBUG > 1:
        print('Molecules in input file:', TOTAL_COMPOUNDS)
        with open(outfile, 'rb') as out_file:
            print('Molecules in results file:', out_file.read().count(b'$$$$'))
        
    if CLEANUP:
        if DEBUG > 0:
            print('Cleaning up', len(temp_files), 'temp files')
        clean_up(temp_files)
    
if __name__ == "__main__":
    main()
