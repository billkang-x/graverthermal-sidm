import os,sys
import time
import SourcePy.evolve as evolve
import SourcePy.record as record

def get_parameters(run_name):
    """
    Function to get user-defined halo and run parameters for new evolution runs,
    given an input run name. Add new sets of run parameters here.

    Parameters
    ----------
    run_name: string
        Name of halo evolution run.

    Returns
    -------
    Tuple of dictionaries containing halo and run parameters that differ from
    their default values defined in `Halo` class.
    """
    if run_name=='default':
        # use default halo parameters
        hdict = {}
        # adjust run parameters from their defaults for illustrative example
        rdict = {'t_end': 25,'save_frequency_rate': 10}
    else:
        raise IOError('Run name {} is not recognized.'.format(run_name))
    return hdict,rdict

def perform_run(run_name,dir_output='Data'):
    """
    Perform halo evolution run.

    Parameters
    ----------
    run_name: string
        Name of halo evolution run, with halo and run parameters determined
        from output of `get_parameters()`.

    dir_output: string, default: 'Data'
        Top-level directory for all halo evolution runs. Output files from each
        run are saved in a subdirectory matching the run name.

    Returns
    -------
    `Halo` class instance created for the evolution run.
    """
    # get halo and run parameters
    hdict,rdict = get_parameters(run_name)

    # set data directory
    dir_data = os.path.join(dir_output,run_name)

    # create halo record
    halorec = record.HaloRecord(dir_data)

    # create halo
    haloevo = evolve.Halo(halorec,**hdict)

    # evolve halo
    start=time.time()
    haloevo.evolve_halo(**rdict)
    end=time.time()

    print('time elapsed for run {} = {}'.format(run_name,end-start))
    return halorec,haloevo

###############################################################################
if __name__=="__main__":
    halorec,haloevo = perform_run('default')
