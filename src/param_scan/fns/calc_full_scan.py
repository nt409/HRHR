from tqdm import tqdm
import numpy as np

from param_scan.fns.calc_sing_ps_run import SinglePSRun, ScanOutput







class ParameterScan:
    """
    Inputs:
    - config: param scan config detailing runs to do and bounds on params
    - seed: random seed for these runs

    Outputs:
    - ScanOutput object
    """
    def __init__(self, config, seed) -> None:
        self.config = config
        self.seed = seed
        

    def run(self):
        """
        Run random scan over uniform dists and save output
        """

        output = self._get_scan_output()
        
        output.save(self.config, self.seed)




    def _get_scan_output(self):

        np.random.seed(self.seed)
        
        scan_output = ScanOutput()

        N_ITS = self.config["n_iterations"]

        for r_ind in tqdm(range(N_ITS)):

            run_index = self.seed*N_ITS + r_ind

            this_run = SinglePSRun(self.config, run_index)

            scan_output.add_new_output(this_run.output)

        return scan_output
    
