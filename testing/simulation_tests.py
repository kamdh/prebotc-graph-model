import numpy as np 
import sys
import unittest

#Must change the directory name depending on directory location for runmodel.py
sys.path.append('/Users/mendoza/Documents/prebotc-graph-model/model/')
import runmodel

#Must change the directory name depending on directory location for doPost.py 
sys.path.append('/Users/mendoza/Documents/prebotc-graph-model/postprocessing_preBotBot/')
import doPost

'''
    This class will test the simulations to ensure that our data is not making large changes when tightening the 
    absolute and relative tolerance of the integrators
'''
class TestSimulations(unittest.TestCase):

    '''If you are passing in any file, make sure you are using the full path
       For testing you can pass in placeholder names for the following paramters:
            runmodel.py: output
            doPost.py: sim,output
        This is because we are just passing the arguments in through the main function instead of them being saved
      sim_param_list: String of parameters seperated by a space 
      post_param_list: String of parameters seperated by a space
    '''
    def __init__(self,testname,sim_param_list, abs_tol, rel_tol, post_param_list):
        super(TestSimulations,self).__init__(testname)
        self.sim_params = sim_param_list.split()
        self.rel_tol = rel_tol
        self.abs_tol = abs_tol
        self.post_params = post_param_list.split()


    '''
        Checks to see if the two integrators were able to generate metrics that were within 10 percent of each other to ensure no drastic metric changes are observed
    '''
    def test_simulation(self):
        #Creates a seperate parameter list, which will have modified integration tolerances
        
        sim_params_newTol = self.sim_params + ['-a',abs_tol,'-r',rel_tol]
        
	
    
        print 'Should Run'
        def_tol_sim = runmodel.main(self.sim_params)
        new_tol_sim = runmodel.main(sim_params_newTol)
        
        (simFn, outFn, trans, sec_flag, spike_thresh, f_sigma, butter_low,
         butter_high, bin_width, cutoff, are_volts) = doPost.parse_args(self.post_params)

        post_data_def = doPost.run(def_tol_sim,trans,sec_flag,spike_thresh,f_sigma,butter_low,butter_high,bin_width,cutoff,are_volts)
        post_data_new = doPost.run(def_tol_sim,trans,sec_flag,spike_thresh,f_sigma,butter_low,butter_high,bin_width,cutoff,are_volts)

        def_lag = post_data_def['mean_phi']
        new_lag = post_data_new['mean_phi']
        lag_diff = abs(def_lag - new_lag) / def_lag
        self.assertTrue(lag_diff < .1)
        def_pop_corr = post_data_def['chi1']
        new_pop_corr = post_data_new['chi1']
        pop_corr_diff = abs(def_pop_corr-new_pop_corr) / def_pop_corr
        self.assertTrue(pop_corr_diff < .1)


	

if __name__ == '__main__':
    params = sys.argv[1]
    abs_tol = sys.argv[2]
    rel_tol = sys.argv[3]
    post_params = sys.argv[4]
    output = sys.arg[5]
    log_file = output
    f = open(log_file,"w")
    suite = unittest.TestSuite()
    suite.addTest(TestSimulations("test_simulation",params,abs_tol,rel_tol,post_params))
    unittest.TextTestRunner(f,verbosity=2).run(suite)
    f.close()
