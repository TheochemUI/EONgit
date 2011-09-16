import logging
logger = logging.getLogger('explorer')
from time import time
import shutil
import StringIO
import os
import sys
import cPickle as pickle

import atoms
import communicator
import config
import displace
import fileio as io
import kdb
import recycling

def get_minmodexplorer():
    if config.akmc_server_side:
        return ServerMinModeExplorer
    else:
        return ClientMinModeExplorer

class Explorer:
    def __init__(self):
        self.wuid_path = os.path.join(config.path_scratch, "wuid")
        self.load_wuid()

    def load_wuid(self):
        try:
            f = open(self.wuid_path)
            self.wuid = int(f.read())
            f.close()
        except IOError:
            self.wuid = 0

    def save_wuid(self):
        f = open(self.wuid_path, 'w')
        f.write("%i\n" % self.wuid)
        f.close()



class MinModeExplorer(Explorer):
    def __init__(self, states, previous_state, state):
        Explorer.__init__(self)
        self.states = states
        self.state = state
        self.previous_state = previous_state
        self.comm = communicator.get_communicator()

        if config.recycling_on: 
            self.nrecycled = 0
            self.recycler = recycling.Recycling(self.states,
                                           self.previous_state, 
                                           self.state, 
                                           config.recycling_move_distance,
                                           config.recycling_save_sugg)

        # If we plan to only displace atoms that moved getting to the current state.
        if config.disp_moved_only and self.state.number != 0:
            moved_atoms = self.recycler.process_atoms
        else:
            moved_atoms = None


        if config.kdb_on:
            self.kdber = kdb.KDB()
            if len(self.state.get_ratetable()) <= 1:
                self.kdber.query(self.state, wait = config.kdb_wait)


        self.reactant = self.state.get_reactant()
        self.displace = displace.DisplacementManager(self.reactant, moved_atoms)

    def explore(self):
        self.register_results()
        if self.state.get_confidence() < config.akmc_confidence:
            self.make_jobs()
        else:
            num_cancelled = self.comm.cancel_state(self.state.number)
            logger.info("cancelled %i workunits from state %i", 
                        num_cancelled, self.state.number)
            #XXX: Do we ever call explore on a completed state twice?
            if config.kdb_on:
                logger.info("Adding relevant processes to kinetic database.")
                for process_id in self.state.get_process_ids():
                    output = self.kdber.add_process(self.state, process_id)
                    logger.debug("kdbaddpr.pl: %s" % output)

    def generate_displacement(self):
        if config.recycling_on and self.state.number is not 0:
            displacement, mode = self.recycler.make_suggestion()
            if displacement:
                self.nrecycled += 1
                return displacement, mode, 'recycling'

        if config.kdb_on:
            displacement, mode = self.kdber.make_suggestion()
            if displacement:
                logger.info('Made a KDB suggestion')
                return displacement, mode, 'kdb'

        displacement, mode = self.displace.make_displacement() 
        return displacement, mode, 'random'


class ClientMinModeExplorer(MinModeExplorer):
    def __init__(self, states, previous_state, state):
        MinModeExplorer.__init__(self, states, previous_state, state)
        job_table_path = os.path.join(config.path_root, "jobs.tbl")
        job_table_columns = [ 'state', 'wuid', 'type']
        self.job_table = io.Table(job_table_path, job_table_columns)

    def make_jobs(self):
        #XXX:what if the user changes the bundle size?
        num_in_buffer = self.comm.get_queue_size()*config.comm_job_bundle_size 
        logger.info("%i searches in the queue" % num_in_buffer)
        num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
        logger.info("making %i searches" % num_to_make)
        
        if num_to_make == 0:
            return
        
        searches = []
        
        invariants = {}

        reactIO = StringIO.StringIO()
        io.savecon(reactIO, self.reactant)
        invariants['reactant_passed.con']=reactIO
        
        ini_changes = [ ('Main', 'job', 'process_search') ]
        invariants['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)

        #Merge potential files into invariants
        invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

        t1 = time()
        for i in range(num_to_make):
            search = {}
            # The search dictionary contains the following key-value pairs:
            # id - CurrentState_WUID
            # displacement - an atoms object containing the point the saddle search will start at
            # mode - an Nx3 numpy array containing the initial mode 
            search['id'] = "%d_%d" % (self.state.number, self.wuid)
            displacement, mode, disp_type = self.generate_displacement()
            self.job_table.add_row( {'state':self.state.number,
                                     'wuid':self.wuid,
                                     'type':disp_type } )

            if displacement:
                dispIO = StringIO.StringIO()
                io.savecon(dispIO, displacement)
                search['displacement_passed.con'] = dispIO
                modeIO = StringIO.StringIO()
                io.save_mode(modeIO, mode)
                search['mode_passed.dat'] = modeIO
                searches.append(search) 
                self.wuid += 1
                # eager write
                self.save_wuid()

        if config.recycling_on and self.nrecycled > 0:
            logger.info("recycled %i saddles" % self.nrecycled)

        try:
            self.comm.submit_jobs(searches, invariants)
            t2 = time()
            logger.info( str(len(searches)) + " searches created") 
            logger.debug( str(num_to_make/(t2-t1)) + " searches per second")
        except:
            logger.exception("Failed to submit searches.")
        self.job_table.write()

    def register_results(self):
        logger.info("registering results")
        t1 = time()
        if os.path.isdir(config.path_jobs_in):
            shutil.rmtree(config.path_jobs_in)  
        os.makedirs(config.path_jobs_in)
        
        #Function used by communicator to determine whether to discard a result
        def keep_result(name):
            # note that all processes are assigned to the current state 
            state_num = int(name.split("_")[0])
            return (state_num == self.state.number and \
                    self.state.get_confidence() < config.akmc_confidence)

            #return (config.debug_register_extra_results or \
            #        state_num == self.state.number or \
            #        self.state.get_confidence() < config.akmc_confidence)

        num_registered = 0
        for result in self.comm.get_results(config.path_jobs_in, keep_result): 
            # The result dictionary contains the following key-value pairs:
            # reactant - an array of strings containing the reactant
            # saddle - an atoms object containing the saddle
            # product - an array of strings containing the product
            # mode - an array of strings conatining the mode
            # results - a dictionary containing the key-value pairs in results.dat
            # id - StateNumber_WUID
            #
            # The reactant, product, and mode are passed as lines of the files because
            # the information contained in them is not needed for registering results
            if config.debug_keep_all_results:
                #XXX: We should only do these checks once to speed things up, 
                #     but at the same time debug options don't have to be fast
                #save_path = os.path.join(config.path_root, "old_searches")
                #if not os.path.isdir(save_path):
                #    os.mkdir(save_path)
                #shutil.copytree(result_path, os.path.join(save_path, i))
                #XXX: This is currently broken by the new result passing 
                #      scheme. Should it be done in communicator?
                pass
            state_num = int(result['name'].split("_")[0])
            id = int(result['name'].split("_")[1]) + result['number']
            searchdata_id = "%d_%d" % (state_num, id)
            # Store information about the search into result_data for the 
            # search_results.txt file in the state directory.
            job_type = self.job_table.get_row('wuid', id)['type']
            result['type'] = job_type
            if job_type == None:
                logger.warning("Could not find search data for search %s" 
                               % searchdata_id)
            else:
                self.job_table.delete_row('wuid', id)
            result['wuid'] = id
            
            #read in the results
            result['results'] = io.parse_results(result['results.dat'])
            if result['results']['termination_reason'] == 0:
                self.state.add_process(result)
            else:
                self.state.register_bad_saddle(result, config.debug_keep_bad_saddles)
            num_registered += 1
            
            if self.state.get_confidence() >= config.akmc_confidence:
                if not config.debug_register_extra_results:
                    break
        
        #Approximate number of searches recieved
        tot_searches = len(os.listdir(config.path_jobs_in)) * config.comm_job_bundle_size
        
        t2 = time()
        logger.info("%i (result) searches processed", num_registered)
        logger.info("Approximately %i (result) searches discarded." % (tot_searches - num_registered))
        if num_registered == 0:
            logger.debug("0 results per second", num_registered)
        else:
            logger.debug("%.1f results per second", (num_registered/(t2-t1)))
            
        self.job_table.write()
        return num_registered

class ServerMinModeExplorer(MinModeExplorer):
    def __init__(self, states, previous_state, state):
        MinModeExplorer.__init__(self, states, previous_state, state)
        #XXX: need to init somehow
        self.search_id = 0

        self.job_queue = []
        self.wuid_to_search_id = {}
        self.process_searches = {}

        if os.path.isfile("explorer.pickle"):
            f = open("explorer.pickle")
            tmp_dict = pickle.load(f)
            f.close()
            self.__dict__.update(tmp_dict)


    def save(self):
        f = open("explorer.pickle", "w")
        d = self.__dict__.copy()
        del d['states']
        del d['previous_state']
        del d['state']
        del d['comm']
        pickle.dump(d, f)
        f.close()

    def explore(self):
        if self.state.get_energy() != None:
            MinModeExplorer.explore(self)
        else:
            if self.comm.get_queue_size() == 0:
                logger.info("submitting initial state minimization job")
                job = {}
                reactant = self.state.get_reactant()
                reactIO = StringIO.StringIO()
                io.savecon(reactIO, reactant)
                job['reactant_passed.con'] = reactIO
                ini_changes = [ ('Main', 'job', 'minimization') ]
                job['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)
                job['id'] = '0'
                invariants = io.load_potfiles(config.path_pot)
                self.comm.submit_jobs([job], invariants)

            while self.comm.get_queue_size() != 0:
                continue

            if not os.path.isdir(config.path_jobs_in):
                os.makedirs(config.path_jobs_in)

            for result in self.comm.get_results(config.path_jobs_in, lambda x: True):
                results_dat = io.parse_results(result['results.dat'])
                energy = results_dat['potential_energy']
                reason = results_dat['termination_reason']
                if reason != 0:
                    logger.fatal("minimization of initial reactant failed")
                    sys.exit(1)
                self.state.set_energy(energy)

    def register_results(self):
        logger.info("registering results")
        t1 = time()
        if os.path.isdir(config.path_jobs_in):
            try:
                shutil.rmtree(config.path_jobs_in)  
            except OSError, msg:
                logger.error("error cleaning up %s: %s", config.path_jobs_in, msg)
        os.makedirs(config.path_jobs_in)

        if not os.path.isdir(config.path_incomplete):
            os.makedirs(config.path_incomplete)
        
        #Function used by communicator to determine whether to keep a result
        def keep_result(name):
            # note that all processes are assigned to the current state 
            state_num = int(name.split("_")[0])
            return (state_num == self.state.number and \
                    self.state.get_confidence() < config.akmc_confidence)
     
            #return (config.debug_register_extra_results or \
            #        state_num == self.state.number or \
            #        self.state.get_confidence() < config.akmc_confidence)

        num_registered = 0
        for result in self.comm.get_results(config.path_jobs_in, keep_result): 
            state_num = int(result['name'].split("_")[0])
            # XXX: doesn't this doesn't give the correct id wrt bundling
            id = int(result['name'].split("_")[1]) + result['number']
            searchdata_id = "%d_%d" % (state_num, id)

            #results = io.parse_results(result['results.dat'])
            #job_type = results['job_type']
            #termination_reason = results['termination_reason']
            #reason = self.job_termination_reasons[job_type][termination_reason]

            search_id = self.wuid_to_search_id[id]
            #logger.info("got result for search_id %i" % search_id)
            self.process_searches[search_id].process_result(result)

            num_registered += 1
            if self.state.get_confidence() >= config.akmc_confidence:
                if not config.debug_register_extra_results:
                    break

        #Approximate number of searches recieved
        tot_searches = len(os.listdir(config.path_jobs_in)) * config.comm_job_bundle_size
        
        t2 = time()
        logger.info("%i (result) searches processed", num_registered)
        logger.info("Approximately %i (result) searches discarded." % (tot_searches - num_registered))
        if num_registered == 0:
            logger.debug("0 results per second", num_registered)
        else:
            logger.debug("%.1f results per second", (num_registered/(t2-t1)))
            
        self.save()
        return num_registered

    def make_jobs(self):
        num_in_buffer = self.comm.get_queue_size()*config.comm_job_bundle_size 
        logger.info("%i searches in the queue" % num_in_buffer)
        num_to_make = max(config.comm_job_buffer_size - num_in_buffer, 0)
        logger.info("making %i searches" % num_to_make)
        
        if num_to_make == 0:
            return
        
        jobs = []
        
        invariants = {}

        #Merge potential files into invariants
        invariants = dict(invariants,  **io.load_potfiles(config.path_pot))

        t1 = time()

        #start new searches
        for i in range(num_to_make):
            job = None
            for ps in self.process_searches.values():
                job = ps.get_job()
                if job:
                    self.wuid_to_search_id[self.wuid] = ps.search_id
                    #logger.info("got job for process search %i" % ps.search_id)
                    break

            if not job:
                displacement, mode, disp_type = self.generate_displacement()
                process_search = ProcessSearch(self.state, displacement, mode, disp_type, self.search_id)
                self.process_searches[self.search_id] = (process_search)
                self.wuid_to_search_id[self.wuid] = self.search_id
                job = process_search.get_job()

                self.search_id += 1


            job['id'] = "%i_%i" % (self.state.number, self.wuid)
            self.wuid += 1
            self.save_wuid()
            jobs.append(job)

        self.save()
        if config.recycling_on and self.nrecycled > 0:
            logger.info("recycled %i saddles" % self.nrecycled)

        try:
            self.comm.submit_jobs(jobs, invariants)
            t2 = time()
            logger.info( str(len(jobs)) + " searches created") 
            logger.debug( str(num_to_make/(t2-t1)) + " searches per second")
        except:
            logger.exception("Failed to submit searches.")



class ProcessSearch:
    def __init__(self, state, displacement, mode, displacement_type, search_id):
        self.state = state
        self.displacement = displacement
        self.displacement_type = displacement_type
        self.mode = mode
        self.search_id = search_id

        #valid statuses are 'not_started', 'running', 'complete', 'unneeded', and 'error'
        self.job_statuses = {
                             'saddle_search':'not_started',
                             'min1':'not_started',
                             'min2':'not_started'
                            }
        self.job_iterations = { 
                               'saddle_search':0,
                               'min1':0,
                               'min2':0,
                              }

        unknown = "unknown_exit_code"
        self.job_termination_reasons = {
                'saddle_search':[ "good", unknown, "no_convex", "high_energy",
                                  "max_concave_iterations", 
                                  "max_iterations", unknown, unknown, unknown, 
                                  unknown, unknown, "potential_failed", ],
                'minimization':[ "good", "max_iterations", "potential_failed", ]}

        self.finished_jobs = []

        self.finished_saddle_name = None
        self.finished_min1_name = None
        self.finished_min2_name = None
        self.finished_reactant_name = None
        self.finished_product_name = None

        self.data = { 
                      'termination_reason':0,
                      'potential_energy_saddle':None,
                      'potential_energy_reactant':None,
                      'potential_energy_product':None,
                      'barrier_reactant_to_product':None,
                      'barrier_product_to_reactant':None,
                      'prefactor_reactant_to_product':config.process_search_default_prefactor,
                      'prefactor_product_to_reactant':config.process_search_default_prefactor,
                      'displacement_saddle_distance':0.0,
                      'force_calls_saddle':0,
                      'force_calls_minimization':0,
                      'force_calls_prefactors':0,
                    }

    def get_job(self):
        #not yet sure how to handle error reporting
        if any( [ s == 'error' for s in self.job_statuses.values() ] ):
            return None

        if self.job_statuses['saddle_search'] == 'not_started':
            self.job_statuses['saddle_search'] = 'running'
            return self.start_search()

        if self.job_statuses['min1'] == 'not_started' and \
           self.job_statuses['saddle_search'] == 'complete':
            self.job_statuses['min1'] = 'running'
            return self.start_minimization('min1')

        if self.job_statuses['min2'] == 'not_started' and \
           self.job_statuses['saddle_search'] == 'complete':
            self.job_statuses['min2'] = 'running'
            return self.start_minimization('min2')

        return None

    def process_result(self, result):
        results_dat = io.parse_results(result['results.dat'])
        result['results.dat'].seek(0)
        job_type = results_dat['job_type']
        termination_code = results_dat['termination_reason']
        termination_reason = self.job_termination_reasons[job_type][termination_code]

        self.save_result(result)
        self.finished_jobs.append(result['name'])

        if termination_reason == 'good':
            if job_type == 'saddle_search':
                logger.info("search_id: %i saddle search complete" % self.search_id)
                self.job_statuses[job_type] = 'complete'
                self.finished_saddle_name = result['name']
                self.finish_search(result)

            elif job_type == 'minimization':
                if self.job_statuses['min1'] not in [ 'complete', 'error' ]:
                    logger.info("search_id: %i minimization 1 complete" % self.search_id)
                    self.job_statuses['min1'] = 'complete'
                    self.finished_min1_name = result['name']
                else:
                    logger.info("search_id: %i minimization 2 complete" % self.search_id)
                    self.job_statuses['min2'] = 'complete'
                    self.finished_min2_name = result['name']
                    self.finish_minimization(result)
            elif job_type == 'minimization':
                if self.job_statuses['min1'] == 'running':
                    self.job_statuses['min1'] = 'incomplete'
                else:
                    self.job_statuses['min2'] = 'incomplete'
        else:
            if job_type == 'saddle_search':
                self.data['termination_reason'] = termination_code
            elif job_type == 'minimization':
                self.data['termination_reason'] = 11
            self.register()

        if all( [ s == 'complete' for s in self.job_statuses.values() ] ):
            logger.info("search_id: %i process search complete" % self.search_id)
            self.register()

    def register(self, repeat=False):
        result = {}
        saddle_result = self.load_result(self.finished_saddle_name)

        if self.finished_reactant_name and not repeat:
            reactant_result = self.load_result(self.finished_reactant_name) 
            result['reactant.con'] = reactant_result['reactant.con']
            product_result = self.load_result(self.finished_product_name)
            result['product.con'] = product_result['reactant.con']


        result['saddle.con'] = saddle_result['saddle.con']
        result['mode.dat'] = saddle_result['mode.dat']
        result['results'] = self.data

        result['results.dat'] = StringIO.StringIO(
                '\n'.join([ "%s %s" % (k,v) for k,v in self.data.items() ]) )
        result['type'] = self.displacement_type
        result['search_id'] = self.search_id

        if self.data['termination_reason'] == 0:
            self.state.add_process(result)
        else:
            self.state.register_bad_saddle(result, config.debug_keep_bad_saddles)

    def save_result(self, result):
        dir_path = os.path.join(config.path_incomplete, result['name'])
        os.makedirs(dir_path)
        for k in result:
            if hasattr(result[k], 'getvalue'):
                fn = os.path.join(dir_path, k)
                f = open(fn, "w")
                f.write(result[k].getvalue())
                f.close()

    def load_result(self, result_name):
        dir_path = os.path.join(config.path_incomplete, result_name)
        result = {}
        for file in os.listdir(dir_path):
            file_path = os.path.join(dir_path, file)
            result[file] = StringIO.StringIO(open(file_path).read())
        return result

    def start_minimization(self, which_min):
        job = {}

        saddle_path = os.path.join(config.path_incomplete, self.finished_saddle_name)

        mode = open(os.path.join(saddle_path, "mode.dat"))
        mode = io.load_mode(mode)

        reactant = open(os.path.join(saddle_path, "saddle.con"))
        reactant = io.loadcon(reactant)

        if which_min == "min2":
            mode = -mode

        reactant.r += config.process_search_minimization_offset*mode

        reactIO = StringIO.StringIO()
        io.savecon(reactIO, reactant)
        job['reactant_passed.con'] = reactIO
        
        ini_changes = [ ('Main', 'job', 'minimization') ]
        job['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)

        return job

    def finish_minimization(self, result):
        #result1 = result
        #result2 = self.load_result(self.finished_min1_name)
        result1 = self.load_result(self.finished_min1_name)
        result2 = result

        atoms1 = io.loadcon(result1['reactant.con'])
        atoms2 = io.loadcon(result2['reactant.con'])

        results_dat1 = io.parse_results(result1['results.dat'])
        results_dat2 = io.parse_results(result2['results.dat'])
        self.data['force_calls_minimization'] += results_dat1['total_force_calls']
        self.data['force_calls_minimization'] += results_dat2['total_force_calls']

        is_reactant = lambda a: atoms.match(a, self.state.get_reactant(), 
                                            config.comp_eps_r, 
                                            config.comp_neighbor_cutoff, True)

        if is_reactant(atoms1):
            reactant_results_dat = results_dat1
            product_results_dat = results_dat2
            self.finished_reactant_name = self.finished_min1_name
            self.finished_product_name = self.finished_min2_name
        elif is_reactant(atoms2):
            reactant_results_dat = results_dat2
            product_results_dat = results_dat1
            self.finished_reactant_name = self.finished_min2_name
            self.finished_product_name = self.finished_min1_name
        else:
            #Not connected
            self.data['termination_reason'] = 6
            self.data['potential_energy_saddle'] = 0.0
            self.data['potential_energy_reactant'] = 0.0
            self.data['potential_energy_product'] = 0.0
            self.data['barrier_reactant_to_product'] = 0.0
            self.data['barrier_product_to_reactant'] = 0.0
            self.register()
            return

        self.data['potential_energy_reactant'] = reactant_results_dat['potential_energy']
        self.data['potential_energy_product'] = product_results_dat['potential_energy']

        self.data['barrier_reactant_to_product'] = self.data['potential_energy_saddle'] - \
                self.data['potential_energy_reactant'] 
        self.data['barrier_product_to_reactant'] = self.data['potential_energy_saddle'] - \
                self.data['potential_energy_product'] 


    def start_search(self):
        job = {}

        dispIO = StringIO.StringIO()
        io.savecon(dispIO, self.displacement)
        job['displacement_passed.con'] = dispIO

        modeIO = StringIO.StringIO()
        io.save_mode(modeIO, self.mode)
        job['mode_passed.dat'] = modeIO
             
        reactIO = StringIO.StringIO()
        io.savecon(reactIO, self.state.get_reactant())
        job['reactant_passed.con'] = reactIO

        ini_changes = [ ('Main', 'job', 'saddle_search') ]
        job['config_passed.ini'] = io.modify_config(config.config_path, ini_changes)

        return job

    def finish_search(self, result):
        results_dat = io.parse_results(result['results.dat'])
        self.data.update(results_dat)
        saddle = io.loadcon(result['saddle.con'])
        barrier = results_dat['potential_energy_saddle'] - self.state.get_energy()
        self.data['potential_energy_reactant'] = self.state.get_energy()
        self.data['barrier_reactant_to_product'] = barrier
        repeat_id = self.state.find_repeat(saddle, barrier)
        if repeat_id != None:
            self.job_statuses['min1'] = 'unneeded'
            self.job_statuses['min2'] = 'unneeded'
            self.register(repeat=True)
