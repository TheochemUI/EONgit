import os
import logging
logger = logging.getLogger('superbasinlist')

import superbasin

class SuperbasinScheme:
    ''' This poorly-named class handles keeping track of which states belong
    to which superbasins, the SuperBasin object of those superbasins, and 
    superbasining criteria. It also expands and merges superbasins'''

    def __init__(self, superbasin_path, states):
        self.path = superbasin_path
        self.states = states

        if not os.path.isdir(self.path):
            logger.warning('superbasin path does not exist, creating %s' % self.path)
            os.makedirs(self.path)
        
        self.superbasins = []
        dirs = os.listdir(self.path)
        self.next_sb_num = 0
        for i in dirs:
            self.next_sb_num = max(self.next_sb_num, int(i))
            path = os.path.join(self.path, i)
            self.superbasins.append(superbasin.Superbasin(path, get_state = states.get_state))
        self.next_sb_num += 1
        self.read_data()

    def get_containing_superbasin(self, state):
        for i in self.superbasins:
            if i.contains_state(state):
                return i
        return None

    def register_transition(self, start_state, end_state):
        raise NotImplementedError()

    def write_data(self):
        raise NotImplementedError()
    
    def read_data(self):
        raise NotImplementedError()

    def __del__(self):
        self.write_data()
    #def get_superbasins(self):
    #    if self.superbasins is not None:
    #        return self.superbasins
    #    else:
    #        dirs = os.listdir(self.path)
    #        for i in dirs:
    #            path = os.path.join(self.path, i)
    #            self.superbasins[int(i)] = superbasin.Superbasin(path)


class TransitionCounting(SuperbasinScheme):
    ''' Implements the transition counting scheme for superbasining '''
    
    def __init__(self, superbasin_path, states, num_transitions):
        self.num_transitions = num_transitions
        SuperbasinScheme.__init__(self,superbasin_path, states)

    def register_transition(self, start_state, end_state):
        logger.debug('Registering transitions')
        if start_state not in self.count:
            self.count[start_state] = {}
            table_path = os.path.join(start_state.path)
            if os.path.isfile(table_path):
                f = open(table_path, 'r')
                for i in f:
                    i = i.strip().split()
                    self.count[start_state][self.states.get_state(int(i[0]))] = i[1]
                f.close()
        if end_state not in self.count[start_state]:
            self.count[start_state][end_state] = 0
        self.count[start_state][end_state] += 1
        
        if self.count[start_state][end_state] >= self.num_transitions:
            logger.debug( "Making basin....")
            self.make_basin([start_state, end_state])

    
    def make_basin(self, merge_states):
        new_sb_states = []
        for i in merge_states:
            sb = self.get_containing_superbasin(i)
            if sb is None:
                new_sb_states.append(i)
            else:
                new_sb_states += sb.states
                sb.delete()
                self.superbasins.remove(sb)
        new_sb_path = os.path.join(self.path, str(self.next_sb_num))
        self.superbasins.append(superbasin.Superbasin(new_sb_path, new_sb_states)) 
        
        logger.info("Created superbasin with states " + str([i.number for i in new_sb_states]))
        self.next_sb_num += 1

    def write_data(self):
        logger.debug('writing')
        for start_state in self.count:
            data_path = os.path.join(start_state.path, 'superbasin_tc')
            f = open(data_path, 'w')
            for end_state in self.count[start_state]:
                print >> f, end_state.number, self.count[start_state][end_state]
            f.close()

    def read_data(self):
        self.count = {}
        logger.debug('reading')
        for i in range(self.states.get_num_states()):
            #TODO: Should each scheme have a unique file?
            state = self.states.get_state(i)
            data_path = os.path.join(state.path, 'superbasin_tc')
            if os.path.isfile(data_path):
                f = open(data_path, 'r')
                self.count[state] = {}
                for i in f:
                    i = i.strip().split()
                    self.count[state][self.states.get_state(int(i[0]))] = int(i[1]) 
                f.close()





#class JeanClaudeThingy(SuperbasinScheme): 

        
        

