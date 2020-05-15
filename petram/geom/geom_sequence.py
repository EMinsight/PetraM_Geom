from __future__ import print_function

import os
import numpy as np
import time
import tempfile

from six.moves.queue import Empty as QueueEmpty

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GeomSequenceOperator')

class GeomSequenceOperator(object):
    def __init__(self):
        self.geom_sequence = []
        self.p = None

    def clean_queue(self):
        self.task_q.close()
        self.q.close()
        self.task_q.cancel_join_thread()
        self.q.cancel_join_thread()
        self.q  = None
        self.task_q  = None        
        
    def add_sequence(self, gui_name, gui_param, geom_name):
        self.geom_sequence.append((gui_name, gui_param, geom_name))


    def run_generator(self, gui,
                      no_mesh=False, finalize=False, filename = '',
                      progressbar = None, create_process = True,
                      process_param = None, start_idx = 0, trash=''):
        
        kwargs = {'PreviewResolution': gui.geom_prev_res,
                  'PreviewAlgorithm': gui.geom_prev_algorithm,
                  'OCCParallel': gui.occ_parallel,
                  'Maxthreads': gui.maxthreads,
                  'SkipFrag': gui.skip_final_frag,
                  'Use1DPreview': gui.use_1d_preview,
                  'UseOCCPreview': gui.use_occ_preview,                  
                  'UseCurvature': gui.use_curvature,
                  'LongEdgeThr': gui.long_edge_thr,
                  'SmallEdgeThr': gui.small_edge_thr,
                  'SmallEdgeSeg': gui.small_edge_seg,
                  'MaxSeg': gui.max_seg}


        p = process_param[0]
        task_q = process_param[1]
        q = process_param[2]
        self.task_q = task_q
        self.q = q

        args = (self.geom_sequence, no_mesh, finalize, filename, start_idx, trash, kwargs)
        
        task_q.put((1, args))
        '''
        p = mp.Process(target = generator,
                       args = (q, self.geom_sequence, no_mesh,
                               finalize, filename, kwargs))
        p.start()
        '''
        logfile = q.get(True)
        dprint1("log file: ", logfile)

        istep = 0
        
        while True:
            try:
                ret = q.get(True, 1)
                if ret[0]: break
                else:
                    dprint1(ret[1])
                if progressbar is not None:
                    istep += 1
                    if istep < progressbar.GetRange():
                        progressbar.Update(istep, newmsg=ret[1])
                    else:
                        print("Goemetry Generator : Step = " + str(istep) + ret[1])
                        
            except QueueEmpty:
                if not p.is_alive():
                    self.clean_queue()
                    if progressbar is not None:                    
                       progressbar.Destroy()
                    assert False, "Child Process Died"
                    break
                time.sleep(1.)                    
                if progressbar is not None:
                    import wx
                    wx.Yield()
                    if progressbar.WasCancelled():
                       if p.is_alive():
                           p.terminate()
                           self.clean_queue()                           
                       progressbar.Destroy()
                       assert False, "Geometry Generation Aborted"
                    
            time.sleep(0.03)

        if ret[1][0] == 'fail':
            p.terminate()
            self.clean_queue()                           
            return False, ret[1][1]
        else:
            self.gui_data, self.objs, brep_file, data, mappings = ret[1]

            if no_mesh:
                ret =  self.gui_data, self.objs, brep_file, None, None, None

            else:
                geom_msh = data[0]
                if geom_msh != '':
                    import gmsh
                    from petram.geom.read_gmsh import read_pts_groups
                    gmsh.open(geom_msh)
                    geom_msh, l, s, v, vcl, esize = data
                    ptx, cells, cell_data = read_pts_groups(gmsh)
                else:
                    geom_msh, l, s, v, vcl, esize, ptx, cells, cell_data = data

                data = ptx, cells, cell_data, l, s, v

                ret = self.gui_data, self.objs, brep_file, data, vcl, esize
            
            return True, ret
