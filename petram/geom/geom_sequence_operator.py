from __future__ import print_function

import os
import numpy as np
import time
import tempfile

import six
if six.PY2:
    import cPickle as pickle
else:
    import pickle

from six.moves.queue import Empty as QueueEmpty

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GeomSequenceOperator')

test_thread = False

class GeomSequenceOperator():
    def __init__(self):
        self.geom_sequence = []
        self._prev_sequence = []

    def __del__(self):
        self.terminate_child()

    def clean_queue(self):
        if self.use_mp:
            self._p.task_q.close()
            self._p.q.close()
            self._p.task_q.cancel_join_thread()
            self._p.q.cancel_join_thread()

    def terminate_child(self):
        if (hasattr(self, "_p") and self._p.is_alive()):
            if self.use_mp:
                self._p.task_q.close()
                self._p.q.close()
                self._p.task_q.cancel_join_thread()
                self._p.q.cancel_join_thread()
                self._p.terminate()
            else:
                self._p.task_q.put((-1, None))
                self._p.task_q.join()
            del self._p

    def create_new_child(self, use_occ, pgb):
        from petram.geom.gmsh_geom_wrapper import (GMSHGeometryGenerator,
                                                   GMSHGeometryGeneratorTH)
        from petram.geom.occ_geom_wrapper import (OCCGeometryGenerator,
                                                  OCCGeometryGeneratorTH)

        self.terminate_child()
        if pgb is None or globals()['test_thread']:
            self.use_mp = False
            if use_occ:
                p = OCCGeometryGeneratorTH()
            else:
                p = GMSHGeometryGeneratorTH()
        else:
            self.use_mp = True
            if use_occ:
                p = OCCGeometryGenerator()
            else:
                p = GMSHGeometryGenerator()

        p.start()
        self._p = p

    def check_create_new_child(self, use_occ, pgb):
        if not hasattr(self, '_prev_sequence'):
            return True

        if len(self.geom_sequence) < len(self._prev_sequence):
            return True
        
        print(self._p.__class__.__name__)
        if pgb is None or globals()['test_thread']:
            if use_occ:
                if self._p.__class__.__name__ != 'OCCGeometryGeneratorTH':
                    return True
            else:
                if self._p.__class__.__name__ != 'GMSHGeometryGeneratorTH':
                   return True
        else:
            if use_occ:
                if self._p.__class__.__name__ != 'OCCGeometryGenerator':
                    return True
            else:
                if self._p.__class__.__name__ != 'GMSHGeometryGenerator':
                   return True

        for k, s in enumerate(self._prev_sequence):
            s_txt1 = pickle.dumps(s)
            s_txt2 = pickle.dumps(self.geom_sequence[k])
            if s_txt1 != s_txt2:
                dprint1("check not passed", s[0])                
                return True
            dprint1("check passed", s[0])
        return False

    def add_sequence(self, gui_name, gui_param, geom_name):
        self.geom_sequence.append((gui_name, gui_param, geom_name))
        
    def clear_sequence(self):
        self.geom_sequence = []        

    def do_run_generator(self, gui, no_mesh=False, finalize=False,
                         filename='', progressbar=None, start_idx=0,
                         trash=''):
        
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

        p = self._p
        task_q = self._p.task_q
        q = self._p.q

        args = (self.geom_sequence, no_mesh, finalize, filename, start_idx, trash, kwargs)

        task_q.put((1, args))

        logfile = q.get(True)
        dprint1("log file: ", logfile)

        istep = 0

        while True:
            try:
                ret = q.get(True, 1)

                if ret[0]:
                    break

                #dprint1(ret[1])
                if progressbar is not None:
                    istep += 1
                    if istep < progressbar.GetRange():
                        progressbar.Update(istep, newmsg=ret[1])
                    else:
                        print("Goemetry Generator : Step = " +
                              str(istep) + ret[1])

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
                    wx.GetApp().Yield()
                    if progressbar.WasCancelled():
                        self.terminate_child()
                        progressbar.Destroy()
                        assert False, "Geometry Generation Aborted"

            time.sleep(0.03)

        if progressbar is not None:
            progressbar.Destroy()

        if ret[1][0] == 'fail':
            self.terminate_child()
            return False, ret[1][1]
        else:
            self.gui_data, self.objs, brep_file, data, mappings = ret[1]

            if no_mesh or data is None:
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


    def run_generator(self, gui, no_mesh=False, finalize=False,
                      filename='', progressbar=None, trash=''):

        if (hasattr(self, "_p") and self._p.is_alive()):
            new_process = self.check_create_new_child(gui.use_occ_preview, progressbar)
        else:
            new_process = True

        if new_process:
            self.terminate_child()
            self.create_new_child(gui.use_occ_preview, progressbar)
            start_idx = 0
        else:
            ll = len(self._prev_sequence)
            start_idx = ll

        self._prev_sequence = self.geom_sequence

        success, dataset = self.do_run_generator(gui, no_mesh=no_mesh,
                                                 finalize=finalize,
                                                 filename=filename,
                                                 progressbar=progressbar,
                                                 start_idx=start_idx,
                                                 trash=trash,)

        if not success:
            print(dataset)  # this is an error message
            # self._p[0].terminate()
            self._prev_sequence = []
            del self._p

        return success, dataset


