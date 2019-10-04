from petram.phys.vtable import VtableElement, Vtable, Vtable_mixin

'''
Custom  Vtable for geometry
'''


class VtableElement_Direction(VtableElement):

    def add_attribute(self, v):
        v = VtableElement.add_attribute(self, v) 
        v['use_normal'] = False
        v['use_fromto_points'] = False
        v['fromto_points_txt1'] = ''
        v['fromto_points_txt2'] = ''
        v['fromto_points_use_dist'] = False
        v['reverse_dir_normal'] = False
        v['reverse_dir_fromto'] = False
        v['use_normalp'] = False
        v['normalp_dest'] = ''
        v['reverse_dir_normalp'] = False        
        
        return v
    
    def panel_param(self, obj, validator = None):
        ret = VtableElement.panel_param(self, obj, validator = validator)
        ret[3] = list(ret[3])
        ret[3][0]['choices'].append('By two points') 
        ret[3][0]['choices'].append('Normal')
        ret[3][0]['choices'].append('Normal to point')
        elp3 = [["Point(from)",  None, 0,  {}],
                ["Point(to)",    None, 0,  {}],
                [None, True, 3,  {"text": "Use distance between points"}],
                [None, True, 3,  {"text": "Reverse direction"}],]             
        elp4 = [ [None, True, 3,  {"text": "Reverse direction"}],
                 [None,  "(Note) Face must be flat-plane", 2,  {}],]
        elp5 = [ ["Point(to)",    None, 0,  {}],
                 [None, True, 3,  {"text": "Reverse direction"}],                    
                 [None,  "(Note) Face must be flat-plane. Length below is mulitplier", 2,  {}],]                
        ret[3].append({'elp': elp3})
        ret[3].append({'elp': elp4})
        ret[3].append({"elp": elp5})

        return ret    

    def get_panel_value(self, obj):
        ret = VtableElement.get_panel_value(self, obj)        
        if obj.use_normal:
            ret[0] = 'Normal'
        elif obj.use_fromto_points:
            ret[0] = 'By two points'
        elif obj.use_normalp:
            ret[0] = 'Normal to point'
        else:
            pass

        ret.append([obj.fromto_points_txt1,
                    obj.fromto_points_txt2,
                    obj.fromto_points_use_dist,
                    obj.reverse_dir_fromto])        
        ret.append([obj.reverse_dir_normal,
                   '(Note) Face must be flat-plane',])
        ret.append([obj.normalp_dest,
                    obj.reverse_dir_normalp,                    
                   '(Note) Face must be flat-plane. Distance below is ignored',])
        return ret
    
    def import_panel_value(self, obj, v):
        obj.use_normal = False
        obj.use_normalp = False        
        obj.use_fromto_points = False
        if v[0] == 'Normal':
            setattr(obj, 'use_m_'+self.name, False)
            obj.use_normal = True
            obj.reverse_dir_normal = v[4][0]
        elif v[0] == 'Normal to point':
            setattr(obj, 'use_m_'+self.name, False)
            obj.use_normalp = True
            obj.normalp_dest = v[5][0]
            obj.reverse_dir_normalp = v[5][1]            
        elif v[0] == 'By two points':
            setattr(obj, 'use_m_'+self.name, False)            
            obj.use_fromto_points = True
            obj.fromto_points_txt1 = v[3][0]
            obj.fromto_points_txt2 = v[3][1]
            obj.fromto_points_use_dist = v[3][2]
            obj.reverse_dir_fromto = v[3][3]            
        else:
            VtableElement.import_panel_value(self, obj, v)

    def make_value_or_expression(self, obj):
        ret = VtableElement.make_value_or_expression(self, obj)
        if obj.use_normal:
             ret = ['normal', obj.reverse_dir_normal]
        elif obj.use_normalp:
             ret = ['normalp', obj.normalp_dest, obj.reverse_dir_normalp]
        elif obj.use_fromto_points:
             ret = ['fromto_points',
                     obj.fromto_points_txt1,
                     obj.fromto_points_txt2,
                     obj.fromto_points_use_dist,
                     obj.reverse_dir_normal]    
        else:
            pass
        return ret
class VtableElement_ToPoint(VtableElement):
    pass

