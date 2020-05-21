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
        v['radial_point_txt'] = ''
        v['radial_axis_txt'] = ''
        v['polar_center_txt'] = ''
        v['use_radial'] = False
        v['use_polar'] = False
        v['use_normalp'] = False
        v['normalp_dest'] = ''
        v['reverse_dir_normalp'] = False
        v['reverse_dir_polar'] = False
        v['reverse_dir_radial'] = False        

        return v

    def panel_param(self, obj, validator=None):
        ret = VtableElement.panel_param(self, obj, validator=validator)
        ret[3] = list(ret[3])
        ret[3][0]['choices'].append('By two points')
        ret[3][0]['choices'].append('Normal')
        ret[3][0]['choices'].append('Normal to point')
        ret[3][0]['choices'].append('Radial')
        ret[3][0]['choices'].append('Polar')
        elp3 = [["Point(from)",  None, 0,  {}],
                ["Point(to)",    None, 0,  {}],
                [None, True, 3,  {"text": "Use distance between points"}],
                [None, True, 3,  {"text": "Reverse direction"}], ]
        elp4 = [[None, True, 3,  {"text": "Reverse direction"}],
                [None,  "(Note) Face must be flat-plane", 2,  {}], ]
        elp5 = [["Point(to)",    None, 0,  {}],
                [None, True, 3,  {"text": "Reverse direction"}],
                [None,  "(Note) Face must be flat-plane. Length below is mulitplier", 2,  {}], ]
        elp6 = [["Axis Dir.",    None, 0,  {}],
                ["Point on axis",    None, 0,  {}],
                [None, True, 3,  {"text": "Reverse direction"}], ]
        elp7 = [["Center",    None, 0,  {}], 
                [None, True, 3,  {"text": "Reverse direction"}], ]        

        ret[3].append({'elp': elp3})
        ret[3].append({'elp': elp4})
        ret[3].append({"elp": elp5})
        ret[3].append({"elp": elp6})
        ret[3].append({"elp": elp7})

        return ret

    def get_panel_value(self, obj):
        ret = VtableElement.get_panel_value(self, obj)
        if obj.use_normal:
            ret[0] = 'Normal'
        elif obj.use_fromto_points:
            ret[0] = 'By two points'
        elif obj.use_normalp:
            ret[0] = 'Normal to point'
        elif obj.use_radial:
            ret[0] = 'Radial'
        elif obj.use_polar:
            ret[0] = 'Polar'
        else:
            pass

        ret.append([obj.fromto_points_txt1,
                    obj.fromto_points_txt2,
                    obj.fromto_points_use_dist,
                    obj.reverse_dir_fromto])
        ret.append([obj.reverse_dir_normal,
                    '(Note) Face must be flat-plane', ])
        ret.append([obj.normalp_dest,
                    obj.reverse_dir_normalp,
                    '(Note) Face must be flat-plane. Distance below is ignored', ])
        ret.append([obj.radial_axis_txt,
                    obj.radial_point_txt, 
                    obj.reverse_dir_radial, ])        
        ret.append([obj.polar_center_txt, 
                    obj.reverse_dir_polar, ])

        return ret

    def import_panel_value(self, obj, v):
        obj.use_normal = False
        obj.use_normalp = False
        obj.use_radial = False
        obj.use_polar = False
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
        elif v[0] == 'Radial':
            setattr(obj, 'use_m_'+self.name, False)           
            obj.use_radial = True
            obj.radial_axis_txt = v[6][0]
            obj.radial_point_txt = v[6][1]
            obj.reverse_dir_radial = v[6][2]            
        elif v[0] == 'Polar':
            setattr(obj, 'use_m_'+self.name, False)
            obj.use_polar = True
            obj.polar_center_txt = v[7][0]
            obj.reverse_dir_polar = v[7][1]
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
        elif obj.use_radial:
            ret = ['radial',
                   obj.radial_axis_txt,
                   obj.radial_point_txt,
                   obj.reverse_dir_radial]            
        elif obj.use_polar:
            ret = ['polar',
                   obj.polar_center_txt, 
                   obj.reverse_dir_polar]
        else:
            pass
        return ret


class VtableElement_Plane(VtableElement):
    def add_attribute(self, v):
        v = VtableElement.add_attribute(self, v)
        v['use_face_parallel'] = False
        v['face_parallel_txt1'] = ''
        v['face_parallel_txt2'] = ''
        v['use_3_points'] = False
        v['by_3_points_txt1'] = ''
        v['by_3_points_txt2'] = ''
        v['by_3_points_txt3'] = ''

        return v

    def panel_param(self, obj, validator=None):
        ret = VtableElement.panel_param(self, obj, validator=validator)
        ret[3] = list(ret[3])
        ret[3][0]['choices'].append('By 3 points')
        ret[3][0]['choices'].append('Face parallel and point')
        elp3 = [["Point 1",  None, 0,  {}],
                ["Point 2",    None, 0,  {}],
                ["Point 3",    None, 0,  {}], ]
        elp4 = [["Face",  None, 0,  {}],
                ["Point",  None, 0,  {}], ]
        ret[3].append({'elp': elp3})
        ret[3].append({'elp': elp4})

        return ret

    def get_panel_value(self, obj):
        ret = VtableElement.get_panel_value(self, obj)
        if obj.use_face_parallel:
            ret[0] = 'Face parallel and point'
        elif obj.use_3_points:
            ret[0] = 'By 3 points'
        else:
            pass

        ret.append([obj.by_3_points_txt1,
                    obj.by_3_points_txt2,
                    obj.by_3_points_txt3, ])
        ret.append([obj.face_parallel_txt1,
                    obj.face_parallel_txt2, ])

        return ret

    def import_panel_value(self, obj, v):
        obj.use_3_points = False
        obj.use_face_parallel = False
        if v[0] == 'Face parallel and point':
            setattr(obj, 'use_m_'+self.name, False)
            obj.use_face_parallel = True
            obj.face_parallel_txt1 = v[4][0]
            obj.face_parallel_txt2 = v[4][1]
        elif v[0] == 'By 3 points':
            setattr(obj, 'use_m_'+self.name, False)
            obj.use_3_points = True
            obj.by_3_points_txt1 = v[3][0]
            obj.by_3_points_txt2 = v[3][1]
            obj.by_3_points_txt3 = v[3][2]
        else:
            VtableElement.import_panel_value(self, obj, v)

    def make_value_or_expression(self, obj):
        ret = VtableElement.make_value_or_expression(self, obj)
        if obj.use_face_parallel:
            ret = ['face_parallel',
                   obj.face_parallel_txt1,
                   obj.face_parallel_txt2]
        elif obj.use_3_points:
            ret = ['3_points',
                   obj.by_3_points_txt1,
                   obj.by_3_points_txt2,
                   obj.by_3_points_txt3, ]
        else:
            ret = ['by_abc', ret]
        return ret


class VtableElement_ToPoint(VtableElement):
    pass
