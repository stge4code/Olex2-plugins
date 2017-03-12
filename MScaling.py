from olexFunctions import OlexFunctions

OV = OlexFunctions()

import os
import htmlTools
import olex
import olx
import olex_core
import re

instance_path = OV.DataDir()

p_path = os.path.dirname(os.path.abspath(__file__))

l = open(os.sep.join([p_path, 'def.txt'])).readlines()
d = {}
for line in l:
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    d[line.split("=")[0].strip()] = line.split("=")[1].strip()

p_name = d['p_name']
p_htm = d['p_htm']
p_img = eval(d['p_img'])
p_scope = d['p_scope']

OV.SetVar('MScaling_plugin_path', p_path)

from PluginTools import PluginTools as PT


class Atom():
    def __init__(self, cx=0, cy=0, cz=0, fx=0, fy=0, fz=0, id=0, name='', au=False, rbond=-1):
        self.cx = cx
        self.cy = cy
        self.cz = cz
        self.fx = fx
        self.fy = fy
        self.fz = fz
        self.id = id
        self.name = name
        self.au = au
        self.rbond = rbond

    def fcords(self):
        return [self.fx, self.fy, self.fz]

    def findID(self, ids):
        if '_' not in self.name:
            for id in ids:
                if self.name == ids[id].name:
                    self.id = id
                    self.au = True

    def fShiftAndScale(self, center, scale):
        self.fx = (self.fx - center.fx) * scale + center.fx
        self.fy = (self.fy - center.fy) * scale + center.fy
        self.fz = (self.fz - center.fz) * scale + center.fz
        return self


    def fApply(self):
        olx.xf.au.SetAtomCrd(self.id, self.fx, self.fy, self.fz)
        # olex.f("xf.au.SetAtomCrd({},{},{},{})".format(self.id, self.fx, self.fy, self.fz))
        return self

    def bondsScale(self, rbonds_scale):
        self.rbond  = self.rbond * rbonds_scale
        OV.cmd("conn {} {}".format(self.rbond, self.name))
        return self

    def __str__(self):
        return "{},{}:({},{},{})".format(self.id, self.name, self.cx, self.cy, self.cz)

    def __repr__(self):
        return "{},{},{}:({},{},{})".format(self.id, self.name, self.au, self.fx, self.fy, self.fz)


class MScaling(PT):
    def __init__(self):
        super(MScaling, self).__init__()
        self.p_name = p_name
        self.p_path = p_path
        self.p_scope = p_scope
        self.p_htm = p_htm
        self.p_img = p_img
        self.deal_with_phil(operation='read')
        self.print_version_date()
        self.setup_gui()
        self.params = OV.GuiParams()
        OV.registerFunction(self.grow, True, "MScaling")
        OV.registerFunction(self.shrink, True, "MScaling")
        OV.registerFunction(self.elongate, True, "MScaling")
        OV.registerFunction(self.shorten, True, "MScaling")
        OV.registerFunction(self.fuse, True, "MScaling")
        OV.registerFunction(self.complete_fragments, True, "MScaling")
        OV.registerFunction(self.set_delta_scale, True, "MScaling")
        OV.registerFunction(self.set_delta_rbonds_scale, True, "MScaling")

        self.cell = self.get_cell()
        self.atoms_pool = {}
        #self.delta_scale = 0.25
        #self.delta_rbonds_scale = 0.25
        self.rbonds = {
            'Q': 0.77,
            'C': 0.77,
            'Co': 1.25,
            'P': 1.1,
            'Cl': 0.99,
            'H': 0.32,
            'Br': 1.14,
            'N': 0.77,
            'Cu': 1.28,
            'S': 1.03,
            'Na': 2.65,
            'Ti': 1.45,
            'Fe': 1.24,
            'V': 1.31,
            'O': 0.66,
            'F': 0.64,
        }

    def set_delta_scale(self, scale):
        if ',' in scale:
            scale = scale.replace(',', '.')
        try:
            if scale == '':
                return
            float(scale)
        except(SyntaxError, NameError, ValueError):
            print('Invalid value for scale provided')
            return
        OV.SetParam('mscaling.delta_scale', scale)
        #print "Scale factor was changed to {}".format(OV.GetParam('mscaling.delta_scale'))

    def set_delta_rbonds_scale(self, rbonds_scale):
        if ',' in rbonds_scale:
            rbonds_scale = rbonds_scale.replace(',', '.')
        try:
            if rbonds_scale == '':
                return
            float(rbonds_scale)
        except(SyntaxError, NameError, ValueError):
            print('Invalid value for scale provided')
            return
        OV.SetParam('mscaling.delta_rbonds_scale', rbonds_scale)
        #print "Bonds scale factor was changed to {}".format(OV.GetParam('mscaling.delta_rbonds_scale'))

    def near(self, tmpf, tmpidf, n=3):

        '''result = True
        if round(tmpf[0], n) != round(tmpidf[0], n): result = result and False
        if round(tmpf[1], n) != round(tmpidf[1], n): result = result and False
        if round(tmpf[2], n) != round(tmpidf[2], n): result = result and False
        '''
        '''
        eps = 0.001
        if abs(tmpf[0] - tmpidf[0]) > eps:
            return False
        elif abs(tmpf[1] - tmpidf[1]) > eps:
            return False
        elif abs(tmpf[2] - tmpidf[2]) > eps:
            return False
        else:
            return True
        '''
        if round(tmpf[0], n) != round(tmpidf[0], n):
            return False
        elif round(tmpf[1], n) != round(tmpidf[1], n):
            return False
        elif round(tmpf[2], n) != round(tmpidf[2], n):
            return False
        else:
            return True

    def findID(self, name, ids):
        for id in ids:
            if name == ids[id][0]: return id
        return 0

    def testAuByName(self):
        if '_' not in self.name:
            return True
        return False

    def get_atoms_list(self, part=None, notype=''):
        model = olex_core.GetRefinementModel(False)
        asym_unit = model['aunit']
        atoms = {}
        for residue in asym_unit['residues']:
            for atom in residue['atoms']:
                if atom['type'] == notype:
                    continue
                if part:
                    if atom['part'] != part:
                        continue
                try:
                    resnum = residue['number']
                except KeyError:
                    resnum = 0
                atoms[atom['aunit_id']] = [atom['label'], atom['crd'][0], atom['part'], resnum, atom['type']]
        return atoms

    def get_rbond(self, name):
        tmptype = re.search(r'^([^0-9]+)*', name).group(0)
        if tmptype in self.rbonds:
            return self.rbonds[tmptype]
        else:
            #out = olex.f("info({})".format(name))
            #print len(out)
            #print re.search(r'Atom[*]+\n([*]+)\n', olex.f("info({})".format(name))).group(0)
            return self.rbonds['Q']

    def rescale(self, scale, rbonds_scale):
        self.refresh()
        atoms_names = olex.f("sel(a)").split()
        if not atoms_names:
            print 'No atoms selected!'
            return
        atoms = []
        ids = self.get_atoms_list()
        tmpf = [float(item) for item in olex.f("ccrd()").split()]
        center = Atom(fx=tmpf[0], fy=tmpf[1], fz=tmpf[2], name="center")
        OV.cmd("sel -u")

        for name in atoms_names:
            if '_' not in name:
                if name not in self.atoms_pool:
                    id = self.findID(name, ids)
                    #tmpf = [float(x) for x in (olx.xf.au.GetAtomCrd(id)).split(',')]
                    tmpf = ids[id][1]
                    self.atoms_pool[name] = Atom(fx=tmpf[0], fy=tmpf[1], fz=tmpf[2], name=name, au=True, id=id, rbond=self.get_rbond(name))
                    atoms.append(self.atoms_pool[name])
                else:
                    atoms.append(self.atoms_pool[name])
        for atom in atoms: atom.fShiftAndScale(center, scale).fApply()

        for atom in atoms: atom.bondsScale(rbonds_scale)
        OV.cmd("fuse")
        OV.cmd("grow")
        return

    def adjust_bonds(self, rbonds_scale):
        self.refresh()
        atoms_names = olex.f("sel(a)").split()
        if not atoms_names:
            print 'No atoms selected!'
            return
        ids = self.get_atoms_list()
        atoms = []
        for name in atoms_names:
            if '_' not in name:
                if name not in self.atoms_pool:
                    id = self.findID(name, ids)
                    #tmpf = [float(x) for x in (olx.xf.au.GetAtomCrd(id)).split(',')]
                    tmpf = ids[id][1]
                    self.atoms_pool[name] = Atom(fx=tmpf[0], fy=tmpf[1], fz=tmpf[2], name=name, au=True, id=id,
                                                 rbond=self.get_rbond(name))
                    atoms.append(self.atoms_pool[name])
                else:
                    atoms.append(self.atoms_pool[name])
        for atom in atoms: atom.bondsScale(rbonds_scale)
        OV.cmd("fuse")
        OV.cmd("grow")

    def grow(self):
        self.rescale(1.0 + OV.GetParam('mscaling.delta_scale'), 1.0 + OV.GetParam('mscaling.delta_rbonds_scale'))
    def shrink(self):
        self.rescale(1.0 - OV.GetParam('mscaling.delta_scale'), 1.0 - OV.GetParam('mscaling.delta_rbonds_scale'))
    def elongate(self):
        self.adjust_bonds(1.0 + OV.GetParam('mscaling.delta_rbonds_scale'))
    def shorten(self):
        self.adjust_bonds(1.0 - OV.GetParam('mscaling.delta_rbonds_scale'))
    def fuse(self):
        OV.cmd("fuse")
    def complete_fragments(self):
        OV.cmd("grow")
    def get_cell(self):
        precell = olex_core.GetRefinementModel(False)['aunit']['cell']
        cell = [precell['a'][0], precell['b'][0], precell['c'][0], precell['alpha'][0],
                precell['beta'][0], precell['gamma'][0]]
        return cell
    def refresh(self):
        cell = self.get_cell()
        for i in range(6):
            if self.cell[i] != cell[i]:
                self.atoms_pool = {}
                self.cell = cell
                return
MScaling_instance = MScaling()
