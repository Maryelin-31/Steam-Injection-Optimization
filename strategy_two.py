import numpy as np
from matplotlib import pyplot as pyplot
from matplotlib.cm import ScalarMappable
import os

from zml import *
from zmlx.config import seepage
from zmlx.seepage_mesh.cube import create_xz
from zml import Dfn2
from zmlx.plt.show_dfn2 import show_dfn2
from zmlx.geometry.seg_point_distance import seg_point_distance
from zmlx.geometry.point_distance import point_distance


from zmlx.fluid.ch4 import create as create_ch4
from zmlx.fluid.h2o_gas import create as create_steam
from zmlx.fluid.conf.O2_gas import create as create_o2
from zmlx.fluid.conf.CO2_gas import create as create_co2
from zmlx.fluid.conf.CO_gas import create as create_co
from zmlx.fluid.h2o import create as create_h2o
from zmlx.fluid.conf.c11h24_liq import create as create_lo
from zmlx.fluid.conf.c22h46_liq import create as create_ho
from zmlx.fluid.kerogen import create as create_kerogen
from zmlx.fluid.char import create as create_char

from zmlx.react import decomposition
from zmlx.react import combustion
from zmlx.react import vapor
from zmlx.react.inh import add_inh

from zmlx.kr.create_krf import create_krf
from zml import Interp1, create_dict
from zmlx.utility.PressureController import PressureController
from zmlx.utility.SeepageCellMonitor import SeepageCellMonitor
from zmlx.utility.SaveManager import SaveManager

def workspace(rate_inj, year_inje):

    def create_mesh():
        mesh = create_xz(x_min=0, dx=0.5, x_max=15,
                          y_min=-1, y_max=0,
                          z_min=0, dz=0.5, z_max=45,)
        # count = 0
        # for cell in mesh.cells:
        #     x, y, z = cell.pos
        #     if 15 <= z <= 30:
        #         count += cell.vol
        #         print((count / 1000) / (24 * 3600))
                # print(cell.vol)
        # print(mesh)
        return mesh
    # create_mesh()
    
    def create_fludefs():
        """
        0. gas = [methane, steam]
        1. h20
        2. lo = light oil
        3. ho = heavy oil
        4. sol= kerogen, char
        
        """
        
        gas = Seepage.FluDef(name='gas')
        gas.add_component(create_ch4(name='ch4'))
        gas.add_component(create_steam(name='steam'))
        
        h2o = create_h2o(name='h2o')
        lo = create_lo(name='lo')
        ho = create_ho(name='ho')
        
        sol = Seepage.FluDef(name='sol')
        sol.add_component(create_kerogen(name='kero'))
        sol.add_component(create_char(name='coke'))
        
        return [gas, h2o, lo, ho, sol]
    
    def create_reactions(temp_max=None):
        
        """
        1. Steam phase transitions
           
        2. Kerogen Decomposition:
            a simplification of doi:10.2172/10169154
            Table 1
        """
        
        result = []
        
        r = vapor.create(vap='steam', wat='h2o', temp_max=temp_max)
        result.append(r)
        
        result.append(r)
        r = decomposition.create(left ='kero',
                                 right=[('ho', 0.663), ('lo', 0.076), ('ch4', 0.046), ('coke', 0.215)], 
                                 temp=563.15, heat=161100.0, rate=4.81e-6)
        result.append(r)
        r = decomposition.create(left ='ho',
                                 right=[('lo', 0.438), ('ch4', 0.217), ('coke', 0.345)], 
                                 temp=623.15, heat=219328.0, rate=2.71e-7)
        # When the solid occupies 80% of the total pores, increase the cracking temperature to limit further decomposition (avoid all pores being occupied by solids)
        add_inh(r, sol='sol', liq=None,
                c=[0, 0.8, 1.0],
                t=[0, 0, 1.0e4])
        result.append(r)
        return result
    
    def create_initial():
        """
        create initial field
        """
    
        def get_initial_t(x, y, z):
            """
            the initial temperature
            """
            if 10 <= z <= 35: 
                return 338.0 + 22.15 - 0.0443 * z  
            else:
                return 300 #isotermal
    
        def get_initial_p(x, y, z):
            """
            the initial pressure
            """
            return 15.0e6 + 5e6 - 1e4 * z
    
        def get_perm(x, y, z):
            """
            the initial permeability
            """
            if 10 <= z <= 35: #impermeable
                return 1.0e-15  
            else:
                return 1.0e-17
            
        def get_initial_s(x, y, z):
            """
            the initial saturation ()
            """
            if 10 <= z <= 35:
                s = {'ch4': 0.08, 'steam': 0,
                     'h2o':0.04, 
                     'lo':0.08,
                     'ho':0.2,
                     'kero':0.6, 'coke':0}
                return s
            else:
                s = {'ch4': 1.0, 'steam': 0,
                     'h2o':0, 
                     'lo':0,
                     'ho':0,
                     'kero':0, 'coke':0}
                return s
    
        def get_fai(x, y, z):
            """
            porosity
            """
            if 10 <= z <= 35:
                return 0.43
            else:
                return 0.01
            
    
        def get_denc(x, y, z):
            """
            density * heat capacity
            """
            if 10 <= z <= 35:
                return 2600 * 1000
            else:
                return 2600 * 2000
    
        def get_heat_cond(x, y, z):
            if 10 <= z <= 35:
                return 2.0
            else:
                return 1.0
        
        
        return {'porosity': get_fai, 'pore_modulus': 100e6, 'p': get_initial_p,
                'temperature': get_initial_t,
                'denc': get_denc, 's': get_initial_s,
                'perm': get_perm, 'heat_cond': get_heat_cond, 'dist': 0.01}
    
    
    "Define Model"
    x, y = create_krf(faic=0.02, n=3.0, k_max=100,
                      s_max=2.0, count=500)
    
    gr = Interp1(x=x, y=y)
    
    kw = create_dict(fludefs=create_fludefs(),
                     reactions=create_reactions(),
                     gr=gr,
                     has_solid=False,)
    
    kw.update(**create_initial())
    
    gravity = [0, 0, -10]
    
    kw.update(create_dict(dt_max=3600.0 * 24.0 * 10.0,
                          gravity=gravity, ))
    
    model = seepage.create(mesh=create_mesh(), **kw)
    
    "Relative Permeability"
    
    vs, vk = create_krf(faic=0.05, n=2.0, count=300)
    model.set_kr(saturation=vs, kr=vk)
    
    "attribute"
    ca = seepage.cell_keys(model)
          
    "Injection"
    # rate_inj = 2.7e-06
    pos_inj = (0, 1.0e3, 22.5)
    id_inj  = model.get_nearest_cell(pos=(0, 1.0e3, 22.5)).index
    cell_inj = model.get_cell(id_inj)
    flu = cell_inj.get_fluid(0).get_component(1)
    flu2 = flu.get_copy()
    flu2.set_attr(seepage.flu_keys(model).temperature, 700)
    model.add_injector(cell=cell_inj, fluid_id=[0, 1], flu=flu2,
                        pos=cell_inj.pos, radi=1.0, opers=[(0, rate_inj)])
    
    # for f in cell_inj.faces:
    #     face_inj = model.get_face(f.index)
        # seepage.set_face(face= face_inj, area = 0)
        # print(face_inj)
    
    "Production"
    
    pos_prd = (15, 1e3, 22.5)
    p_prod  = 5e6
    id_prod = model.get_nearest_cell(pos_prd).index
    virtual_cell = seepage.add_cell(model, pos=pos_prd, porosity=1.0e5, pore_modulus=100e6, vol=1.0,
                                    temperature=350,
                                    p=p_prod, 
                                    s=((1.0, 0), 0, 0, 0, (0, 0)))
    seepage.add_face(model, virtual_cell, model.get_cell(id_prod),
                     heat_cond=0, perm=1.0e-14,
                     area=0.0, #the area 0 indicates the well is closed
                     length=1.0)
    pre_ctrl = PressureController(virtual_cell, t=[0, 1e10], p=[p_prod, p_prod])
    monitor  = SeepageCellMonitor(get_t=lambda: seepage.get_time(model), cell=(virtual_cell, pre_ctrl))
    
    
    "Time step Strategy"
    seepage.set_dv_relative(model, 0.5) # The ratio of the distance traveled by each time step to the grid size
    seepage.set_dt(model, 0.01) # initial value for time step
    seepage.set_dt_max(model, 24 * 3600 * 7)  # Maximum value of time step <one week> # Maximum value of time step <one week>
    seepage.set_dt_min(model, 3600)  # Minimum step size is 1 hour
    
    "save"
    
    def cell_mass(cell):
        fluid = []
        total = []
        for i in range(cell.fluid_number):
            total.append(cell.get_fluid(i).mass)
            if cell.get_fluid(i).component_number == 0:
                fluid.append(cell.get_fluid(i).mass)
            else:
                for j in range(cell.get_fluid(i).component_number):
                    fluid.append(cell.get_fluid(i).get_component(j).mass)               
        saturation = [i / sum(total) for i in fluid]
        return saturation
    
    
    def save_wt(path):
        name = os.path.basename(__file__)
        result_folder = os.path.join(os.getcwd(), f'data_{name}', f'Results_rate{rate_inj}_year{year_inje}', 'wt_cells')
        
        if os.path.exists(f'result_folder'):
            import shutil
            shutil.rmtree(f'result_folder')            
        os.makedirs(result_folder, exist_ok=True)
        
        SavePath = os.path.join(result_folder, path)
        with open(SavePath, 'w') as file:
            for cell in model.cells:
                x, y, z = cell.pos
                satu = cell_mass(cell)
                satu_str = ' '.join(str(i) for i in satu)
                file.write(f'{x} {y} {z} {satu_str}\n')
                
    def save_mass(path):
        name = os.path.basename(__file__)
        result_folder = os.path.join(os.getcwd(), f'data_{name}', f'Results_rate{rate_inj}_year{year_inje}', 'Results_cells')
        
        if os.path.exists(f'result_folder'):
            import shutil
            shutil.rmtree(f'result_folder')            
        os.makedirs(result_folder, exist_ok=True)
        
        SavePath = os.path.join(result_folder, path)
        with open(SavePath, 'w') as file:
            for cell in model.cells:
                x, y, z = cell.pos
                temp = cell.get_attr(seepage.cell_keys(model).temperature)
                pres = cell.pre
                file.write(f'{x} {y} {z} '
                            f'{cell.get_fluid(0).get_component(0).mass} {cell.get_fluid(0).get_component(1).mass} '
                            f'{cell.get_fluid(1).mass} '
                            f'{cell.get_fluid(2).mass} '
                            f'{cell.get_fluid(3).mass} '
                            f'{cell.get_fluid(4).get_component(0).mass} {cell.get_fluid(4).get_component(1).mass} '
                            f'{temp} {pres}\n')
    
    def mass(i, j=None):
        """
        input 
        i = fluid
        j = component
        
        output 
        contour plot saturation
        """    
        X = []
        Z = [] 
        total = []
        fluid = []
        for cell in model.cells:
            x, y, z = cell.pos
            X.append(x)
            Z.append(z)
            total.append(cell.fluid_vol)
            if cell.get_fluid(i).component_number == 0 :
                fluid.append(cell.get_fluid(i).vol)
            else:
                fluid.append(cell.get_fluid(i).get_component(j).vol)
        
        
        vol_frac = [ i / j for i, j in zip(fluid, total)]
        vol_frac = vol_frac[:len(vol_frac) - 1]
        vol_frac = np.array(vol_frac)
        vol_frac = np.transpose(vol_frac.reshape(200, 200))
        fig, ax = pyplot.subplots()
        plot = ax.contourf(vol_frac, 20, extent=[0, 100, 0, 100], cmap='coolwarm', antialiased=True)
        ax.set_xlabel('x, m')
        ax.set_ylabel('y, m') 
        ax.set_ylim(100, 0)
        cbar = fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap), ax=ax)
        cbar.set_label('vol_sat')  # Set label for the colorbar
        
    def temperature():
        """
        input 
        i = fluid
        j = component
        
        output 
        contour plot saturation
        """    
        X = []
        Z = [] 
        temp = []
        for cell in model.cells:
            x, y, z = cell.pos
            temp.append(cell.get_attr(seepage.cell_keys(model).temperature))
        
        
        temp = temp[:len(temp) - 1]
        temp = np.array(temp)
        temp = np.transpose(temp.reshape(30, 90))
        fig, ax = pyplot.subplots()
        # plot = ax.contourf(temp, 20, extent=[0, 100, 0, 100], cmap='coolwarm', antialiased=True, vmin=np.min(350), vmax=np.max(1000))
        plot = ax.contourf(temp, 20, extent=[0, 15, 0, 45], cmap='coolwarm', antialiased=True)
        ax.set_xlabel('x, m')
        ax.set_ylabel('y, m') 
        ax.set_ylim(45, 0)
        ax.set_aspect('equal')
        cbar = fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap), ax=ax)
        cbar.set_label('Temperature, K')  # Set label for the colorbar
        
        pyplot.show()
    
    def iterate():
        name = os.path.basename(__file__)
        
        data_path = os.path.join(os.getcwd(), f'data_{name}')
        resu_path = os.path.join(data_path, f'Results_rate{rate_inj}_year{year_inje}')
        os.makedirs(os.path.join(os.getcwd(), data_path, f'Results_rate{rate_inj}_year{year_inje}'), exist_ok=True)
        folder = os.path.join(os.getcwd(), data_path, f'Results_rate{rate_inj}_year{year_inje}')
        
        solver = ConjugateGradientSolver(tolerance=1.0e-20)
        for step in range(100000000):
            seepage.iterate(model, solver=solver)
            pre_ctrl.update(seepage.get_time(model))
            monitor.update(dt=3600.0 * 24.0)
            
            
            time = seepage.get_time(model)
            if time > 3600 * 24 * 365 * year_inje:
                for f in cell_inj.faces:
                    face_inj = model.get_face(f.index)
                    seepage.set_face(face=face_inj, area=0.0)
            if time > 3600 * 24 * 365 * year_inje:
                virtual_face = model.get_face(model.face_number - 1)
                seepage.set_face(face=virtual_face, area= 1.0)
            
            if time > 3600 * 24 * 365 * 20:
                print(f'time Finish = {seepage.get_time(model) / 3600*34*365}')
                break
            path = f'time_{time}.txt'
            if step % 100 == 0:
                # temperature()
                save_mass(path)
                save_wt(path)
                monitor.save(os.path.join(folder, f'prod.txt'))
                print(f'time = {time/ (3600 * 24 * 365)} , step = {step}')
                
                
            
    iterate()        


rate = [2.7e0, 2.7e-2, 2.7e-4, 2.6e-6]
year = [4    , 6    , 8      , 10]

for i in rate:
    for j in year:
        workspace(rate_inj=i, year_inje=j)













    

