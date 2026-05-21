# import numpy as np 
# import matplotlib.pyplot as plt
# import meshio
import pyvista as pv
import pandas as pd

# def plot_error(mesh, results, analytical_sol): 
#
# def analytical_sol(x,y): 
#     return 0.5 * (x**2 - 1) * (y**2 - 1)

# def output_vtu(mesh,solution):
    # points = mesh.points
    # cells = mesh.cells
    # point_data = solution

    # mesh = meshio.Mesh(points = points, cells=cells, point_data=point_data)
    # mesh.write('../output/solution_2d_poisson_dirichlet.vtu')

if __name__ == '__main__':
    print('Reading in Mesh')
    mesh = pv.read('../data/meshes/quader5.vtk')
    print('Reading in solution')
    solution = pd.read_csv('../output/poisson_bpx.csv')
    mesh['u'] = solution['u'].values
    print('Saving mesh')
    output_file = '../output/solution_bpx.vtk'
    mesh.save(output_file)
    print(f'VTK file at {output_file}')


