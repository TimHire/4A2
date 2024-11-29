#
#   plot_coord
#                               
#   Script to plot a mesh created using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_coord.py casename"

# Import modules and functions
from routines import *

def main():

    nblocks = {"naca":5, "bump":1, "bend":1, "turbine_c":5}

    # Construct full filename to read the grid data
    filename = 'out_coord_' + sys.argv[-1] + '.bin'

    # Read the case from file
    g = read_case(filename)

    blocks = []
    for n in range(nblocks[sys.argv[-1]]):
        b = cut_block(g, n)
        blocks.append(b)

    x_min = float('inf')
    y_min = float('inf')
    x_max = float('-inf')
    y_max = float('-inf')

    for b in blocks:
        x_min = min(x_min, np.min(b['x']))
        x_max = max(x_max, np.max(b['x']))
        y_min = min(y_min, np.min(b['y']))
        y_max = max(y_max, np.max(b['y']))

    # Open figure window and set the axes to be equal
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('x / m'); ax.set_ylabel('y / m');
    ax.set_aspect('equal',adjustable='box'); ax.tick_params(direction='in')
    
    for i in range(nblocks[sys.argv[-1]]):
        # Plot the mesh coordinates to show the cells
        ax.plot(blocks[i]['x'],blocks[i]['y'],color=cols[i,:],linewidth=0.5)
        ax.plot(np.transpose(blocks[i]['x']),np.transpose(blocks[i]['y']),color=cols[i,:],
            linewidth=0.5)

        # Draw the boundary of the block
        #plot_bound(ax,blocks[i])

    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_min,y_max)
    # Show all the plots
    plt.show()

    
main()


