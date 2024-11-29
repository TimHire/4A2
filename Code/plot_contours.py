#
#   plot_contours
#                               
#   Script to plot a converged flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"

# Import modules and functions
from routines import *

def main():

    nblocks = {"naca":5, "bump":1, "bend":1, "turbine_c":5}
    inletblock = {"naca":4, "bump":1, "bend":1, "turbine_c":4}

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'





    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)
    
    blocks = []
    
    # When presenting results all values should be non-dimensionalised. Two
    # variables of interest might be:
    #    1. Static pressure coefficient, (p - p_ref) / (pstag_ref - p_ref)
    #    2. Mach number, v / (ga * rgas * t)**0.5

    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!
    #g = calc_secondary(av,l[i])

    # Use the "cut_i", "mass_av" AND "area_av" functions to calculate the
    # reference pressures at the inlet plane and therefore the static pressure
    # coefficient
    # INSERT

    for n in range(nblocks[sys.argv[-1]]):
        b = cut_block(g, n)
        b = calc_secondary(av, b)
        blocks.append(b)



    c = cut_i(blocks[inletblock[sys.argv[-1]] - 1], 0)         # Select the right block which contains the inlet
    p_ref, l = area_av(c,'p')
    pstag_ref, mass = mass_av(c,'pstag')

    for n in range(len(blocks)):
        cp = (blocks[n]['p'] - p_ref) / (pstag_ref - p_ref)
        cp0 = (blocks[n]['pstag'] - pstag_ref) / (pstag_ref - p_ref)
        blocks[n]['cp'] = cp
        blocks[n]['cp0'] = cp0

    #g['cp'] = (g['p'] - p_ref) / (pstag_ref - p_ref)


    # Specify the parameters to plot
    fieldnames = ['cp', 'mach'];
    colnames = ['Static pressure coefficient','Mach number']


        # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):
        # Open figure window
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();

        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box'); ax.axis('off')

        global_min = float('inf')
        global_max = float('-inf')
        x_min = float('inf')
        y_min = float('inf')
        x_max = float('-inf')
        y_max = float('-inf')

        for b in blocks:
            global_min = min(global_min, np.min(b[name]))
            global_max = max(global_max, np.max(b[name]))
            x_min = min(x_min, np.min(b['x']))
            x_max = max(x_max, np.max(b['x']))
            y_min = min(y_min, np.min(b['y']))
            y_max = max(y_max, np.max(b['y']))


        for nn in range(len(blocks)):
            # Plot filled contour levels
            hc = ax.pcolormesh(blocks[nn]['x'],blocks[nn]['y'],blocks[nn][name],shading='gouraud', vmin=global_min, vmax=global_max)

            # Add Mach = 1 contours
            if name == 'mach':
                ax.contour(blocks[nn]['x'],blocks[nn]['y'],blocks[nn]['mach'],[1.0],colors='w',
                linewidths=0.5)


            # Draw the walls of the block
            plot_wall(ax,blocks[nn])
        # Add colorbar with variable name
        colorbar(hc,colnames[n])
        ax.set_xlim(x_min,x_max)
        ax.set_ylim(y_min,y_max)

    # Show all the plots
    plt.show()

    
main()


