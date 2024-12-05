#
#   plot_guess                       
#                               
#   Script to plot an initial flowfield guess created using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_guess.py casename"

# Import modules and functions
from routines import *

def main():
    nblocks = {"naca":5, "bump":1, "bend":1, "turbine_c":5, "turbine_h":1}
    # Construct full filenames to read the guess data
    filename = 'out_guess_' + sys.argv[-1] + '.bin'
    inname = 'input_' + sys.argv[-1] + '.txt'

    # Read the case from file
    av = read_settings(inname)
    g = read_case(filename)
    
    blocks = []
    for n in range(nblocks[sys.argv[-1]]):
        b = cut_block(g, n)
        b = calc_secondary(av, b)
        blocks.append(b)

    # Open figure window and open four subplots
    fig,ax = plt.subplots(2,2,sharex=True,sharey=True,figsize=[14.4,7.2]); 
    fig.tight_layout()

    # Set subplot aspect ratios as equal and remove axes labels
    ax = ax.flatten()
    for a in ax:
        a.set_aspect('equal',adjustable='box'); a.axis('off')

    # Plot the primary flow variables to show the guess
    fieldnames = ['ro','roe','rovx','rovy']
    for n,name in enumerate(fieldnames):
        for nn in range(len(blocks)):
        	# Plot filled contour levels
        	hc = ax[n].pcolormesh(blocks[nn]['x'],blocks[nn]['y'],blocks[nn][name],shading='gouraud')

  		# Add colorbar with variable name
        	colorbar(hc,name)

        	# Draw the walls of the block
        	plot_wall(ax[n],blocks[nn])

    # Show all the plots
    plt.show()

    
main()


