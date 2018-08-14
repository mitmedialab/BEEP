from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import matplotlib.pyplot as plt
import sys
import pandas as pd


##########################################
# BEEP - Base Editing Evaluation Program #
##########################################


#Argument: sequence as a string
#Returns: reverse complement of sequence as a string
def rev_comp(seq):
	return str(Seq(seq).reverse_complement())

#Arguments: the trace metadata, the position of the edited base in the trace file, window around edited base to display
#Returns: list of lists with chromatogram values for each base at each position within window 
def get_window_vals(trace, pos, window):
	#stores values of each base within a specified	
	G_vals = []
	A_vals = []
	T_vals = []
	C_vals = []

	#centers the window 
	for i in xrange(pos - (window/2), pos + ((window+3)/2)):
	#Store the chromatogram values for each base at that position
		base_dict = {'G': trace['DATA9'][i], 'A': trace['DATA10'][i], 'T': trace['DATA11'][i], 'C': trace['DATA12'][i]}
		G_vals.append(base_dict['G'])
		A_vals.append(base_dict['A'])
		T_vals.append(base_dict['T'])
		C_vals.append(base_dict['C'])

	return [G_vals, A_vals, T_vals, C_vals]

#Arguments: List of values for each base to plot, name of sample for title
#Returns: Chromatogram of trace values for each base, labelled with title and legend for base colors
def plot_trace(base_val_list, sample_name):
	#Plot together on one plot with Geneious base colors
	g_plot, = plt.plot(base_val_list[0], color = 'y', label = 'G') #G
	a_plot, = plt.plot(base_val_list[1], color = 'r', label = 'A') #A
	t_plot, = plt.plot(base_val_list[2], color = 'g', label = 'T') #T
	c_plot, = plt.plot(base_val_list[3], color = 'b', label = 'C') #C
	plt.legend([a_plot, c_plot, g_plot, t_plot], ['A', 'C', 'G', 'T'])
	plt.title(sample_name)
	return plt

#Arguments: directory of files, name of file, spacer sequence, position of base in spacer, base conversion, control or sample file, from csv or just one file
#Returns: Ratio of edited base value to total value (edited + nonedited values)
def get_ratio(directory, file_name, spacer, base_pos_in_spacer, base_change, samp_type, multiple):
	sample = SeqIO.read('./'+directory+'/'+file_name, 'abi')

	#Channels for bases: G, A, T, C
	channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']
	channelLabels = 'GATC'

	trace = defaultdict(list)

	for c in channels:
		trace[c] = sample.annotations['abif_raw'][c]

	#Base calls for file
	call = sample.annotations['abif_raw']['PBAS1']
	#These are the positions that correspond to each base call in the original trace
	peaks_vals = sample.annotations['abif_raw']['PLOC1']
	#Find the position in the base call of the first base in the spacer
	pos_in_call = call.find(spacer)
	#Check in reverse complement as well
	if (pos_in_call == -1):
		try:
			pos_in_call = call.find(rev_comp(spacer))
			strand = "-"
		except ValueError:
			print("Cannot find spacer" + spacer + "in " + file_name + ".")
	else:
		strand = "+"
	#Go back to find the exact position of the edited base in the original trace 
	if (strand == "+"):
		actual_pos = peaks_vals[pos_in_call+(base_pos_in_spacer - 1)]
	else:
		#make sure we get correct base in orientation
		actual_pos = peaks_vals[pos_in_call+(len(spacer)-base_pos_in_spacer)]
	
	#Plot trace file within specified window,
	window_to_plot = 200
	if(multiple == False):
		if(samp_type == 'ctrl'):
			ctrl_plt = plot_trace(get_window_vals(trace, actual_pos, window_to_plot), 'Control: ' + file_name)
			ctrl_plt.show()
		else:
			samp_plt = plot_trace(get_window_vals(trace, actual_pos, window_to_plot), 'Sample: ' + file_name)
			samp_plt.show()
	#Save plot for each ab1 file in directory
	else:
		new_plt = plot_trace(get_window_vals(trace, actual_pos, window_to_plot), file_name)
		new_plt.savefig('./'+directory+'/'+file_name.split(".")[0]+'.png')
		new_plt.close()

	#Values of bases within small window
	edited_pos_vals = get_window_vals(trace, actual_pos, 5)

	#Max of values of each base within small window to get appropriate editing ratio
	base_dict = {'G': max(edited_pos_vals[0]), 'A': max(edited_pos_vals[1]), 'T': max(edited_pos_vals[2]), 'C': max(edited_pos_vals[3])}

	#Calculate edited base ratio, depending on strand
	if(strand == "+"):
		ratio = float(base_dict[base_change[1]])/(float(base_dict[base_change[0]]) + float(base_dict[base_change[1]]))
	else:
		#For (-) strand, need to get reverse complement of conversion bases
		ratio = float(base_dict[rev_comp(base_change[1])])/(float(base_dict[rev_comp(base_change[0])]) + float(base_dict[rev_comp(base_change[1])]))
	
	return ratio

#Input: directory of files, name of control file, name of sample file, spacer sequence, position of base in spacer, base conversion, from csv (True) or just one file (False)
#Returns: Normalized efficiency (sample ratio multiplied by normalization factor, calculated from control ratio)
def get_efficiency(directory, control_file, sample_file, spacer, be_position, base_change, multiple = False):
	#Call ratio function for control and sample
	ctrl_ratio = get_ratio(directory, control_file, spacer, int(be_position), base_change, 'ctrl', multiple)
	samp_ratio = get_ratio(directory, sample_file, spacer, int(be_position), base_change, 'samp', multiple)
	#Normalization factor
	norm_factor = (1.0 - ctrl_ratio)
	eff = (samp_ratio)*100.0*norm_factor
	efficiency = "%.2f" % eff

	return efficiency

def main():
	#Handle csv input with multiple samples
	#Example: python beep.py ./folder_with_ab1s/be_template.csv
	if(sys.argv[1].endswith(".csv")):
		df = pd.read_csv(sys.argv[1])
		efficiencies = []
		multiple = True
		for index, row in df.iterrows():
			directory = row['Directory']
			control_file = row['Control']
			sample_file = row['Sample']
			spacer = row['Spacer']
			be_position = row['Base Position in Spacer']
			base_change = row['Conversion']
			#Efficiency for each sample
			efficiencies.append(get_efficiency(directory, control_file, sample_file, spacer, be_position, base_change, multiple))
		#Create new column of efficiencies
		df['Efficiency'] = efficiencies

		df.to_csv(sys.argv[1], index = False)

	#One sample
	#Example: python beep.py be3samples PC14-Mine.ab1 PC04-Mine.ab1 GTGTCTGTGTGGGTGAGTGA 5 CT
	else:
		directory = sys.argv[1]
		control_file = sys.argv[2]
		sample_file = sys.argv[3]
		spacer = sys.argv[4] 
		be_position = sys.argv[5]
		base_change = sys.argv[6]

		efficiency = get_efficiency(directory, control_file, sample_file, spacer, be_position, base_change)
		print("")
		print('Efficiency: {}%'.format(efficiency))
		print("")

if __name__ == "__main__": 
	main()



