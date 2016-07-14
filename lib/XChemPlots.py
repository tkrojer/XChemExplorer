import os,sys

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

import XChemMain

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np


class summary_plot(object):
    def __init__(self,datasource,overview_axes):
        self.datasource=datasource
        self.overview_axes=overview_axes

    def update_overview(self):

        db_dict=XChemMain.get_datasource_summary(self.datasource)

        columns_in_plot = [ 'Crystals mounted',
                            'Data Collection',
                            'Maps',
                            'PANDDA',
                            'Refinement',
                            'Comp Chem'             ]


#        self.overview_figure, self.overview_axes = plt.subplots()
#        ax0, ax1, ax2, ax3 = self.overview_axes.flat

#        in example
#        fig, ax = plt.subplots()

        N = 5
        menMeans = (20, 35, 30, 35, 27)
        menStd = (2, 3, 4, 1, 2)

        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        rects1 = self.overview_axes.bar(ind, menMeans, width, color='r', yerr=menStd)

        womenMeans = (25, 32, 34, 20, 25)
        womenStd = (3, 5, 2, 3, 3)
        rects2 = self.overview_axes.bar(ind + width, womenMeans, width, color='y', yerr=womenStd)

        # add some text for labels, title and axes ticks
        self.overview_axes.set_ylabel('Scores')
        self.overview_axes.set_title('Scores by group and gender')
        self.overview_axes.set_xticks(ind + width)
        self.overview_axes.set_xticklabels(('G1', 'G2', 'G3', 'G4', 'G5'))

#        ax0 = self.overview_axes.flat
#        ax.plot

#        ax0.set_title('Data Collection - Outcome')
#        ax0.set_ylabel("n")
#        outcome=[]
#        frequency=[]
#        for key in db_dict:
#            if key.startswith('no_data_collection'):
#                outcome.append(key.replace('no_data_collection_',''))
#                frequency.append(db_dict[key])
#        y_pos = np.arange(len(outcome))
#        try:
#            ax0.bar(y_pos, frequency, width=0.15)
#        except ValueError:
#            pass
#        ax0.set_xticks(np.arange(len(outcome)) + 0.15/2)
#        ax0.set_xticklabels(outcome, rotation=0)
#        print 'outcome',outcome
#        print 'frequency',frequency


#        N = 5
#        menMeans = (20, 35, 30, 35, 27)
#        womenMeans = (25, 32, 34, 20, 25)
#        menStd = (2, 3, 4, 1, 2)
#        womenStd = (3, 5, 2, 3, 3)
#        ind = np.arange(N)    # the x locations for the groups
#        width = 0.35       # the width of the bars: can also be len(x) sequence
#
#        p1 = plt.bar(ind, menMeans, width, color='r', yerr=menStd)
#        p2 = plt.bar(ind, womenMeans, width, color='y',
#                bottom=menMeans, yerr=womenStd)
#
#        plt.ylabel('Scores')
#        plt.title('Scores by group and gender')
#        plt.xticks(ind + width/2., ('G1', 'G2', 'G3', 'G4', 'G5'))
#        plt.yticks(np.arange(0, 81, 10))
#        plt.legend((p1[0], p2[0]), ('Men', 'Women'))
#
#        plt.show()

