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
#        if os.path.isfile(self.settings['data_source']):
#            sqlite_query = ('Select '
#                            'CrystalName,DataCollectionOutcome,RefinementOutcome,DataProcessingResolutionHigh,RefinementRfree '
#                            'FROM mainTable')
#            query=XChemDB.data_source(self.settings['data_source']).execute_statement(sqlite_query)
#            Rfree_present=[]
#            Rfree_missing=[]
#            Resolution_present=[]
#            Resolution_missing=[]
#            data_collection_outcome={}
#            refinement_outcome={}
#
#
#            for item in query:
#
#                if not str(item[1]).replace('Failed - ','') in data_collection_outcome:
#                    if str(item[1])=='None':
#                        data_collection_outcome['pending']=1
#                    else:
#                        data_collection_outcome[str(item[1]).replace('Failed - ','')]=1
#                else:
#                    if str(item[1])=='None':
#                        data_collection_outcome['pending']+=1
#                    else:
#                        data_collection_outcome[str(item[1]).replace('Failed - ','')]+=1
#
#                if not str(item[2]) in refinement_outcome:
#                    if str(item[2]) != 'None':
#                        refinement_outcome[str(item[2])]=1
#                else:
#                    refinement_outcome[str(item[2])]+=1
#
#                if str(item[3])=='None':
#                    Resolution_missing.append(0)
#                else:
##                    if isinstance(float(item[3]),float):
#                    try:
#                        Resolution_present.append(float(item[3]))
#                    except ValueError:
#                        pass
#
#                if str(item[4])=='None':
#                    Rfree_missing.append(0)
#                else:
#                    try:
#                        Rfree_present.append(float(item[4]))
#                    except ValueError:
#                        pass

        db_dict=XChemMain.get_datasource_summary(self.datasource)


        ax0, ax1, ax2, ax3 = self.overview_axes.flat

        ax0.set_title('Data Collection - Outcome')
        ax0.set_ylabel("Frequency")
        outcome=[]
        frequency=[]
        for key in db_dict:
            if key.startswith('no_data_collection'):
                outcome.append(key.replace('no_data_collection_',''))
                frequency.append(db_dict[key])
            y_pos = np.arange(len(outcome))
            try:
                ax0.bar(y_pos, frequency, width=0.15)
            except ValueError:
                pass
            ax0.set_xticks(np.arange(len(outcome)) + 0.15/2)
            ax0.set_xticklabels(outcome, rotation=0)

#            ax1.set_title('Data Collection - Resolution')
#            try:
#                ax1.hist((Resolution_missing, Resolution_present), bins=20, color=("red", "green"), label=("missing","analysed"))
#            except ValueError:
#                pass
#            ax1.set_xlabel("Resolution")
#            ax1.legend(prop={'size': 10})
##
#            ax2.set_title('Map Analysis - Outcome')
#            ax2.set_ylabel("Frequency")
#            outcome=[]
#            frequency=[]
#            for key in refinement_outcome:
#                outcome.append(key)
#                frequency.append(refinement_outcome[key])
#            y_pos = np.arange(len(outcome))
#            try:
#                ax2.bar(y_pos, frequency, width=0.15)
#            except ValueError:
#                pass
#            ax2.set_xticks(np.arange(len(outcome)) + 0.15/2)
#            ax2.set_xticklabels(outcome, rotation=0)
#
#            ax3.set_title('Refinement - Rfree')
#            try:
#                ax3.hist((Rfree_missing, Rfree_present), bins=20, color=("red", "green"), label=("missing","analysed"))
#            except ValueError:
#                pass
#            ax3.set_xlabel("Rfree")
#            ax3.legend(prop={'size': 10})
#
#
#            ax = self.overview_figure.add_subplot(111)
#            ax.hist((Rfree_missing, Rfree_present), bins=20, color=("red", "green"), label=("missing","analysed"))
#            ax.legend(prop={'size': 10})
#            ax.set_title('bars with legend')
#            ax.set_xlabel("Rfree")
#            ax.set_ylabel("Frequency")
#        self.overview_canvas.draw()
