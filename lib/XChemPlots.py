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

        category = (    'Crystals mounted',
                        'Data Collection',
                        'Maps',
                        'PANDDA',
                        'Refinement',
                        'Comp Chem'             )


        N = len(category)
        Success = np.array([    db_dict['nr_samples'],
                                db_dict['nr_data_collection_success'],
                                3,
                                4,
                                5,
                                6   ])

        Pending = np.array([    0,
                                db_dict['nr_data_collection_pending'],
                                4,
                                3,
                                2,
                                1   ])

        Failure = np.array([    db_dict['nr_samples_failed_to_mount'],
                                db_dict['nr_data_collection_failed'],
                                2,
                                3,
                                2,
                                3   ])

#        Success = np.array([1, 2, 3, 4, 5, 6])
#        Pending = np.array([6, 5, 4, 3, 2, 1])
#        Failure = np.array([2, 3, 2, 3, 2, 3])


        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        p1 = self.overview_axes.bar(ind, Success, width, color='g')
        p2 = self.overview_axes.bar(ind, Pending, width, color='y', bottom=Success)
        p3 = self.overview_axes.bar(ind, Failure, width, color='r', bottom=Pending+Success)


        # add some text for labels, title and axes ticks
        self.overview_axes.set_ylabel('N')
        self.overview_axes.set_title('Overview')
        self.overview_axes.set_xticks(ind + width)
        self.overview_axes.set_xticklabels(category)

        self.overview_axes.legend((p1[0], p2[0], p3[0]), ('Success', 'Pending', 'Failure'))
