import gobject
import sys
import os
import glob
import pickle
import subprocess

#from matplotlib.figure import Figure
#from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas

# last commit: 05/12/2015

# SParkle libraries
sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
sys.path.append(os.getenv('CCP4')+'/lib/python2.7/site-packages')

import XChemDB
import XChemRefine

import XChemUtils

# libraries from COOT
import pygtk, gtk, pango
import coot
# had to adapt the original coot_utils.py file
# otherwise unable to import the original file without complaints about missing modules etc.
# modified file is now in $XChemExplorer_DIR/lib
import coot_utils_XChem


class GUI(object):

    """
    main class which opens the actual GUI
    """

    def __init__(self):

        print '=====>',os.getcwd()
        # read in settings file from XChemExplorer to set the relevant paths
        self.settings = pickle.load(open("XChemExplorer_settings.pkl","rb"))
        self.refine_model_directory=self.settings['refine_model_directory']
        self.database_directory=self.settings['database_directory']
        self.data_source=self.settings['data_source']
        self.db=XChemDB.data_source(self.data_source)

        # checking for external software packages
        self.external_software=XChemUtils.external_software().check()

        self.selection_criteria =   {   'Show All Datasets':                'RefinementPDB_latest is not null',
                                        'Show Analysis Pending Only':       "RefinementOutcome='Analysis Pending'",
                                        'Show Datasets Under Refinement':   "RefinementOutcome='Refinement Ongoing'",
                                        'Show Confirmed Ligands':           "RefinementOutcome='Ligand Confirmed'",
                                        'SHow Final Structures':            "RefinementOutcome='Structure finished'"   }

        # this decides which samples will be looked at
        self.selection_mode = ''

        # the Folder is kind of a legacy thing because my inital idea was to have separate folders
        # for Data Processing and Refinement
        self.project_directory = self.settings['initial_model_directory']
        self.Serial=0
        self.Refine=None
        self.index = -1
        self.Todo=[]

        self.xtalID=''
        self.compoundID=''
#        self.datasetOutcome=''

        self.pdb_style='refine.pdb'
        self.mtz_style='refine.mtz'

        # some COOT settings
        coot.set_map_radius(15)
        coot.set_colour_map_rotation_for_map(0)
        coot.set_colour_map_rotation_on_read_pdb_flag(0)

        self.QualityIndicators = {  'R':                            '-',
                                    'RRfree':                       '-',
                                    'RRfreeColor':                  'gray',
                                    'Resolution':                   'n/a',
                                    'ResolutionColor':              'gray',
                                    'MolprobityScore':              'n/a',
                                    'MolprobityScoreColor':         'gray',
                                    'RamachandranOutliers':         'n/a',
                                    'RamachandranOutliersColor':    'gray',
                                    'RamachandranFavored':          'n/a',
                                    'RamachandranFavoredColor':     'gray',
                                    'LigandCCcolor':                'gray',
                                    'LigandCC':                     'n/a'   }

        # default refmac parameters
        self.RefmacParams={ 'HKLIN': '',        'HKLOUT': '',
                            'XYZIN': '',        'XYZOUT': '',
                            'LIBIN': '',        'LIBOUT': '',
                            'TLSIN': '',        'TLSOUT': '',
                            'TLSADD': '',
                            'NCYCLES': '10',
                            'MATRIX_WEIGHT': 'AUTO',
                            'BREF':   '    bref ISOT\n',
                            'TLS':    '',
                            'NCS':    '',
                            'TWIN':   ''    }

        refinement_cycle=[1,2,3,4,5,6,7,8,9,10]
        Rfree=  [0.1,0.2,0.14,0.12,0.1,0.1,0.1,0.1,0.1,0.1]
        Rcryst= [0.3,0.4,0.14,0.12,0.2,0.1,0.3,0.1,0.2,0.1]
#        self.canvas = FigureCanvas(self.update_plot(refinement_cycle,Rfree,Rcryst))  # a gtk.DrawingArea
#        self.canvas.set_size_request(300, 150)



    def StartGUI(self):

        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        self.window.set_default_size(400, 1000)
        self.window.set_title("XChemExplorer")
        self.vbox = gtk.VBox()

        # choose here which subset of samples should be looked at
        self.vbox.add(gtk.Label('Select Samples'))
        self.hbox_select_samples=gtk.HBox()
        self.cb_select_samples = gtk.combo_box_new_text()
        self.cb_select_samples.connect("changed", self.set_selection_mode)
        for citeria in self.selection_criteria:
            self.cb_select_samples.append_text(citeria)
        self.hbox_select_samples.add(self.cb_select_samples)
        self.select_samples_button = gtk.Button(label="GO")
        self.select_samples_button.connect("clicked",self.get_samples_to_look_at)
        self.hbox_select_samples.add(self.select_samples_button)
        self.vbox.pack_start(self.hbox_select_samples)

        # next comes a section which displays some global quality indicators
        # a combination of labels and textview widgets, arranged in a table

        self.RRfreeLabel = gtk.Label('R/Rfree')
        self.RRfreeValue = gtk.Label(self.QualityIndicators['RRfree'])
        self.RRfreeBox = gtk.EventBox()
        self.RRfreeBox.add(self.RRfreeValue)

        self.ResolutionLabel = gtk.Label('Resolution')
        self.ResolutionValue = gtk.Label(self.QualityIndicators['Resolution'])
        self.ResolutionBox = gtk.EventBox()
        self.ResolutionBox.add(self.ResolutionValue)

        self.MolprobityScoreLabel = gtk.Label('MolprobityScore')
        self.MolprobityScoreValue = gtk.Label(self.QualityIndicators['MolprobityScore'])
        self.MolprobityScoreBox = gtk.EventBox()
        self.MolprobityScoreBox.add(self.MolprobityScoreValue)

        self.RamachandranOutliersLabel = gtk.Label('Ramachandran Outliers')
        self.RamachandranOutliersValue = gtk.Label(self.QualityIndicators['RamachandranOutliers'])
        self.RamachandranOutliersBox = gtk.EventBox()
        self.RamachandranOutliersBox.add(self.RamachandranOutliersValue)

        self.RamachandranFavoredLabel = gtk.Label('Ramachandran Favored')
        self.RamachandranFavoredValue = gtk.Label(self.QualityIndicators['RamachandranFavored'])
        self.RamachandranFavoredBox = gtk.EventBox()
        self.RamachandranFavoredBox.add(self.RamachandranFavoredValue)

        self.LigandCCLabel = gtk.Label('Ligand CC')
        self.LigandCCValue = gtk.Label(self.QualityIndicators['LigandCC'])
        self.LigandCCBox = gtk.EventBox()
        self.LigandCCBox.add(self.LigandCCValue)

        self.table = gtk.Table(4, 2, False)
        self.table.attach(self.RRfreeLabel,                 0, 1, 0, 1)
        self.table.attach(self.RRfreeBox,                   1, 2, 0, 1)
        self.table.attach(self.ResolutionLabel,             0, 1, 1, 2)
        self.table.attach(self.ResolutionBox,               1, 2, 1, 2)
        self.table.attach(self.MolprobityScoreLabel,        0, 1, 2, 3)
        self.table.attach(self.MolprobityScoreBox,          1, 2, 2, 3)
        self.table.attach(self.RamachandranOutliersLabel,   0, 1, 3, 4)
        self.table.attach(self.RamachandranOutliersBox,     1, 2, 3, 4)
        self.table.attach(self.RamachandranFavoredLabel,    0, 1, 4, 5)
        self.table.attach(self.RamachandranFavoredBox,      1, 2, 4, 5)
        self.table.attach(self.LigandCCLabel,               0, 1, 5, 6)
        self.table.attach(self.LigandCCBox,                 1, 2, 5, 6)
        self.vbox.add(self.table)

        # compound picture
        self.pic = gtk.gdk.pixbuf_new_from_file(os.path.join(os.getenv('XChemExplorer_DIR'),'image','NO_COMPOUND_IMAGE_AVAILABLE.png'))
        self.image = gtk.Image()
        self.image.set_from_pixbuf(self.pic)
        self.vbox.pack_start(self.image)

        # define main buttons
        self.PREVbutton = gtk.Button(label="<<<")
        self.NEXTbutton = gtk.Button(label=">>>")
        self.NOREFINEMENTFAILEDbutton = gtk.Button(label="Refinement Failed")
        self.NOLIGANDFAILEDbutton = gtk.Button(label="No Ligand")
        self.REFINEbutton = gtk.Button(label="Refine")
        self.VALIDATEbutton = gtk.Button(label="validate structure")
        self.DEPOSITbutton = gtk.Button(label="prepare for deposition")
#        self.LigandConfidenceButton = gtk.Button(label="ligand confidence")
        self.RefinementParamsButton = gtk.Button(label="refinement parameters")
#        self.AdjustDataProcessingButton = gtk.Button(label="adjust data processing")
        self.merge_ligand_button=gtk.Button(label="Merge Ligand")
        self.place_ligand_here_button=gtk.Button(label="Place Ligand here")
        self.fit_ligand_to_density_button=gtk.Button(label='Fit Ligand to Density')

        self.cb = gtk.combo_box_new_text()
        self.cb.connect("changed", self.ChooseXtal)
        self.vbox.add(self.cb)

        # Buttons
        self.hbox1 = gtk.HBox()
        self.hbox1.pack_start(self.PREVbutton)
        self.hbox1.pack_start(self.NEXTbutton)
        self.vbox.pack_start(self.hbox1)

        # --- functions that are called when buttons are clicked ---
        self.PREVbutton.connect("clicked", self.ChangeXtal,-1)
        self.NEXTbutton.connect("clicked", self.ChangeXtal,+1)

        # --- modeling shortcuts ---
        self.hbox_for_modeling=gtk.HBox()
        self.hbox_for_modeling.add(self.place_ligand_here_button)
        self.place_ligand_here_button.connect("clicked",self.place_ligand_here)
        self.hbox_for_modeling.add(self.merge_ligand_button)
        self.merge_ligand_button.connect("clicked",self.merge_ligand_into_protein)
        self.hbox_for_modeling.add(self.fit_ligand_to_density_button)
        self.fit_ligand_to_density_button.connect("clicked",self.fit_ligand)
        self.vbox.add(self.hbox_for_modeling)

        # --- things concerning the dataset outcome ---
        self.hbox_refinemnt_outcome=gtk.HBox()
        self.NOREFINEMENTFAILEDbutton.connect("clicked",self.update_data_source,"Refinement Failed")
        self.hbox_refinemnt_outcome.add(self.NOREFINEMENTFAILEDbutton)
        self.NOLIGANDFAILEDbutton.connect("clicked",self.update_data_source,"No Ligand Bound")
        self.hbox_refinemnt_outcome.add(self.NOLIGANDFAILEDbutton)
        self.vbox.pack_start(self.hbox_refinemnt_outcome)

        # Refinement History
#        self.pic2 = gtk.gdk.pixbuf_new_from_file(os.getenv('XChemExplorer_DIR')+'/image/NO_REFINEMENT_HISTORY_AVAILABLE.png')
#        self.image2 = gtk.Image()
#        self.image2.set_from_pixbuf(self.pic2)
#        self.vbox.pack_start(self.image2)

        # not sure why, but adding self.canvas to the window makes it forget all other global variables
#        self.vbox.add(self.canvas)


        # --- Refine button ---
        self.REFINEbutton.connect("clicked",self.REFINE)
        self.vbox.pack_start(self.REFINEbutton)

        # --- Refinement parameters ---
        self.RefinementParamsButton.connect("clicked",self.RefinementParams)
        self.vbox.pack_start(self.RefinementParamsButton)

        # --- CANCEL button ---
        self.CANCELbutton = gtk.Button(label="CANCEL")
        self.CANCELbutton.connect("clicked", self.CANCEL)
        self.vbox.add(self.CANCELbutton)

        self.window.add(self.vbox)
        self.window.show_all()

    def CANCEL(self,widget):
        self.window.destroy()


    def ChangeXtal(self,widget,data=None):
        self.index = self.index + data
        if self.index < 0:
            self.index = 0
        if self.index >= len(self.Todo):
            self.index = len(self.Todo)
        self.cb.set_active(self.index)

    def ChooseXtal(self, widget):
        self.xtalID = str(widget.get_active_text())
        for n,item in enumerate(self.Todo):
            if str(item[0]) == self.xtalID:
                self.index = n
        self.xtalID=str(self.Todo[self.index][0])
        if str(self.Todo[self.index][0]) != None:
            self.compoundID=str(self.Todo[self.index][1])
        else:
            self.compoundID=''
        self.RefreshData()

    def update_data_source(self,widget,data=None):              # update and move to next xtal
        outcome_dict={'RefinementOutcome': data}
        self.db.update_data_source(self.xtalID,outcome_dict)
        self.index+=1
        if self.index >= len(self.Todo):
            self.index = len(self.Todo)
        self.cb.set_active(self.index)

    def RefreshData(self):
        # initialize Refinement library
        self.Refine=XChemRefine.Refine(self.project_directory,self.xtalID,self.compoundID)
        self.Serial=self.Refine.GetSerial()

        self.QualityIndicators=XChemUtils.ParseFiles(self.project_directory,self.xtalID).UpdateQualityIndicators()
        # if the structure was previously refined, try to read the parameters
        if self.Serial > 1: self.RefmacParams=self.Refine.ParamsFromPreviousCycle(self.Serial-1)
# disabled for the moment until matplotlib display function is implemented
        print self.Refine.GetRefinementHistory()
#        try:
#            self.pic2 = gtk.gdk.pixbuf_new_from_file(os.path.join(self.project_directory,self.xtalID,'RefinementHistory.png'))
#            self.image2.set_from_pixbuf(self.pic2)
#        except gobject.GError:
#            self.pic2 = gtk.gdk.pixbuf_new_from_file(os.path.join(os.getenv('XChemExplorer_DIR'),'image','NO_REFINEMENT_HISTORY_AVAILABLE.png'))
#            self.image2.set_from_pixbuf(self.pic2)

        # update pdb & maps
        # - get a list of all molecules which are currently opened in COOT
        # - remove all molecules/ maps before loading a new set
        if len(coot_utils_XChem.molecule_number_list()) > 0:
            for item in coot_utils_XChem.molecule_number_list():
#                if not coot.molecule_name(item).endswith('reference.pdb'): coot.close_molecule(item)
                coot.close_molecule(item)
        coot.set_nomenclature_errors_on_read("ignore")
        # read protein molecule after ligand so that this one is the active molecule
        #coot.handle_read_draw_molecule_with_recentre(Pdb[0],0)
        # read ligand pdb+cif
        coot.handle_read_draw_molecule_with_recentre(os.path.join(self.project_directory,self.xtalID,self.pdb_style),0)
        for item in coot_utils_XChem.molecule_number_list():
            if coot.molecule_name(item).endswith(self.pdb_style):
                coot.set_show_symmetry_master(1)    # master switch to show symmetry molecules
                coot.set_show_symmetry_molecule(item,1) # show symm for model
        if os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.compoundID+'.pdb')):
            coot.handle_read_draw_molecule_with_recentre(os.path.join(self.project_directory,self.xtalID,self.compoundID+'.pdb'),0)
            coot.read_cif_dictionary(os.path.join(self.project_directory,self.xtalID,self.compoundID+'.cif'))

        # read fofo maps
        # - read ccp4 map: 0 - 2fofc map, 1 - fofc.map
        if os.path.isfile(os.path.join(self.project_directory,self.xtalID,'2fofc.map')):
            coot.set_default_initial_contour_level_for_difference_map(3)
            coot.handle_read_ccp4_map(os.path.join(self.project_directory,self.xtalID,'fofc.map'),1)
            coot.set_default_initial_contour_level_for_map(1)
            coot.handle_read_ccp4_map(os.path.join(self.project_directory,self.xtalID,'2fofc.map'),0)
            coot.set_last_map_colour(0,0,1)
        else:
            # try to open mtz file with same name as pdb file
            coot.auto_read_make_and_draw_maps(os.path.join(self.project_directory,self.xtalID,self.mtz_style))

        # update Quality Indicator table
        self.RRfreeValue.set_label(self.QualityIndicators['R']+' / '+self.QualityIndicators['RRfree'])
        self.RRfreeBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RRfreeColor']))
        self.ResolutionValue.set_label(self.QualityIndicators['Resolution'])
        self.ResolutionBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['ResolutionColor']))
        self.MolprobityScoreValue.set_label(self.QualityIndicators['MolprobityScore'])
        self.MolprobityScoreBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['MolprobityScoreColor']))
        self.RamachandranOutliersValue.set_label(self.QualityIndicators['RamachandranOutliers'])
        self.RamachandranOutliersBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RamachandranOutliersColor']))
        self.RamachandranFavoredValue.set_label(self.QualityIndicators['RamachandranFavored'])
        self.RamachandranFavoredBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RamachandranFavoredColor']))
        self.LigandCCValue.set_label(self.QualityIndicators['LigandCC'])
        self.LigandCCBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['LigandCCcolor']))

        print 'HHHHHHHHH',self.compoundID
        try:
            pic = gtk.gdk.pixbuf_new_from_file(os.path.join(self.project_directory,self.xtalID,self.compoundID+'.png'))
            self.image.set_from_pixbuf(pic)
        except gobject.GError:
            pic = gtk.gdk.pixbuf_new_from_file(os.path.join(os.getenv('XChemExplorer_DIR'),'image','NO_COMPOUND_IMAGE_AVAILABLE.png'))
            self.image.set_from_pixbuf(pic)


    def REFINE(self,widget):
        outcome_dict={'RefinementOutcome': 'Refinement Ongoing'}
        self.db.update_data_source(self.xtalID,outcome_dict)
        self.Refine.RunRefmac(self.Serial,self.RefmacParams,self.external_software)
        self.index+=1
        if self.index >= len(self.Todo):
            self.index = len(self.Todo)
        self.cb.set_active(self.index)

    
    def RefinementParams(self,widget):
        print '\n==> XCE: changing refinement parameters'
        self.RefmacParams=self.Refine.RefinementParams(self.RefmacParams)

    def set_selection_mode(self,widget):
        for criteria in self.selection_criteria:
            if criteria==widget.get_active_text():
                self.selection_mode=self.selection_criteria[criteria]
                break

    def get_samples_to_look_at(self,widget):
        # first remove old samples if present
        if len(self.Todo) != 0:
            for n,item in enumerate(self.Todo):
                self.cb.remove_text(0)
        self.Todo=[]
        self.Todo=self.db.get_samples_for_coot(self.selection_mode)
        for item in self.Todo:
            print item
            self.cb.append_text('%s' %item[0])

    def update_plot(self,refinement_cycle,Rfree,Rcryst):
        fig = Figure(figsize=(3, 2), dpi=75)
        Plot = fig.add_subplot(111)
        Plot.set_ylim([0,max(Rcryst+Rfree)])
        Plot.set_xlabel('Refinement Cycle')
        Plot.plot(refinement_cycle,Rfree,label='Rfree')
        Plot.plot(refinement_cycle,Rcryst,label='Rcryst')
        Plot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                ncol=2, mode="expand", borderaxespad=0.)
        return fig

    def place_ligand_here(self,widget):
        print 'hallo'
        coot.move_molecule_here(<molecule_number>)

    def merge_ligand_into_protein(self,widget):
        print 'ok'

    def fit_ligand(self,widget):
        print 'fit'

if __name__=='__main__':
    GUI().StartGUI()