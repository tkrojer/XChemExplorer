import gobject
import sys
import os
import glob
import pickle
import subprocess

# SParkle libraries
sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
sys.path.append(os.getenv('CCP4')+'/lib/python2.7/site-packages')

import XChemDB
import XChemRefine

import SPutils
#import SPmodel
#import SPrefine
#import SPdatabase
#import SPvalidate
#import prepare4deposition
#import SPProcess

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

        # read in settings file from XChemExplorer to set the relevant paths
        self.settings = pickle.load(open("XChemExplorer_settings.pkl","rb"))

        print self.settings
#        self.external_software=self.settings['external_software']
        self.refine_model_directory=self.settings['refine_model_directory']
        self.database_directory=self.settings['database_directory']
        self.data_source_file=self.settings['data_source']

        self.db=XChemDB.data_source(self.data_source_file)


        self.selection_criteria =   [   'All Datasets',
                                        'Analysis Pending',
                                        'Datasets under refinement',
                                        'ligand confirmed',
                                        'finished models'   ]
        # this decides which samples will be looked at
        self.selection_mode = ''


        # Todo list:
        # 0: sampleID
        # 1: compoundID


        # the Folder is kind of a legacy thing because my inital idea was to have separate folders
        # for Data Processing and Refinement
        self.project_directory = self.settings['initial_model_directory']
#        self.DataPath = DataPath
        self.Serial=0
        self.Refine=None
#        self.Action = Action
#        self.ActionPath = ''    # self.ActionPath is were stuff is actually done; depending on Action it's either the ProjectPath or DataPath
#        self .Todo=[]
#        tmp=[]
#        for folders in glob.glob(self.refine_model_directory+'/*'):
#            sample=folders[folders.rfind('/')+1:folders.rfind('.')]
#            if os.path.isfile(os.path.join(self.refine_model_directory,sample,refine.pdb)):
#                tmp.append(os.path.join(self.refine_model_directory,sample,refine.pdb)
#            if os.path.isfile(os.path.join(self.refine_model_directory,sample,refine.mtz)):
#                tmp.append(os.path.join(self.refine_model_directory,sample,refine.mtz)


#        self.User = User
#        self.Password = Password
        self.index = -1
        self.Todo=[]
#        self.Database = SPdatabase.initDB().whichDatabase()

        self.xtalID=''
        self.compoundID=''
#        self.datasetOutcome=''

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
        self.RefmacParams={ 'HKLIN': '', 'HKLOUT': '',
                            'XYZIN': '', 'XYZOUT': '',
                            'LIBIN': '', 'LIBOUT': '',
                            'TLSIN': '', 'TLSOUT': '',
                            'TLSADD': '',
                            'NCYCLES': '10',
                            'MATRIX_WEIGHT': 'AUTO',
                            'BREF':   '    bref ISOT\n',
                            'TLS':    '',
                            'NCS':    '',
                            'TWIN':   ''    }



    def StartGUI(self):

        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        self.window.set_default_size(400, 1200)
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
#        self.RRfreeBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("gray"))

        self.ResolutionLabel = gtk.Label('Resolution')
        self.ResolutionValue = gtk.Label(self.QualityIndicators['Resolution'])
        self.ResolutionBox = gtk.EventBox()
        self.ResolutionBox.add(self.ResolutionValue)
#        self.ResolutionBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("gray"))

        self.MolprobityScoreLabel = gtk.Label('MolprobityScore')
        self.MolprobityScoreValue = gtk.Label(self.QualityIndicators['MolprobityScore'])
        self.MolprobityScoreBox = gtk.EventBox()
        self.MolprobityScoreBox.add(self.MolprobityScoreValue)
#        self.MolprobityScoreBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("gray"))

        self.RamachandranOutliersLabel = gtk.Label('Ramachandran Outliers')
        self.RamachandranOutliersValue = gtk.Label(self.QualityIndicators['RamachandranOutliers'])
        self.RamachandranOutliersBox = gtk.EventBox()
        self.RamachandranOutliersBox.add(self.RamachandranOutliersValue)
#        self.RamachandranOutliersBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("gray"))

        self.RamachandranFavoredLabel = gtk.Label('Ramachandran Favored')
        self.RamachandranFavoredValue = gtk.Label(self.QualityIndicators['RamachandranFavored'])
        self.RamachandranFavoredBox = gtk.EventBox()
        self.RamachandranFavoredBox.add(self.RamachandranFavoredValue)
#        self.RamachandranFavoredBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("gray"))

        self.LigandCCLabel = gtk.Label('Ligand CC')
        self.LigandCCValue = gtk.Label(self.QualityIndicators['LigandCC'])
        self.LigandCCBox = gtk.EventBox()
        self.LigandCCBox.add(self.LigandCCValue)
#        self.LigandCCBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("gray"))

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
        self.pic = gtk.gdk.pixbuf_new_from_file(os.getenv('XChemExplorer_DIR')+'/image/NO_COMPOUND_IMAGE_AVAILABLE.png')
        self.image = gtk.Image()
        self.image.set_from_pixbuf(self.pic)
        self.vbox.pack_start(self.image)

        # define main buttons
        self.PREVbutton = gtk.Button(label="<<<")
        self.NEXTbutton = gtk.Button(label=">>>")
#        self.PREPAREforREFINEMENTbutton = gtk.Button(label="Prepare for Refinement")
#        self.NODATACOLLECTIONFAILEDbutton = gtk.Button(label="Data Collection Failed")
#        self.NODATAPROCESSINGFAILEDbutton = gtk.Button(label="Data Processing Failed")
        self.NOREFINEMENTFAILEDbutton = gtk.Button(label="Refinement Failed")
        self.NOLIGANDFAILEDbutton = gtk.Button(label="No Ligand")
        self.REFINEbutton = gtk.Button(label="Refine")
#        self.PDBREDObutton = gtk.Button(label="Refine @ PDBREDO")
#        self.MODELbutton = gtk.Button(label="create/update Model")
        self.VALIDATEbutton = gtk.Button(label="validate structure")
        self.DEPOSITbutton = gtk.Button(label="prepare for deposition")
#        self.LigandConfidenceButton = gtk.Button(label="ligand confidence")
        self.RefinementParamsButton = gtk.Button(label="refinement parameters")
#        self.AdjustDataProcessingButton = gtk.Button(label="adjust data processing")

        self.cb = gtk.combo_box_new_text()
        self.cb.connect("changed", self.ChooseXtal)
#        for item in self.Todo:
#            self.cb.append_text('%s' %item[0])
        self.vbox.add(self.cb)

        # Buttons
        self.hbox1 = gtk.HBox()
        self.hbox1.pack_start(self.PREVbutton)
        self.hbox1.pack_start(self.NEXTbutton)
        self.vbox.pack_start(self.hbox1)

        # --- functions that are called when buttons are clicked ---
        self.PREVbutton.connect("clicked", self.ChangeXtal,-1)
        self.NEXTbutton.connect("clicked", self.ChangeXtal,+1)

        # things concerning the dataset outcome

        self.hbox5 = gtk.HBox()
        self.hbox5.pack_start(self.NOREFINEMENTFAILEDbutton)
        self.NOREFINEMENTFAILEDbutton.connect("clicked",self.update_data_source,"Failed - refinement fails")
        self.vbox.pack_start(self.hbox5)

        self.hbox6 = gtk.HBox()
        self.hbox6.pack_start(self.NOLIGANDFAILEDbutton)
        self.NOLIGANDFAILEDbutton.connect("clicked",self.update_data_source,"Failed - no ligand")
        self.vbox.pack_start(self.hbox6)

#        elif self.Action == 'Refine':

        # Refinement History
        self.pic2 = gtk.gdk.pixbuf_new_from_file(os.getenv('XChemExplorer_DIR')+'/image/NO_REFINEMENT_HISTORY_AVAILABLE.png')
        self.image2 = gtk.Image()
        self.image2.set_from_pixbuf(self.pic2)
        self.vbox.pack_start(self.image2)

        self.hbox2a = gtk.HBox()
        self.hbox2a.pack_start(self.REFINEbutton)
        self.REFINEbutton.connect("clicked",self.REFINE)
        self.vbox.pack_start(self.hbox2a)

#            self.hbox2b = gtk.HBox()
#            self.hbox2b.pack_start(self.PDBREDObutton)
#            self.PDBREDObutton.connect("clicked",self.PDBREDO)
#            self.vbox.pack_start(self.hbox2b)

        self.hbox2c = gtk.HBox()
        self.hbox2c.pack_start(self.RefinementParamsButton)
        self.RefinementParamsButton.connect("clicked",self.RefinementParams)
        self.vbox.pack_start(self.hbox2c)

#            self.hbox3a = gtk.HBox()
#            self.hbox3a.pack_start(self.MODELbutton)
#            self.MODELbutton.connect("clicked",self.CreateUpdateModel)
#            self.vbox.pack_start(self.hbox3a)

        self.hbox4 = gtk.HBox()
        self.hbox4.pack_start(self.NOLIGANDFAILEDbutton)
        self.NOLIGANDFAILEDbutton.connect("clicked",self.update_data_source,"Failed - no ligand")
        self.vbox.pack_start(self.hbox4)

#        self.hbox6 = gtk.HBox()
#        self.hbox6.pack_start(self.VALIDATEbutton)
#        self.VALIDATEbutton.connect("clicked",self.Validate)
#        self.vbox.pack_start(self.hbox6)

#            self.hbox7 = gtk.HBox()
#            self.hbox7.pack_start(self.DEPOSITbutton)
#            self.DEPOSITbutton.connect("clicked",self.Prepare4Deposition)
#            self.vbox.pack_start(self.hbox7)


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
        self.xtalID = widget.get_active_text()
        for n,item in enumerate(self.Todo):
            if item[0] == self.xtalID:
                self.index = n
#        print '\n==> cootSP: current xtal = %s' %widget.get_active_text()
        # update current xtal
        self.xtalID=self.Todo[self.index][0]
#        if self.Todo[self.index][2]==None:
#            self.compoundID='None'
#        else:
#            self.compoundID=self.Todo[self.index][1]
#        self.datasetOutcome=self.Todo[self.index][2]
        self.RefreshData()

    def update_data_source(self,widget,data=None):
        print 'hallo'
#        data_source(self.data_source_file).update_data_source_from_coot()

    def RefreshData(self):
        # initialize Refinement library
        self.Refine=XChemRefine.Refine(self.project_directory,self.xtalID,self.compoundID)
        self.Serial=self.Refine.GetSerial()

        self.QualityIndicators=SPutils.ParseFiles(self.project_directory,self.xtalID).UpdateQualityIndicators()
        PDBin='refine.pdb'
        # if the structure was previously refined, try to read the parameters
        if self.Serial > 1: self.RefmacParams=self.Refine.ParamsFromPreviousCycle(self.Serial-1)
        self.Refine.GetRefinementHistory()
        try:
            self.pic2 = gtk.gdk.pixbuf_new_from_file(os.path.join(self.project_directory,self.xtalID,'RefinementHistory.png'))
            self.image2.set_from_pixbuf(self.pic2)
        except gobject.GError:
            self.pic2 = gtk.gdk.pixbuf_new_from_file(os.getenv('XChemExplorer_DIR')+'/image/NO_REFINEMENT_HISTORY_AVAILABLE.png')
            self.image2.set_from_pixbuf(self.pic2)

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
        coot.handle_read_draw_molecule_with_recentre(os.path.join(self.project_directory,self.xtalID,'refine.pdb'),0)
        for item in coot_utils_XChem.molecule_number_list():
            if coot.molecule_name(item).endswith('refine.pdb'):
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
            coot.auto_read_make_and_draw_maps(os.path.join(self.project_directory,self.xtalID,'refine.mtz'))

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

        try:
            pic = gtk.gdk.pixbuf_new_from_file(os.path.join(self.project_directory,self.xtalID,self.compoundID+'.png'))
            self.image.set_from_pixbuf(pic)
        except gobject.GError:
            pic = gtk.gdk.pixbuf_new_from_file(os.getenv('XChemExplorer_DIR')+'/image/NO_COMPOUND_IMAGE_AVAILABLE.png')
            self.image.set_from_pixbuf(pic)


    def REFINE(self,widget):
        self.Refine.RunRefmac(self.Serial,self.RefmacParams)
    
    def RefinementParams(self,widget):
        print '\n==> cootSP: changing refinement parameters'
        self.RefmacParams=self.Refine.RefinementParams(self.RefmacParams)

#    def Prepare4Deposition(self,widget):
#        print '\n==> cootSP: current model = %s' %SPdatabase.scarab().GetModels(self.Todo[self.index][3])[0][0]
#        print '\n==> cootSP: starting prepare4deposition GUI...'
#        print '\n==> cootSP: ... be patient, we\'re asking many questions...\n'
#        prepare4deposition.Prepare4Deposition(SPdatabase.initDB().GetModels(self.Todo[self.index][3])[0][0])

#    def Validate(self,widget):
#        SPvalidate.Validate('SParkle',self.ProjectPath+'/'+self.xtalID)

    def set_selection_mode(self,widget):
        self.selection_mode = widget.get_active_text()

    def get_samples_to_look_at(self,widget):
        print 'hallo'
        self.Todo=self.db.get_samples_for_coot(self.selection_mode)

        try:
            # get todo_list from datasource
            self.Todo=self.db.get_samples_for_coot(self.selection_mode)
            print self.Todo
        except:
            # if not available, then parse file system
            for samples in glob.glob(self.refine_model_directory+'/*'):
                xtalID=samples[samples.rfind('/')+1:]
                print xtalID
#
                if os.path.isfile(os.path.join(self.refine_model_directory,xtalID,'refine.pdb')):
#                 if os.path.isfile(os.path.join(self.refine_model_directory,xtalID,'*cif')):
                    compoundID=None
                    self.Todo.append([xtalID,compoundID,None])
        for item in self.Todo:
            self.cb.append_text('%s' %item[0])


if __name__=='__main__':
    GUI().StartGUI()