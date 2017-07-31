# last edited: 28/07/2017 - 15:00

import gobject
import sys
import os
import pickle

from matplotlib.figure import Figure

# XCE libraries
sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
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

        ###########################################################################################
        # read in settings file from XChemExplorer to set the relevant paths
        self.settings = pickle.load(open(".xce_settings.pkl","rb"))
        remote_qsub_submission=self.settings['remote_qsub']
        self.database_directory=self.settings['database_directory']
        self.xce_logfile=self.settings['xce_logfile']
        self.data_source=self.settings['data_source']
        self.db=XChemDB.data_source(self.data_source)

        # checking for external software packages
        self.external_software=XChemUtils.external_software(self.xce_logfile).check()
        self.external_software['qsub_remote']=remote_qsub_submission

        # the Folder is kind of a legacy thing because my inital idea was to have separate folders
        # for Data Processing and Refinement
        self.project_directory = self.settings['initial_model_directory']
        self.reference_directory = self.settings['reference_directory']
        self.Serial=0
        self.Refine=None

        self.xtalID=''
        self.compoundID=''
        self.spider_plot=''
        self.refinement_folder=''

        self.pdb_style='refine.pdb'
        self.mtz_style='refine.mtz'

        # stores imol of currently loaded molecules and maps
        self.mol_dict = {   'protein':  -1,
                            'ligand':   -1,
                            '2fofc':    -1,
                            'fofc':     -1,
                            'event':    -1  }

        # two dictionaries which are flushed when a new crystal is loaded
        # and which contain information to update the data source if necessary
        self.db_dict_mainTable={}
        self.db_dict_panddaTable={}

        ###########################################################################################
        # some COOT settings
        coot.set_map_radius(17)
        coot.set_colour_map_rotation_for_map(0)
#        coot.set_colour_map_rotation_on_read_pdb_flag(21)

        self.QualityIndicators = {  'RefinementRcryst':                         '-',
                                    'RefinementRfree':                          '-',
                                    'RefinementRfreeTraficLight':               'gray',
                                    'RefinementResolution':                     '-',
                                    'RefinementResolutionTL':                   'gray',
                                    'RefinementMolProbityScore':                '-',
                                    'RefinementMolProbityScoreTL':              'gray',
                                    'RefinementRamachandranOutliers':           '-',
                                    'RefinementRamachandranOutliersTL':         'gray',
                                    'RefinementRamachandranFavored':            '-',
                                    'RefinementRamachandranFavoredTL':          'gray',
                                    'RefinementRmsdBonds':                      '-',
                                    'RefinementRmsdBondsTL':                    'gray',
                                    'RefinementRmsdAngles':                     '-',
                                    'RefinementRmsdAnglesTL':                   'gray',
                                    'RefinementMatrixWeight':                   '-'   }

        # default refmac parameters
        self.RefmacParams={ 'HKLIN':            '',                 'HKLOUT': '',
                            'XYZIN':            '',                 'XYZOUT': '',
                            'LIBIN':            '',                 'LIBOUT': '',
                            'TLSIN':            '',                 'TLSOUT': '',
                            'TLSADD':           '',
                            'NCYCLES':          '10',
                            'MATRIX_WEIGHT':    'AUTO',
                            'BREF':             '    bref ISOT\n',
                            'TLS':              '',
                            'NCS':              '',
                            'TWIN':             ''    }



    def StartGUI(self):

        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        self.window.set_default_size(400, 800)
        self.window.set_title("XChemExplorer")
        self.vbox = gtk.VBox()                      # this is the main container

        #################################################################################
        # --- PDB file selection ---
        # checking for pdb files in reference directory
        referenceFiles=[]
        for files in glob.glob(os.path.join(self.reference_directory,'*.pdb')):
            pdbFile=files[files.rfind('/')+1:]
            referenceFiles.append(pdbFile)
        frame = gtk.Frame(label='Select PDB file')
        hbox=gtk.HBox()
        self.cb_select_pdb = gtk.combo_box_new_text()
        self.cb_select_pdb.connect("changed", self.set_selection_mode)
        for pdbFile in referenceFiles:
            self.cb_select_pdb.append_text(pdbFile)
        hbox.add(self.cb_select_pdb)

        self.load_pdb_file_button = gtk.Button(label="Load")
#        self.load_pdb_file_button.connect("clicked",self.get_samples_to_look_at)
        self.load_pdb_file_button.connect("clicked",self.load_pdb_file)
        hbox.add(self.load_pdb_file_button)
        frame.add(self.load_pdb_file_button)
        self.vbox.pack_start(frame)

        # SPACER
        self.vbox.add(gtk.Label(' '))

        frame = gtk.Frame(label='MTZ file to refine against')
        hbox=gtk.HBox()
        self.mtzFree=''
        self.mtzFree_label=gtk.Label()
        hbox.add(self.mtzFree_label)
        self.vbox.pack_start(frame)

        # SPACER
        self.vbox.add(gtk.Label(' '))

        frame = gtk.Frame(label='MTZ file after refinement')
        hbox=gtk.HBox()
        self.mtzRefine_label=gtk.Label()
        hbox.add(self.mtzRefine_label)
        self.vbox.pack_start(frame)



        # SPACER
        self.vbox.add(gtk.Label(' \n \n \n '))


        #################################################################################
        # --- ground state mean map ---
        # checking for ground state mean map in reference folder
        self.ground_state_mean_map=''
        self.ground_state_mean_map_info="map not available; please run 'prepare ground state map'"
        for files in glob.glob(os.path.join(self.reference_directory,'pandda_reference','processed_datasets','*','*ground-state-mean-map.native.ccp4')):
            if os.path.isfile(files):
                self.ground_state_mean_map=files
                xtal=files.split('/')[len(files.split('/'))-2]
                self.ground_state_mean_map_info='ground state mean map: '+xtal
                break

        frame = gtk.Frame(label='Load ground-state-mean-map file')
        hbox=gtk.HBox()

        frameLabel=gtk.Frame()
        self.status_label=gtk.Label()
        frameLabel.add(gtk.Label(self.ground_state_mean_map_info))
        hbox.add(frameLabel)

        self.load_ground_state_map_button = gtk.Button(label="Load")
        self.load_ground_state_map_button.connect("clicked",self.load_ground_state_map)
        hbox.add(self.load_ground_state_map_button)
        frame.add(self.load_ground_state_map_button)
        self.vbox.pack_start(frame)

        if self.ground_state_mean_map=='':
            self.load_ground_state_map_button.set_sensitive(False)

        # SPACER
        self.vbox.add(gtk.Label(' \n \n \n '))





        #################################################################################
        # --- status window ---
        frame=gtk.Frame()
        self.status_label=gtk.Label()
        frame.add(self.status_label)
        self.vbox.pack_start(frame)

        #################################################################################
        # --- Refinement Statistics ---
        # next comes a section which displays some global quality indicators
        # a combination of labels and textview widgets, arranged in a table

        RRfreeLabel_frame=gtk.Frame()
        self.RRfreeLabel = gtk.Label('R/Rfree')
        RRfreeLabel_frame.add(self.RRfreeLabel)
        self.RRfreeValue = gtk.Label(self.QualityIndicators['RefinementRcryst']+'/'+self.QualityIndicators['RefinementRfree'])
        RRfreeBox_frame=gtk.Frame()
        self.RRfreeBox = gtk.EventBox()
        self.RRfreeBox.add(self.RRfreeValue)
        RRfreeBox_frame.add(self.RRfreeBox)

        ResolutionLabel_frame=gtk.Frame()
        self.ResolutionLabel = gtk.Label('Resolution')
        ResolutionLabel_frame.add(self.ResolutionLabel)
        self.ResolutionValue = gtk.Label(self.QualityIndicators['RefinementResolution'])
        ResolutionBox_frame=gtk.Frame()
        self.ResolutionBox = gtk.EventBox()
        self.ResolutionBox.add(self.ResolutionValue)
        ResolutionBox_frame.add(self.ResolutionBox)

        MolprobityScoreLabel_frame=gtk.Frame()
        self.MolprobityScoreLabel = gtk.Label('MolprobityScore')
        MolprobityScoreLabel_frame.add(self.MolprobityScoreLabel)
        self.MolprobityScoreValue = gtk.Label(self.QualityIndicators['RefinementMolProbityScore'])
        MolprobityScoreBox_frame=gtk.Frame()
        self.MolprobityScoreBox = gtk.EventBox()
        self.MolprobityScoreBox.add(self.MolprobityScoreValue)
        MolprobityScoreBox_frame.add(self.MolprobityScoreBox)

        RamachandranOutliersLabel_frame=gtk.Frame()
        self.RamachandranOutliersLabel = gtk.Label('Rama Outliers')
        RamachandranOutliersLabel_frame.add(self.RamachandranOutliersLabel)
        self.RamachandranOutliersValue = gtk.Label(self.QualityIndicators['RefinementRamachandranOutliers'])
        RamachandranOutliersBox_frame=gtk.Frame()
        self.RamachandranOutliersBox = gtk.EventBox()
        self.RamachandranOutliersBox.add(self.RamachandranOutliersValue)
        RamachandranOutliersBox_frame.add(self.RamachandranOutliersBox)

        RamachandranFavoredLabel_frame=gtk.Frame()
        self.RamachandranFavoredLabel = gtk.Label('Rama Favored')
        RamachandranFavoredLabel_frame.add(self.RamachandranFavoredLabel)
        self.RamachandranFavoredValue = gtk.Label(self.QualityIndicators['RefinementRamachandranFavored'])
        RamachandranFavoredBox_frame=gtk.Frame()
        self.RamachandranFavoredBox = gtk.EventBox()
        self.RamachandranFavoredBox.add(self.RamachandranFavoredValue)
        RamachandranFavoredBox_frame.add(self.RamachandranFavoredBox)

        rmsdBondsLabel_frame=gtk.Frame()
        self.rmsdBondsLabel = gtk.Label('rmsd(Bonds)')
        rmsdBondsLabel_frame.add(self.rmsdBondsLabel)
        self.rmsdBondsValue = gtk.Label(self.QualityIndicators['RefinementRmsdBonds'])
        rmsdBondsBox_frame=gtk.Frame()
        self.rmsdBondsBox = gtk.EventBox()
        self.rmsdBondsBox.add(self.rmsdBondsValue)
        rmsdBondsBox_frame.add(self.rmsdBondsBox)

        rmsdAnglesLabel_frame=gtk.Frame()
        self.rmsdAnglesLabel = gtk.Label('rmsd(Angles)')
        rmsdAnglesLabel_frame.add(self.rmsdAnglesLabel)
        self.rmsdAnglesValue = gtk.Label(self.QualityIndicators['RefinementRmsdAngles'])
        rmsdAnglesBox_frame=gtk.Frame()
        self.rmsdAnglesBox = gtk.EventBox()
        self.rmsdAnglesBox.add(self.rmsdAnglesValue)
        rmsdAnglesBox_frame.add(self.rmsdAnglesBox)

        MatrixWeightLabel_frame=gtk.Frame()
        self.MatrixWeightLabel = gtk.Label('Matrix Weight')
        MatrixWeightLabel_frame.add(self.MatrixWeightLabel)
        self.MatrixWeightValue = gtk.Label(self.QualityIndicators['RefinementMatrixWeight'])
        MatrixWeightBox_frame=gtk.Frame()
        self.MatrixWeightBox = gtk.EventBox()
        self.MatrixWeightBox.add(self.MatrixWeightValue)
        MatrixWeightBox_frame.add(self.MatrixWeightBox)

        outer_frame = gtk.Frame()
        hbox = gtk.HBox()

        frame = gtk.Frame()
        self.table_left  = gtk.Table(8, 2, False)
        self.table_left.attach(RRfreeLabel_frame,                 0, 1, 0, 1)
        self.table_left.attach(ResolutionLabel_frame,             0, 1, 1, 2)
        self.table_left.attach(MolprobityScoreLabel_frame,        0, 1, 2, 3)
        self.table_left.attach(RamachandranOutliersLabel_frame,   0, 1, 3, 4)
        self.table_left.attach(RamachandranFavoredLabel_frame,    0, 1, 4, 5)
        self.table_left.attach(rmsdBondsLabel_frame,              0, 1, 5, 6)
        self.table_left.attach(rmsdAnglesLabel_frame,             0, 1, 6, 7)
        self.table_left.attach(MatrixWeightLabel_frame,           0, 1, 7, 8)
        self.table_left.attach(RRfreeBox_frame,                   1, 2, 0, 1)
        self.table_left.attach(ResolutionBox_frame,               1, 2, 1, 2)
        self.table_left.attach(MolprobityScoreBox_frame,          1, 2, 2, 3)
        self.table_left.attach(RamachandranOutliersBox_frame,     1, 2, 3, 4)
        self.table_left.attach(RamachandranFavoredBox_frame,      1, 2, 4, 5)
        self.table_left.attach(rmsdBondsBox_frame,                1, 2, 5, 6)
        self.table_left.attach(rmsdAnglesBox_frame,               1, 2, 6, 7)
        self.table_left.attach(MatrixWeightBox_frame,             1, 2, 7, 8)
        frame.add(self.table_left)
        hbox.add(frame)

        outer_frame.add(hbox)
        self.vbox.add(outer_frame)

        button = gtk.Button(label="Show MolProbity to-do list")
        button.connect("clicked",self.show_molprobity_to_do)
        self.vbox.add(button)
        self.vbox.pack_start(frame)

        # SPACER
        self.vbox.add(gtk.Label(' '))

        # --- refinement & options ---
        self.hbox_for_refinement=gtk.HBox()
        self.REFINEbutton = gtk.Button(label="Refine")
        self.RefinementParamsButton = gtk.Button(label="refinement parameters")
        self.REFINEbutton.connect("clicked",self.REFINE)
        self.hbox_for_refinement.add(self.REFINEbutton)
        self.RefinementParamsButton.connect("clicked",self.RefinementParams)
        self.hbox_for_refinement.add(self.RefinementParamsButton)
        self.vbox.add(self.hbox_for_refinement)

        # --- CANCEL button ---
        self.CANCELbutton = gtk.Button(label="CANCEL")
        self.CANCELbutton.connect("clicked", self.CANCEL)
        self.vbox.add(self.CANCELbutton)

        self.window.add(self.vbox)
        self.window.show_all()

    def CANCEL(self,widget):
        self.window.destroy()


    def RefreshData(self):

        # initialize Refinement library
        self.Refine=XChemRefine.Refine(self.project_directory,self.xtalID,self.compoundID,self.data_source)
        self.Serial=XChemRefine.GetSerial(self.project_directory,self.xtalID)
#        self.Serial=self.Refine.GetSerial()
        if self.Serial==1:
            # i.e. no refinement has been done; data is probably straight out of dimple
            if os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.pdb_style)):
                print '==> XCE: updating quality indicators in data source for '+self.xtalID
                XChemUtils.parse().update_datasource_with_PDBheader(self.xtalID,self.data_source,os.path.join(self.project_directory,self.xtalID,self.pdb_style))
                XChemUtils.parse().update_datasource_with_phenix_validation_summary(self.xtalID,self.data_source,'')   # '' because file does not exist
            elif os.path.isfile(os.path.join(self.project_directory,self.xtalID,'dimple.pdb')):
                print '==> XCE: updating quality indicators in data source for '+self.xtalID
                XChemUtils.parse().update_datasource_with_PDBheader(self.xtalID,self.data_source,os.path.join(self.project_directory,self.xtalID,'dimple.pdb'))
                XChemUtils.parse().update_datasource_with_phenix_validation_summary(self.xtalID,self.data_source,'')   # '' because file does not exist

        # all this information is now updated in the datasource after each refinement cycle
        self.QualityIndicators=self.db.get_db_dict_for_sample(self.xtalID)

        #########################################################################################
        # history
        # if the structure was previously refined, try to read the parameters
#        self.hbox_for_info_graphics.remove(self.canvas)
        if self.Serial > 1:
            self.RefmacParams=self.Refine.ParamsFromPreviousCycle(self.Serial-1)
#            refinement_cycle,Rfree,Rcryst=self.Refine.GetRefinementHistory()
#            self.canvas = FigureCanvas(self.update_plot(refinement_cycle,Rfree,Rcryst))
#        else:
#            self.canvas = FigureCanvas(self.update_plot([0],[0],[0]))  # a gtk.DrawingArea
#        self.canvas.set_size_request(190, 190)
#        self.hbox_for_info_graphics.add(self.canvas)
#        self.canvas.show()

        #########################################################################################
        # update pdb & maps

        #########################################################################################
        # delete old PDB and MAP files
        # - get a list of all molecules which are currently opened in COOT
        # - remove all molecules/ maps before loading a new set
        if len(coot_utils_XChem.molecule_number_list()) > 0:
            for item in coot_utils_XChem.molecule_number_list():
                coot.close_molecule(item)

        #########################################################################################
        # read new PDB files
        # read protein molecule after ligand so that this one is the active molecule
        coot.set_nomenclature_errors_on_read("ignore")
        if os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.compoundID+'.pdb')):
            coot.set_colour_map_rotation_on_read_pdb(0)
            imol=coot.handle_read_draw_molecule_with_recentre(os.path.join(self.project_directory,self.xtalID,self.compoundID+'.pdb'),0)
            self.mol_dict['ligand']=imol
            coot.read_cif_dictionary(os.path.join(self.project_directory,self.xtalID,self.compoundID+'.cif'))
        if not os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.pdb_style)):
            os.chdir(os.path.join(self.project_directory,self.xtalID))

        if self.refinementProtocol=='pandda':
            if os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.pdb_style.replace('.pdb','')+'.split.ground-state.pdb')):
                os.chdir(os.path.join(self.project_directory,self.xtalID))
                coot.set_colour_map_rotation_on_read_pdb(0)
                color_wheel_rotation=160/float(imol+2)
                coot.set_colour_map_rotation_on_read_pdb(color_wheel_rotation)
                imol=coot.handle_read_draw_molecule_with_recentre(os.path.join(self.project_directory,self.xtalID,self.pdb_style.replace('.pdb','')+'.split.ground-state.pdb'),0)
                coot.set_colour_by_molecule(imol)
                coot.set_mol_active(imol,0)
            if os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.pdb_style.replace('.pdb','')+'.split.bound-state.pdb')):
                os.chdir(os.path.join(self.project_directory,self.xtalID))
                coot.set_colour_map_rotation_on_read_pdb(0)
                color_wheel_rotation=21/float(imol+2)
                coot.set_colour_map_rotation_on_read_pdb(color_wheel_rotation)
                imol=coot.handle_read_draw_molecule_with_recentre(os.path.join(self.project_directory,self.xtalID,self.pdb_style.replace('.pdb','')+'.split.bound-state.pdb'),0)
                self.mol_dict['protein']=imol
            else:
                self.go_to_next_xtal()
        else:
            if os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.pdb_style)):
                os.chdir(os.path.join(self.project_directory,self.xtalID))
                imol=coot.handle_read_draw_molecule_with_recentre(os.path.join(self.project_directory,self.xtalID,self.pdb_style),0)
            elif os.path.isfile(os.path.join(self.project_directory,self.xtalID,'dimple.pdb')):
                os.chdir(os.path.join(self.project_directory,self.xtalID))
                imol=coot.handle_read_draw_molecule_with_recentre(os.path.join(self.project_directory,self.xtalID,'dimple.pdb'),0)
            else:
                self.go_to_next_xtal()
            self.mol_dict['protein']=imol

        for item in coot_utils_XChem.molecule_number_list():
            if coot.molecule_name(item).endswith(self.pdb_style.replace('.pdb','')+'.split.bound-state.pdb') or coot.molecule_name(item).endswith(self.pdb_style):
                coot.set_show_symmetry_master(1)    # master switch to show symmetry molecules
                coot.set_show_symmetry_molecule(item,1) # show symm for model

        #########################################################################################
        # read fofo maps
        # - read ccp4 map: 0 - 2fofc map, 1 - fofc.map
        # read 2fofc map last so that one can change its contour level
        if os.path.isfile(os.path.join(self.project_directory,self.xtalID,'2fofc.map')):
            coot.set_colour_map_rotation_on_read_pdb(0)
            coot.set_default_initial_contour_level_for_difference_map(3)
            coot.handle_read_ccp4_map(os.path.join(self.project_directory,self.xtalID,'fofc.map'),1)
            coot.set_default_initial_contour_level_for_map(1)
            coot.handle_read_ccp4_map(os.path.join(self.project_directory,self.xtalID,'2fofc.map'),0)
            coot.set_last_map_colour(0,0,1)
        else:
            # try to open mtz file with same name as pdb file
            coot.set_default_initial_contour_level_for_map(1)
#            if not os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.mtz_style)):
#                os.chdir(os.path.join(self.project_directory,self.xtalID))
#                if not os.path.isfile('REFINEMENT_IN_PROGRESS'):
#                    if os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.xtalID+'-pandda-input.mtz')):
#                        os.symlink(self.xtalID+'-pandda-input.mtz',self.mtz_style)
#                    elif os.path.isfile(os.path.join(self.project_directory,self.xtalID,'dimple.mtz')):
#                        os.symlink('dimple.mtz',self.mtz_style)
            if os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.mtz_style)):
                coot.auto_read_make_and_draw_maps(os.path.join(self.project_directory,self.xtalID,self.mtz_style))
            elif os.path.isfile(os.path.join(self.project_directory,self.xtalID,'dimple.mtz')):
                coot.auto_read_make_and_draw_maps(os.path.join(self.project_directory,self.xtalID,'dimple.mtz'))

        #########################################################################################
        # update Quality Indicator table
        try:
            self.RRfreeValue.set_label(  str(round(float(self.QualityIndicators['RefinementRcryst']),3)) +' / '+str(round(float(self.QualityIndicators['RefinementRfree']),3)))
        except ValueError:
            self.RRfreeValue.set_label('-')

        try:
            self.RRfreeBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RefinementRfreeTraficLight']))
        except ValueError:
            pass
        self.ResolutionValue.set_label(self.QualityIndicators['RefinementResolution'])
        try:
            self.ResolutionBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RefinementResolutionTL']))
        except ValueError:
            pass
        self.MolprobityScoreValue.set_label(self.QualityIndicators['RefinementMolProbityScore'])
        try:
            self.MolprobityScoreBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RefinementMolProbityScoreTL']))
        except ValueError:
            pass
        self.RamachandranOutliersValue.set_label(self.QualityIndicators['RefinementRamachandranOutliers'])
        try:
            self.RamachandranOutliersBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RefinementRamachandranOutliersTL']))
        except ValueError:
            pass
        self.RamachandranFavoredValue.set_label(self.QualityIndicators['RefinementRamachandranFavored'])
        try:
            self.RamachandranFavoredBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RefinementRamachandranFavoredTL']))
        except ValueError:
            pass
        self.rmsdBondsValue.set_label(self.QualityIndicators['RefinementRmsdBonds'])
        try:
            self.rmsdBondsBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RefinementRmsdBondsTL']))
        except ValueError:
            pass
        self.rmsdAnglesValue.set_label(self.QualityIndicators['RefinementRmsdAngles'])
        try:
            self.rmsdAnglesBox.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(self.QualityIndicators['RefinementRmsdAnglesTL']))
        except ValueError:
            pass
        self.MatrixWeightValue.set_label(self.QualityIndicators['RefinementMatrixWeight'])

        try:
            pic = gtk.gdk.pixbuf_new_from_file(os.path.join(self.project_directory,self.xtalID,self.compoundID+'.png'))
        except gobject.GError:
            pic = gtk.gdk.pixbuf_new_from_file(os.path.join(os.getenv('XChemExplorer_DIR'),'image','NO_COMPOUND_IMAGE_AVAILABLE.png'))
        self.pic = pic.scale_simple(190, 190, gtk.gdk.INTERP_BILINEAR)
        self.image.set_from_pixbuf(self.pic)

    def go_to_next_xtal(self):
        self.index+=1
        if self.index >= len(self.Todo):
            self.index = len(self.Todo)
        self.cb.set_active(self.index)


    def REFINE(self,widget):

        #######################################################
        if not os.path.isdir(os.path.join(self.project_directory,self.xtalID,'cootOut')):
            os.mkdir(os.path.join(self.project_directory,self.xtalID,'cootOut'))
        # create folder for new refinement cycle
        os.mkdir(os.path.join(self.project_directory,self.xtalID,'cootOut','Refine_'+str(self.Serial)))

        #######################################################
        # write PDB file
        # now take protein pdb file and write it to newly create Refine_<serial> folder
        # note: the user has to make sure that the ligand file was merged into main file
        for item in coot_utils_XChem.molecule_number_list():
            if coot.molecule_name(item).endswith(self.pdb_style.replace('.pdb','')+'.split.bound-state.pdb') or coot.molecule_name(item).endswith(self.pdb_style):
                coot.write_pdb_file(item,os.path.join(self.project_directory,self.xtalID,'cootOut','Refine_'+str(self.Serial),'refine.modified.pdb'))
                break
            elif coot.molecule_name(item).endswith('dimple.pdb'):
                coot.write_pdb_file(item,os.path.join(self.project_directory,self.xtalID,'cootOut','Refine_'+str(self.Serial),'refine.modified.pdb'))
                break


        #######################################################
        if self.refinementProtocol=='pandda':
            XChemRefine.panddaRefine(self.project_directory,self.xtalID,self.compoundID,self.data_source).RunQuickRefine(self.Serial,self.RefmacParams,self.external_software,self.xce_logfile)
        else:
            print 'sorry, not programmed yet'

#            self.Refine.RunRefmac(self.Serial,self.RefmacParams,self.external_software,self.xce_logfile)

        self.index+=1
        if self.index >= len(self.Todo):
#            self.index = len(self.Todo)
            self.index = 0
        self.cb.set_active(self.index)


    def RefinementParams(self,widget):
        print '\n==> XCE: changing refinement parameters'
        self.RefmacParams=XChemRefine.RefineParams(self.project_directory,self.xtalID,self.compoundID,self.data_source).RefmacRefinementParams(self.RefmacParams)

    def set_selection_mode(self,widget):
        self.selection_mode=widget.get_active_text()


    def load_pdb_file(self,widget):
        # first remove all pdb files
        if len(coot_utils_XChem.molecule_number_list()) > 0:
            for item in coot_utils_XChem.molecule_number_list():
                if coot.molecule_name(item).endswith('.pdb'):
                    coot.close_molecule(item)

        pdbFile=self.cb_select_pdb.get_active_text()

        if os.path.isfile(os.path.join(self.project_directory,pdbFile)):
            os.chdir(os.path.join(self.project_directory,self.xtalID))
            coot.set_colour_map_rotation_on_read_pdb(0)
            color_wheel_rotation=160/float(imol+2)
            coot.set_colour_map_rotation_on_read_pdb(color_wheel_rotation)
            imol=coot.handle_read_draw_molecule_with_recentre(os.path.join(self.project_directory,pdbFile),0)
            coot.set_colour_by_molecule(imol)
            coot.set_mol_active(imol,0)

        mtzFree=pdbFile.replace('.pdb','')+'.free.mtz'
        if os.path.isfile(os.path.join(self.project_directory,mtzFree)):
            self.mtzFree=os.path.join(self.project_directory,mtzFree)
            self.mtzFree_label.set_text(mtzFree)
        else:
            self.mtzFree_label.set_text('cannot find %s in reference directory' %mtzFree)

        mtzRefine=pdbFile.replace('.pdb','')+'.refine.mtz'
        if os.path.isfile(os.path.join(self.project_directory,mtzRefine)):
            self.mtzRefine=os.path.join(self.project_directory,mtzRefine)
            self.mtzRefine_label.set_text(mtzFree)
        else:
            self.mtzRefine_label.set_text('cannot find %s in reference directory' %mtzRefine)

    def load_ground_state_map(self,widget):
        # first remove all ground state maps files
        if len(coot_utils_XChem.molecule_number_list()) > 0:
            for item in coot_utils_XChem.molecule_number_list():
                if 'ground-state-mean-map' in coot.molecule_name(item):
                    coot.close_molecule(item)

        coot.set_colour_map_rotation_on_read_pdb(0)
        coot.handle_read_ccp4_map((self.ground_state_mean_map),0)
        for imol in coot_utils_XChem.molecule_number_list():
            if self.ground_state_mean_map in coot.molecule_name(imol):
                coot.set_contour_level_in_sigma(imol,2)
                coot.set_last_map_colour(0.74,0.44,0.02)

    def show_molprobity_to_do(self,widget):
        if os.path.isfile(os.path.join(self.project_directory,self.xtalID,'Refine_'+str(self.Serial-1),'molprobity_coot.py')):
            print '==> XCE: running MolProbity Summary for',self.xtalID
            coot.run_script(os.path.join(self.project_directory,self.xtalID,'Refine_'+str(self.Serial-1),'molprobity_coot.py'))
        else:
            print '==> XCE: cannot find '+os.path.join(self.project_directory,self.xtalID,'Refine_'+str(self.Serial-1),'molprobity_coot.py')

    def refinementProtocolCallback(self, widget):
        if widget.get_active():
            self.refinementProtocol='pandda'
            self.PREVbuttonSite.set_sensitive(True)
            self.NEXTbuttonSite.set_sensitive(True)
        else:
            self.refinementProtocol='refmac'
            self.PREVbuttonSite.set_sensitive(False)
            self.NEXTbuttonSite.set_sensitive(False)





if __name__=='__main__':
    GUI().StartGUI()