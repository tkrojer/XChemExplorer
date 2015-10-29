import sys
import os
import glob
sys.path.append('/usr/lib64/python2.6/site-packages/gtk-2.0')
import gtk
#import pygtk
#import pango
#sys.path.append(os.getenv('TINSEL_DIR')+'/lib')
#import GetAimless
from datetime import datetime
import threading
import gobject 
gtk.gdk.threads_init()
import time
import math
#sys.path.append('/dls/labxchem/data/2015/lb13385-1/processing/_PROCESSING_TEST/BLing/lib')
sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
from XChemUtils import process
from XChemUtils import parse
from XChemUtils import queue
from XChemUtils import mtztools

class TInsel:
    def __init__(self):
        # Directories
        self.CurrentDir=os.getcwd()
        counter=0
        for n,char in enumerate(self.CurrentDir):
            if char=='/': counter+=1
            if counter==6:
                self.ProjectDir=self.CurrentDir[:n]
                break
        self.BeamlineDir=self.ProjectDir+'/processing/beamline'
        self.InitialModelDir=self.ProjectDir+'/processing/analysis/initial_model'
        self.RefineModelDir=self.ProjectDir+'/processing/analysis/refine_model'
        self.FindHitsDir=self.ProjectDir+'/processing/analysis/find_hits'
        self.DatabaseDir=self.ProjectDir+'/processing/database'
        self.reference_directory=self.ProjectDir+'/processing/reference'
        self.reference_file_root='reference'

        # Targets
        self.VisitList=[]
        self.TargetList=[]
        self.Target=''
        for dir in glob.glob(self.BeamlineDir+'/*'):
            self.VisitList.append(os.path.realpath(dir))
            for target in glob.glob(os.path.realpath(dir)+'/processed/*'):
                if target[target.rfind('/')+1:] not in ['results','README-log','edna-latest.html']:
                    if target[target.rfind('/')+1:] not in self.TargetList: self.TargetList.append(target[target.rfind('/')+1:])

        # some lists
        self.LOGfiles = []
        self.Images = []
        self.Runs = []
        self.ExperimentList=[]

        # default settings
        self.XChemDatabase=''


#        self.VisitList=SelectTarget.GetTarget().GetVisitIDs()
        self.Visit=''
        self.VisitDir=''
        self.BestResoFile=[]


        self.DatasetList=[]
        self.DatasetOutcomeList=[]
        self.TreeviewList=[]
        self.inital_model_treeview=[]

        # flag which gets set to '1' when TInsel does something
        self.Active=0
        self.load_initial_models_first_time=0
        # text for statusbar
        self.status_bar_message=''
        self.reference_file_root=''


        # main window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        # it seems(!) that there is a frequent X-windows error in case
        # a certain window size is requested
        # self.window.set_default_size(1920, 1200)
        self.window.set_title("XChemExplorer")
        vboxMain = gtk.VBox()

        # Target selection
        hboxTarget=gtk.HBox()
        hboxTarget.add(gtk.Label('Target'))
        self.cbTarget = gtk.combo_box_new_text()
        self.cbTarget.connect("changed", self.ChooseTarget)
        for item in self.TargetList:
                self.cbTarget.append_text(item)
        hboxTarget.add(self.cbTarget)
        vboxMain.pack_start(hboxTarget,False)

        # Data Collection
        data_collection_vbox=gtk.VBox(False,3)
        self.data_collection_scrolled=gtk.ScrolledWindow()
        self.data_collection_scrolled.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
        data_collection_vbox.add(self.data_collection_scrolled)
        data_collection_button_hbox=gtk.HBox(True,5)
        halign_data_collection_button_hbox= gtk.Alignment(0, 0, 0, 0)
        get_data_collection_button=gtk.Button(label="\nRead Data\n")
        get_data_collection_button.connect("clicked",self.button_task,"read_data")
        data_collection_button_hbox.add(get_data_collection_button)
        write_files_button = gtk.Button(label="\nWrite Files\n")
        write_files_button.connect("clicked", self.button_task,"write_files")
        data_collection_button_hbox.add(write_files_button)
        halign_data_collection_button_hbox.add(data_collection_button_hbox)
        data_collection_vbox.pack_start(halign_data_collection_button_hbox,False,False,3)

        # Initial Model
        initial_model_vbox=gtk.VBox(False,3)

        toggle_dimple_run_button = gtk.CheckButton('(de)-select all sample for dimple')
        toggle_dimple_run_button.connect("toggled", self.initial_model_toggle_dimple_run)
        initial_model_vbox.pack_start(toggle_dimple_run_button,False,False,3)

        set_reference_hbox=gtk.HBox(True,3)
        halign_initial_model_button_hbox= gtk.Alignment(0, 0, 0, 0)
        initial_model_choose_reference=gtk.combo_box_new_text()
        initial_model_choose_reference.connect("changed", self.initial_model_choose_reference)
        for files in glob.glob(self.reference_directory+'/*'):
            if files.endswith('.pdb'):
                reference_root=files[files.rfind('/')+1:files.rfind('.')]
                initial_model_choose_reference.append_text(reference_root)
        set_reference_hbox.add(initial_model_choose_reference)
        set_reference_button=gtk.Button(label="\nset new reference if applicable\n")
        set_reference_button.connect("clicked",self.button_task,"initial_refinement_set_reference")
        set_reference_hbox.add(set_reference_button)
        halign_initial_model_button_hbox.add(set_reference_hbox)
        initial_model_vbox.pack_start(halign_initial_model_button_hbox,False,False,3)


        self.initial_model_scrolled=gtk.ScrolledWindow()
        self.initial_model_scrolled.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
        initial_model_vbox.add(self.initial_model_scrolled)
        initial_model_button_hbox=gtk.HBox(True,5)
        halign_initial_model_button_hbox= gtk.Alignment(0, 0, 0, 0)
        get_samples_initial_model_button=gtk.Button(label="\nLoad Samples\n")
        get_samples_initial_model_button.connect("clicked",self.button_task,"load_samples")
        initial_model_button_hbox.add(get_samples_initial_model_button)
        run_dimple_button=gtk.Button(label="\nRun Dimple\n")
        run_dimple_button.connect("clicked",self.button_task,"run_dimple")
        initial_model_button_hbox.add(run_dimple_button)
        refresh_initial_model_button=gtk.Button(label="\nRefresh\n")
        refresh_initial_model_button.connect("clicked",self.button_task,"refresh_initial_refinement")
        initial_model_button_hbox.add(refresh_initial_model_button)
        halign_initial_model_button_hbox.add(initial_model_button_hbox)
        initial_model_vbox.pack_start(halign_initial_model_button_hbox,False,False,3)

        # PANDDAs
        PANDDAsTemp=gtk.Label('under construction')

        # Summary
        SummaryTemp=gtk.Label('under construction')

        # Queue Control
        queue_control_vbox=gtk.VBox(False,3)
        self.queue_control_scrolled=gtk.ScrolledWindow()
        self.queue_control_scrolled.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
        queue_control_vbox.add(self.queue_control_scrolled)
        queue_control_button_hbox=gtk.HBox(True,5)
        halign_queue_control_button_hbox= gtk.Alignment(0, 0, 0, 0)
        get_jobs_in_queue_button=gtk.Button(label="\nGet Jobs\n")
        get_jobs_in_queue_button.connect("clicked",self.button_task,"get_queue_jobs")
        queue_control_button_hbox.add(get_jobs_in_queue_button)
        remove_queue_jobs_button=gtk.Button(label="\nRemove Jobs\n")
        remove_queue_jobs_button.connect("clicked",self.button_task,"remove_queue_jobs")
        queue_control_button_hbox.add(remove_queue_jobs_button)
        refresh_queue_button=gtk.Button(label="\nRefresh\n")
        refresh_queue_button.connect("clicked",self.button_task,"refresh_queue_jobs")
        queue_control_button_hbox.add(refresh_queue_button)

        halign_queue_control_button_hbox.add(queue_control_button_hbox)
        queue_control_vbox.pack_start(halign_queue_control_button_hbox,False,False,3)


        # Settings
        SettingsTemp=gtk.Label('under construction')

        # Create a new notebook, place the position of the tabs
        self.notebook = gtk.Notebook()
        self.notebook.set_tab_pos(gtk.POS_TOP)
        self.notebook.show()
        self.show_tabs = True
        self.show_border = True

        Tab=['Data Collection','Initial Model','PANDDAs','Summary','Queue Control','Settings']
        NotebookPage = [    data_collection_vbox,
                            initial_model_vbox,
                            PANDDAsTemp,
                            SummaryTemp,
                            queue_control_vbox,
                            SettingsTemp]
        for i in range(6):
            label = gtk.Label(Tab[i])
            self.notebook.append_page(NotebookPage[i],label)
        vboxMain.add(self.notebook)

        hbox = gtk.HBox()

        CANCELbutton = gtk.Button(label="CANCEL")
        CANCELbutton.connect("clicked",self.button_task,"cancel")
        hbox.add(CANCELbutton)

        vboxMain.pack_start(hbox, False)

        # Statusbar
        hboxStatus = gtk.HBox()
        self.status_bar = gtk.Statusbar()
        self.context_id=self.status_bar.get_context_id('progress update')
        hboxStatus.add(self.status_bar)
        # Progressbar
        self.progress_bar = gtk.ProgressBar()
        hboxStatus.add(self.progress_bar)
        vboxMain.pack_start(hboxStatus, False)

        # False makes that buttons have the minimal size!

        self.window.add(vboxMain)
        self.window.show_all()


    def WRITEFILES(self):       
        csv_header =  'SampleID,Visit,Date,Outcome,'
        csv_header += 'Program,Run,SpaceGroup,UnitCell,Resolution Overall,Resolution Low,Resolution High,'
        csv_header += 'RmergeOverall,RmergeLow,RmergeHigh,IsigOverall,IsigLow,IsigHigh,'
        csv_header += 'CompletenessOverall,CompletenessLow,CompletenessHigh,MultiplicityOverall,MultiplicityLow,MultiplicityHigh,'
        csv_header += 'Path to Logfile\n' 

        datasetOutcome=''
        outCSV=csv_header

        fraction=0
        fractionIncrement=1/float(len(self.DatasetOutcomeList))
        self.status_bar_message='writing files'
        gobject.idle_add(self.set_status_bar_message, self.context_id)

        for n,button in enumerate(sorted(self.DatasetOutcomeList)):
            active_radios = [r for r in button[1].get_group() if r.get_active()][0]
            datasetOutcome = active_radios.get_label()
#            outCSV=outCSV+button[0]+','+datasetOutcome
            Visit='n/a'
            Timestamp='n/a'
            Run='n/a'           
            if datasetOutcome.startswith('success'):
                for row in self.TreeviewList:
                    # just make sure again that there really is a table
                    # in case the user forgot to set it not to failed
                    if len(row)==2 and row[0]==button[0]:
#                        print len(row),row,button[0]
                        tree_selection = row[1].get_selection()
                        (model, pathlist) = tree_selection.get_selected_rows()
                        for path in pathlist:
#                            print path[0]
                            for dataset in self.DatasetList:
                                if dataset[0]==button[0]:
#                                    print dataset[path[0]+1]
                                    if dataset[path[0]+1] != None:
                                        Logfile=dataset[path[0]+1]
                                        for Experiment in self.ExperimentList:
                                            if Experiment[0]==button[0]:
                                                for runs in Experiment[2]:
                                                    if runs[0][2] in dataset[path[0]+1]:
                                                        Visit=runs[0][2]
                                                        Timestamp=runs[0][1]
                                                        Run=runs[0][0]
#                                        outCSV=outCSV+Visit+','+Timestamp+','+datasetOutcome+','
                                        Aimless=parse().GetAimlessLog(dataset[path[0]+1])
                                        csvIN =  button[0]+','+Visit+','+Timestamp+','+datasetOutcome+','
                                        csvIN += Aimless['AutoProc']+','+Run+','+Aimless['SpaceGroup']+','+Aimless['UnitCell']+','
                                        csvIN += Aimless['ResolutionLow']+'-'+Aimless['ResolutionHigh']+','+Aimless['ResolutionLow']+'-'+Aimless['ResolutionLowInnerShell']+','
                                        csvIN += Aimless['ResolutionHighOuterShell']+'-'+Aimless['ResolutionHigh']+','
                                        csvIN += Aimless['RmergeOverall']+','+Aimless['RmergeLow']+','+Aimless['RmergeHigh']+','
                                        csvIN += Aimless['IsigOverall']+','+Aimless['IsigLow']+','+Aimless['IsigHigh']+','
                                        csvIN += Aimless['CompletenessOverall']+','+Aimless['CompletenessLow']+','+Aimless['CompletenessHigh']+','
                                        csvIN += Aimless['MultiplicityOverall']+','+Aimless['MultiplicityLow']+','+Aimless['MultiplicityHigh']+','
                                        csvIN += dataset[path[0]+1]+'\n'
#                                        print csvIN
#                                        print outCSV
                                        outCSV += csvIN
                                        if not os.path.isdir(self.InitialModelDir+'/'+button[0]+'/autoprocessing'):
                                            os.mkdir(self.InitialModelDir+'/'+button[0])
                                            os.mkdir(self.InitialModelDir+'/'+button[0]+'/autoprocessing')
                                            os.chdir(self.InitialModelDir+'/'+button[0]+'/autoprocessing')
                                            self.status_bar_message='copying '+str(button[0])
                                            gobject.idle_add(self.set_status_bar_message, self.context_id)
                                            if 'fast_dp' in Logfile:
                                                os.system('cp -R %s .' %Logfile[:Logfile.find('fast_dp')+7])
                                                os.chdir(self.InitialModelDir+'/'+button[0])
                                                os.system("ctruncate -hklin autoprocessing/fast_dp/fast_dp.mtz -hklout autoprocessing/fast_dp/ctruncate.mtz -colin '/*/*/[IMEAN,SIGIMEAN]' > autoprocessing/fast_dp/ctruncate.log")
                                                os.system('ln -s autoprocessing/fast_dp/aimless.log %s.log' %button[0])
                                                os.system('ln -s autoprocessing/fast_dp/ctruncate.mtz %s.mtz' %button[0])
                                            else:
                                                os.system('cp -R %s .' %Logfile[:Logfile.find('LogFiles')-1])
                                                os.chdir(self.InitialModelDir+'/'+button[0])
                                                for datafile in glob.glob('autoprocessing/*/DataFiles/*'):
                                                    if datafile.endswith('free.mtz'):
                                                        os.system('ln -s %s %s.mtz' %(datafile,button[0]))                                                            
                                                        break
                                                for logfile in glob.glob('autoprocessing/*/LogFiles/*'):
                                                    if logfile.endswith('aimless.log'):
                                                        os.system('ln -s %s %s.log' %(logfile,button[0]))
                                                        break



#            os.mkdir(self.CurrentDir+'/../analysis/initial_model/'+item[0])
#            os.chdir(self.CurrentDir+'/../analysis/initial_model/'+item[0])



            else:
                for Experiment in self.ExperimentList:
                    if Experiment[0]==button[0]:
                        for runs in Experiment[2]:
                            Visit=runs[0][2]
                            Timestamp=runs[0][1]
                            Run=runs[0][0]
                            csvIN =  button[0]+','+Visit+','+Timestamp+','+datasetOutcome+'\n'
                            outCSV += csvIN
                            break
                            
            fraction += fractionIncrement
            gobject.idle_add(self.set_progress_bar_fraction, fraction)


        print outCSV
        f=open(os.path.join(self.DatabaseDir,'TInsel.csv'),'w')
        f.write(outCSV)
        f.close()

        self.status_bar_message='idle'
        gobject.idle_add(self.set_status_bar_message, self.context_id)
        gobject.idle_add(self.set_progress_bar_fraction, 0.0)
        self.Active=0





    def ChooseTarget(self,widget):
        self.Target=widget.get_active_text()


    def button_task(self,widget,data=None):
        if self.Target != '' and data=='read_data' and self.Active==0:
            self.Active=1
            threading.Thread(target=self.AutoprocessingResults, args=()).start()
        elif self.Target != '' and data=='load_samples' and self.Active==0 \
             and self.load_initial_models_first_time==0:
            self.Active=1
            threading.Thread(target=self.InitialRefinement, args=()).start()
            self.load_initial_models_first_time=1
        elif self.Target != '' and data=='write_files' and self.Active==0:
            self.Active=1
            threading.Thread(target=self.WRITEFILES, args=()).start()
        elif self.Target != '' and data=='run_dimple' and self.Active==0:
            self.Active=1
            threading.Thread(target=self.RunDimple, args=()).start()
        elif self.Target != '' and data=='refresh_inital_refinement' and self.Active==0:
            self.initial_model_liststore.clear()
            self.initial_model_scrolled.remove(self.initial_model_treeview)
            self.Active=1
            threading.Thread(target=self.InitialRefinement, args=()).start()
        elif data=='get_queue_jobs' and self.Active==0:
            self.Active=1
            threading.Thread(target=self.show_queue_jobs, args=()).start()    
        elif data=='refresh_queue_jobs' and self.Active==0:
            self.Active=1
            self.queue_jobs_liststore.clear()
            self.queue_control_scrolled.remove(self.queue_jobs_treeview)
            threading.Thread(target=self.show_queue_jobs, args=()).start()    
        elif data=='remove_queue_jobs' and self.Active==0:
            self.Active=1
            threading.Thread(target=self.remove_queue_jobs, args=()).start()            
        elif data=='cancel':
            self.window.destroy()
            quit()
        elif data=='initial_refinement_set_reference' and self.Active==0:
            self.Active=1
            try:
                if self.reference_file_root != '':
                    self.initial_model_liststore.clear()
                    self.initial_model_scrolled.remove(self.initial_model_treeview)
                    self.Active=1
                    threading.Thread(target=self.InitialRefinement, args=()).start()
            except AttributeError:
                print 'no samples loaded'
        else:
            buff = '-> please select a target first'
            self.status_bar.push(self.context_id, buff)


    def set_status_bar_message(self, data):
        self.status_bar.push(data, self.status_bar_message)

    def set_progress_bar_fraction(self, fraction):
        if fraction >= 1:
            fraction=1
        self.progress_bar.set_fraction(fraction)

    def samples_in_directory(self,folder):
        samples=0
        for sample in glob.glob(folder+'/*'):
            samples+=1
        return samples


    def AutoprocessingResults(self):
        # assumption: if a sample in TInsel.csv exists and dataset_outcome=success
        #             then we will ignore it! Even if another maybe better run exists somewhere
        samples_to_ignore=[]
        if os.path.isfile(os.path.join(self.DatabaseDir,'TInsel.csv')):
            for line in open(os.path.join(self.DatabaseDir,'TInsel.csv')):
                if line.split(',')[3].startswith('success'):
                    samples_to_ignore.append(line.split(',')[0])

        self.ExperimentList=[]
        fraction=0
        number_of_visits=len(self.VisitList)
        number_of_current_visit=1
        for VisitDir in sorted(self.VisitList):
            Visit=VisitDir.split('/')[5]
            self.status_bar_message=(   'step 1 of 2: searching for logfiles in '
                                        'visit '+str(number_of_current_visit)+
                                        ' of '+str(number_of_visits)+' ('+Visit+')' )
            gobject.idle_add(self.set_status_bar_message, self.context_id)
            NumberXtalsInVisit=self.samples_in_directory(VisitDir+'/processed/'+self.Target)
            for xtals in sorted(glob.glob(VisitDir+'/processed/'+self.Target+'/*')):
                fractionIncrement=1/float(NumberXtalsInVisit)
                Runs=[]
                Logfiles=[]
                Images=[]
                Sample=xtals[xtals.rfind('/')+1:]
                if Sample in samples_to_ignore:
                    continue
                fraction += fractionIncrement
                gobject.idle_add(self.set_progress_bar_fraction, fraction)
                for runs in glob.glob(xtals+'/*'):
                    run=runs[runs.rfind('/')+1:]
                    timestamp=datetime.fromtimestamp(os.path.getmtime(runs)).strftime('%Y-%m-%d %H:%M:%S')
                    Runs.append([(run,timestamp,Visit)])
                    for (path, dirs, files) in os.walk(runs):
                        if 'edna' in dirs:
                            dirs.remove('edna')
                        if 'auto_mrbump' in dirs:
                            dirs.remove('auto_mrbump')
                        if 'fast_ep' in dirs:
                            dirs.remove('fast_ep')
                        if 'multi-xia2' in dirs:
                            dirs.remove('multi-xia2')
                        for item in files:
                            if item.endswith('aimless.log') and self.Target in path:
                                Logfiles.append(path+'/'+item)
                                continue
                for image in glob.glob(VisitDir+'/jpegs/'+self.Target+'/'+Sample+'/*'):
                    if image.endswith('t.png'):
                        Images.append(image)
                    if image.endswith('thumb.jpeg'):
                        Images.append(image)
                    if image.endswith('_.png'):
                        Images.append(image)
                # need to check if crystal was measured during another visit
                # doing this the really ugly way though...
                temp=[]
                found=0
                for Experiment in self.ExperimentList:
                    if Experiment[0]==Sample and not Experiment[1]==Visit:
#                        if Sample=='ATAD2A-x416': print Images
                        Experiment[1].extend([Visit])
                        Experiment[2].extend(Runs)
                        Experiment[3].extend(Logfiles)
                        Experiment[4].extend(sorted(Images))
                        temp.append([Sample,Experiment[1],Experiment[2],Experiment[3],Experiment[4]])
                        found=1
                    else:
                        temp.append(Experiment)
                if not found:
                    temp.append([Sample,[Visit],Runs,Logfiles,sorted(Images)])
                self.ExperimentList=temp
            fraction=0
            number_of_current_visit+=1

        DatasetOutcome = [  "success",
                            "Failed - centring failed",
                            "Failed - no diffraction",
                            "Failed - processing barfs",
                            "Failed - loop empty",
                            "Failed - low resolution",
                            "Failed - no X-rays",
                            "Failed - unknown" ]

        DiffractionDataColumnName = [   'Visit',
                                        'Program',
                                        'Run',
                                        'SpaceGroup',
                                        'Unit Cell',
                                        'Resolution\nOverall',
                                        'Resolution\nInner Shell',
                                        'Resolution\nOuter Shell',
                                        'Rmerge\nOverall',
                                        'Rmerge\nInner Shell',
                                        'Rmerge\nOuter Shell',
                                        'Mn(I/sig(I))\nOverall',
                                        'Mn(I/sig(I))\nInner Shell',
                                        'Mn(I/sig(I))\nOuter Shell',
                                        'Completeness\nOverall',
                                        'Completeness\nInner Shell',
                                        'Completeness\nOuter Shell',
                                        'Multiplicity\nOverall',
                                        'Multiplicity\nInner Shell',
                                        'Multiplicity\nOuter Shell' ]

#        out=''
#        for Experiment in ExperimentList:
#            out=out+str(Experiment)+'\n'
#        f=open('test.txt','w')
#        f.write(out)
#        f.close()
#        print self.ExperimentList

        fraction=0
        fractionIncrement=1/float(len(self.ExperimentList))
        self.status_bar_message='step 2 of 2: building widgets'
        gobject.idle_add(self.set_status_bar_message, self.context_id)
        table = gtk.Table(len(self.ExperimentList), 2, False)
        for n, Experiment in enumerate(self.ExperimentList):
            Sample=Experiment[0]
#            if Sample=='ATAD2A-x553': print Experiment
            # a nested list that in the end will look like this:
            # [ [[xtal1],[ButtonObject1,ButtonObject2,...]],[[xtal2],[ButtonObject1,ButtonObject2,...]],...]
            DatasetListTemp=[]
            DatasetListTemp.append(Sample)
            DatasetOutcomeListTemp=[]
            DatasetOutcomeListTemp.append(Sample)
            TreeviewListTemp=[]
            TreeviewListTemp.append(Sample)

#            hboxSample=gtk.HBox()
#            hboxSampleID=gtk.HBox()
#            hboxSampleID.add(gtk.Label(Sample))
#            hboxSample.add(hboxSampleID)
            table.attach(gtk.Label(Sample), 0, 1, n, n+1)

            vboxOverview=gtk.VBox()
            Space=gtk.Label('')
            vboxOverview.add(Space)

            # Table with autoprocessing results
            # first sort list of logfiles so that the one with highest resolution comes first
            # doing it already here so if files do not meet criteria or if there is no logfile
            # then datasetoutcome='Failed unknown'
            temp=[]
            LogFilesOnly=[]
            for logfile in Experiment[3]:
                if not logfile==None and not logfile==Sample:
                    LogFilesOnly.append(logfile)
                    for line in open(logfile):
                        if line.startswith('High resolution limit'): temp.append(line.split()[3])
            SortedFiles=[]
            if temp != []:
                Min = temp.index(min(temp))
                for u,item in enumerate(LogFilesOnly):
                    if u == Min:
                        SortedFiles.append(item)
                        DatasetListTemp.append(item)
                        self.BestResoFile.append([Sample[0],item])
                for item in LogFilesOnly:
                    if item not in SortedFiles:
                        SortedFiles.append(item)
                        DatasetListTemp.append(item)
            else:
                SortedFiles=[None]
                DatasetListTemp.append(None)
            self.DatasetList.append(DatasetListTemp)



            # top row: radio buttons for XDS outcome
            hboxOutcome = gtk.HBox()
            for index,outcome in enumerate(DatasetOutcome):
                if index==0:
                    OutComeReferenceButton=gtk.RadioButton(None,outcome)
                    DatasetOutcomeListTemp.append(OutComeReferenceButton)
#                    hboxOutcome.add(OutComeReferenceButton)
                    hboxOutcome.pack_start(OutComeReferenceButton,False,False,12)
                else:
                    OutcomeButton=gtk.RadioButton(OutComeReferenceButton,outcome)
                    DatasetOutcomeListTemp.append(OutcomeButton)
#                    hboxOutcome.add(OutcomeButton)
                    hboxOutcome.pack_start(OutcomeButton,False,False,12)
            vboxOverview.add(hboxOutcome)
            Space=gtk.Label('')
            vboxOverview.add(Space)

            if SortedFiles == [None]:
                for q,button in enumerate(DatasetOutcomeListTemp):
                    if q==8:
                        button.set_active(True)
                    

            # Picture container
            for run in Experiment[2]:
                hboxImage=gtk.HBox()
                image_header=gtk.Label('   ----- '+run[0][0]+' ('+run[0][2]+', '+run[0][1]+') -----')
                image_header.set_alignment(xalign=0,yalign=0.5)
                vboxOverview.add(image_header)
                for image in Experiment[4]:
                    if run[0][0] in image and run[0][2] in image:
                        imageGTK = gtk.Image()
                        if os.path.getsize(image) != 0:
                            pic = gtk.gdk.pixbuf_new_from_file(image)
                            scaled_pic=pic.scale_simple(256,192,gtk.gdk.INTERP_BILINEAR)
                            imageGTK.set_from_pixbuf(scaled_pic)
                            hboxImage.pack_start(imageGTK,False,False,2)
                vboxOverview.pack_start(hboxImage)

            # Experiment[0]     = xtalID
            # Experiment[1]     = visitList
            # Experiment[2]     = runs      ; [0] run, [1] timestamp, [2] visit
            # Experiment[3]     = logfiles
            # Experiment[4]     = images


            liststore = gtk.ListStore(str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str)
            if SortedFiles != [None]:
                for files in SortedFiles:
                    for runs in Experiment[2]:
                        visit='?'
                        timestamp='n/a'
                        if runs[0][2] in files:
                            visit=runs[0][2]
                            timestamp=runs[0][1]
                            break
                        else:
                            continue
                    if files != None:
                        Aimless=parse().GetAimlessLog(files)
                        rowIN=(    visit,
                                   Aimless['AutoProc'],
                                   Aimless['Run'],
                                   Aimless['SpaceGroup'],
                                   Aimless['UnitCell'],
                                   Aimless['ResolutionLow']+'-'+Aimless['ResolutionHigh'],
                                   Aimless['ResolutionLow']+'-'+Aimless['ResolutionLowInnerShell'],
                                   Aimless['ResolutionHighOuterShell']+'-'+Aimless['ResolutionHigh'],
                                   Aimless['RmergeOverall'],
                                   Aimless['RmergeLow'],
                                   Aimless['RmergeHigh'],
                                   Aimless['IsigOverall'],
                                   Aimless['IsigLow'],
                                   Aimless['IsigHigh'],
                                   Aimless['CompletenessOverall'],
                                   Aimless['CompletenessLow'],
                                   Aimless['CompletenessHigh'],
                                   Aimless['MultiplicityOverall'],
                                   Aimless['MultiplicityLow'],
                                   Aimless['MultiplicityHigh'],
                                   Aimless['Alert']  )
                        liststore.append(rowIN)
            else:
                rowIN=('###','###','###','###','###','###','###','###','###','###','###','###','###','###','###','###','###','###','###','###','#FF0000')
                liststore.append(rowIN)   
            treeview =gtk.TreeView(liststore)
            TreeviewListTemp.append(treeview)

                
            for p,ColumnName in enumerate(DiffractionDataColumnName):
                renderer_text = gtk.CellRendererText()
                column_text = gtk.TreeViewColumn(ColumnName)
                column_text.set_sort_column_id(p)
                treeview.append_column(column_text)
                column_text.pack_start(renderer_text,False)
                column_text.add_attribute(renderer_text,"text",p)
                column_text.add_attribute(renderer_text,"cell_background",20)

            # always select the first row (which should be the one with the highest resolution)
            tree_selection = treeview.get_selection()
            tree_selection.select_path(0)

            vboxOverview.add(treeview)
            table.attach(vboxOverview, 1, 2, n, n+1)
            self.DatasetOutcomeList.append(DatasetOutcomeListTemp)
            self.TreeviewList.append(TreeviewListTemp)
            fraction += fractionIncrement
            gobject.idle_add(self.set_progress_bar_fraction, fraction)

        self.data_collection_scrolled.add_with_viewport(table)
        self.window.show_all()
        self.status_bar_message='idle'
        gobject.idle_add(self.set_status_bar_message, self.context_id)
        gobject.idle_add(self.set_progress_bar_fraction, 0.0)
        self.Active=0

    def get_reference_file_list(self):
        # check available reference files
        reference_file_list=[]
        for files in glob.glob(self.reference_directory+'/*'):
            if files.endswith('.pdb'):
                reference_root=files[files.rfind('/')+1:files.rfind('.')]
                if os.path.isfile(self.reference_directory+'/'+reference_root+'.mtz'):
                    mtz_reference=mtztools(self.reference_directory+'/'+reference_root+'.mtz').get_all_values_as_dict()
                    spg_reference=mtz_reference['spacegroup']
                    unitcell_reference=mtz_reference['unitcell']
                    lattice_reference=mtz_reference['bravais_lattice']
                    unitcell_volume_reference=mtz_reference['unitcell_volume']
                    if self.reference_file_root=='':
                        reference_file_list.append([reference_root,
                                                    spg_reference,
                                                    unitcell_reference,
                                                    lattice_reference,
                                                    unitcell_volume_reference])
                    else:
                        if reference_root==self.reference_file_root:
                            reference_file_list.append([reference_root,
                                                        spg_reference,
                                                        unitcell_reference,
                                                        lattice_reference,
                                                        unitcell_volume_reference])
        return reference_file_list


    def InitialRefinement(self):
        reference_file_list=self.get_reference_file_list()
        liststore_reference_files = gtk.ListStore(str)
        for item in reference_file_list:
            liststore_reference_files.append([item[0]])

#        quit()

#                    mtz=mtztools(self.reference_directory+'/'+reference_root+'.mtz')
#                    print mtz.calc_unitcell_volume_from_mtz()
#                    print mtz.get_bravais_lattice_from_mtz()
#                    print mtz.get_spg_from_mtz()
#                    print mtz.get_spg_number_from_mtz()
#                    print mtz.get_unit_cell_from_mtz()


        fraction=0
        fractionIncrement=1/float(self.samples_in_directory(self.InitialModelDir))
        self.status_bar_message='parsing initial_model directory'
        gobject.idle_add(self.set_status_bar_message, self.context_id)
        sample_list=[]
        for n,sample_dir in enumerate(glob.glob(self.InitialModelDir+'/*')):
            run_dimple=True
            resolution_high=''
            Rcryst='pending'
            Rfree='pending'
            spg_autoproc=''
            unitcell_autoproc=''
            spg_reference=''
            unitcell_reference=''
            unitcell_difference=''
            reference=''
            alert='#E0E0E0'
            sample=sample_dir[sample_dir.rfind('/')+1:]
            if os.path.isfile(sample_dir+'/'+sample+'.mtz'):
                mtz_autoproc=mtztools(sample_dir+'/'+sample+'.mtz').get_all_values_as_dict()
                resolution_high=mtz_autoproc['resolution_high']
                spg_autoproc=mtz_autoproc['spacegroup']
                unitcell_autoproc=mtz_autoproc['unitcell']
                lattice_autoproc=mtz_autoproc['bravais_lattice']
                unitcell_volume_autoproc=mtz_autoproc['unitcell_volume']
                # check which reference file is most similar
                for o,reference_file in enumerate(reference_file_list):
                    unitcell_difference=round((math.fabs(reference_file[4]-unitcell_volume_autoproc)/reference_file[4])*100,1)
                    # reference file is accepted when different in unitcell volume < 5%
                    # and both files have the same lattice type
                    if unitcell_difference < 5 and lattice_autoproc==reference_file[3]:
                        spg_reference=reference_file[1]
                        unitcell_reference=reference_file[2]
                        reference=reference_file[0]
                        break
            if os.path.isdir(sample_dir+'/Dimple'):
                    if os.path.isfile(sample_dir+'/Dimple/dimple/final.pdb'):
                        pdb=parse().PDBheader(sample_dir+'/Dimple/dimple/final.pdb')
                        Rcryst=pdb['Rcryst']
                        Rfree=pdb['Rfree']
                        alert=pdb['Alert']
                    elif os.path.isfile(sample_dir+'/dimple_run_in_progress'):
                        Rcryst='in progress'
                        Rfree='in progress'
                        alert='#00CCFF'

            sample_list.append( [ sample,
                                  run_dimple,
                                  resolution_high,
                                  Rcryst,
                                  Rfree,
                                  spg_autoproc,
                                  spg_reference,
                                  unitcell_difference,
                                  unitcell_autoproc,
                                  unitcell_reference,
                                  reference,
                                  alert ] )
                                 
            fraction += fractionIncrement
            gobject.idle_add(self.set_progress_bar_fraction, fraction)

          
        InitalModelTreeviewList=[]
        self.initial_model_liststore = gtk.ListStore(str,bool,str,str,str,str,str,str,str,str,str,str)
        for Sample in sample_list:
            self.initial_model_liststore.append(Sample)
        self.initial_model_treeview =gtk.TreeView(self.initial_model_liststore)
        InitalModelTreeviewList.append(self.initial_model_treeview)
        InitialModelColumnName = [  'SampleID',
                                    'Run\nDimple',
                                    'Resolution',
                                    'Rcryst',
                                    'Rfree',
                                    'Space Group\nautoprocessing',
                                    'Space Group\nreference',
                                    'Difference\nUnit Cell Volume (%)',
                                    'Unit Cell\nautoprocessing',
                                    'Unit Cell\nreference',
                                    'Reference File'    ]
        for n,ColumnName in enumerate(InitialModelColumnName):
            if n == 1:
                col = gtk.TreeViewColumn(ColumnName)
                self.initial_model_treeview.append_column( col)
                cell = gtk.CellRendererToggle()
                cell.set_property( "activatable", True)
                col.pack_start( cell, expand=False)
                col.set_attributes( cell, active=1, cell_background=11)
                col.set_sort_column_id( 1)
                cell.connect('toggled', self._editable_toggled, 1)
            elif n == 10:
                renderer_combo = gtk.CellRendererCombo()
                renderer_combo.set_property("editable", True)
                renderer_combo.set_property("model", liststore_reference_files)
                renderer_combo.set_property("text-column", 0)
                renderer_combo.set_property("has-entry", False)
                renderer_combo.connect("edited", self.on_combo_changed)
                column_combo = gtk.TreeViewColumn(ColumnName, renderer_combo, text=10)
                self.initial_model_treeview.append_column(column_combo)
                column_combo.add_attribute(renderer_combo,"cell_background",11)
            else:
                renderer_text = gtk.CellRendererText()
                column_text = gtk.TreeViewColumn(ColumnName)
                self.initial_model_treeview.append_column(column_text)
                column_text.pack_start(renderer_text,False)
                column_text.add_attribute(renderer_text,"text",n)
                column_text.set_sort_column_id(n)
                column_text.add_attribute(renderer_text,"cell_background",11)

#        self.initial_model_treeview.set_grid_lines(True)             
        self.initial_model_treeview.set_grid_lines(gtk.TREE_VIEW_GRID_LINES_BOTH)
        self.initial_model_scrolled.add(self.initial_model_treeview)
        self.window.show_all()
        self.status_bar_message='idle'
        gobject.idle_add(self.set_status_bar_message, self.context_id)
        gobject.idle_add(self.set_progress_bar_fraction, 0.0)
        self.reference_file_root=''
        self.Active=0

    def on_combo_changed(self, widget, path, text):
        self.initial_model_liststore[path][10] = text

    def _editable_toggled( self, w, row, column):
        self.initial_model_liststore[row][column] = not self.initial_model_liststore[row][column]

    def RunDimple(self):
        self.Active=1
        SampleList=[]
        tree_model=self.initial_model_treeview.get_model()
        for row in tree_model:
            SampleList.append(list(row))
#        for treeview in self.inital_model_treeview:
#            tree_model=treeview.get_model()
#            for row in tree_model:
#                SampleList.append(list(row))

        fraction=0
        fractionIncrement=1/float(len(SampleList))
        for item in SampleList:
            dimple_commands={   'project_directory': self.InitialModelDir,
                                'delete_old': item[1],
                                'xtalID': item[0],
                                'compoundID': '',
                                'smiles': '',
                                'reference': self.reference_directory+'/'+self.reference_file_root  }

            self.status_bar_message='sening jobs to cluster: '+item[0]
            gobject.idle_add(self.set_status_bar_message, self.context_id)
            process(dimple_commands).dimple()
            fraction += fractionIncrement
            gobject.idle_add(self.set_progress_bar_fraction, fraction)

        self.status_bar_message='idle'
        gobject.idle_add(self.set_status_bar_message, self.context_id)
        gobject.idle_add(self.set_progress_bar_fraction, 0.0)
        self.Active=0


    def show_queue_jobs(self):
        job_list=queue().jobs_in_queue()
        if job_list == []:
            job_list=[['---','---',False,'---','---']]
        self.queue_jobs_liststore = gtk.ListStore(str,str,bool,str,str)
        for job in job_list:
            self.queue_jobs_liststore.append(job)
        self.queue_jobs_treeview =gtk.TreeView(self.queue_jobs_liststore)
        queue_column_name = [   'n',
                                'JobID',
                                'delete',
                                'Name',
                                'User'  ]
        for n,column_name in enumerate(queue_column_name):
            if n != 2:
                renderer_text = gtk.CellRendererText()
                column_text = gtk.TreeViewColumn(column_name)
                column_text.set_min_width(120)
                self.queue_jobs_treeview.append_column(column_text)
                column_text.pack_start(renderer_text,False)
                column_text.add_attribute(renderer_text,"text",n)
            else:
                col = gtk.TreeViewColumn(column_name)
                self.queue_jobs_treeview.append_column( col)
                cell = gtk.CellRendererToggle()
#                cell.set_min_width(50)
                cell.set_property( "activatable", True)
                col.pack_start( cell, expand=False)
                col.set_attributes( cell, active=2)
                cell.connect('toggled', self._editable_toggled, 1)
        self.queue_control_scrolled.add(self.queue_jobs_treeview)
        self.window.show_all()
        self.Active=0

    def remove_queue_jobs(self):
        job_list=[]
        tree_model=self.queue_jobs_treeview.get_model()
        for row in tree_model:
            job_list.append(list(row))
        if job_list != []:
            for job in job_list:
                if job[2]==True:
                    print job[1]
        self.Active=0

    def initial_model_choose_reference(self, widget):
        self.reference_file_root=widget.get_active_text()


    def initial_model_toggle_dimple_run(self, widget):
        try:
            tree_model=self.initial_model_treeview.get_model()
            if widget.get_active():
                for n,row in enumerate(tree_model):
                    self.initial_model_liststore[n][1] = True
            else:
                for n,row in enumerate(tree_model):
                    self.initial_model_liststore[n][1] = False                
        except AttributeError:
            widget.set_active(False)

# NOTE: to refresh, try liststore.clear()


#        PDBinfo = { 'Rcryst':         'n/a',
#                    'Rfree':          'n/a',
#                    'SpaceGroup':     'n/a',
#                    'UnitCell':       'n/a',
#                    'ResolutionHigh': 'n/a',
#                    'Alert':          '#E0E0E0' }

#                if os.path.isdir('Dimple')
#                    if final.pdb
#                        # get R/Rfree -> color accordingly

        # Read smiles strings from EXCEL sheet

#        vboxInitialRefinement=gtk.VBox()



#    def Settings(self):
#            vboxSettings=gtk.VBox()





if __name__ == '__main__':
    startGUI=TInsel()
    gtk.main()
