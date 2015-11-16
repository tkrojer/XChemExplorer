import pygtk, gtk, pango
import os
import glob
import sys
import getpass
#sys.path.append('/usr/local/coot/current/lib/python2.7/site-packages')
import coot
#sys.path.append('/usr/local/coot/SoakProc/lib')
#sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
import coot_utils_XChem


class Refine(object):

    def __init__(self,ProjectPath,xtalID,compoundID):
        self.ProjectPath = ProjectPath
        self.xtalID = xtalID
        self.compoundID = compoundID
        self.prefix = 'refine'

    def PrepareForRefinement(self,DataPath):
        # if user decides to keep it all in one directory: never mind!
        if self.ProjectPath==DataPath:
            return None
        # check if xtalID folder exists
        # if no: create folder & link free.mtz + cif/pdb/png file of ligands
        if os.path.isdir(self.ProjectPath+'/'+self.xtalID):
            coot.info_dialog(self.xtalID+' is already under refinement!')
            return None
        else:
            os.mkdir(self.ProjectPath+'/'+self.xtalID)
            os.system('ln -s '+DataPath+'/'+self.xtalID+'/*.pdb '+self.ProjectPath+'/'+self.xtalID)
            os.system('ln -s '+DataPath+'/'+self.xtalID+'/*.mtz '+self.ProjectPath+'/'+self.xtalID)
            os.system('ln -s '+DataPath+'/'+self.xtalID+'/*.map '+self.ProjectPath+'/'+self.xtalID)
            os.system('ln -s '+DataPath+'/'+self.xtalID+'/'+self.compoundID+'.cif '+self.ProjectPath+'/'+self.xtalID)
            os.system('ln -s '+DataPath+'/'+self.xtalID+'/'+self.compoundID+'.png '+self.ProjectPath+'/'+self.xtalID)

    def GetSerial(self):
        # check if there were already previous refinements
        # if no: create a folder Refine_1
        # if yes: create a folder Refine_<max+1>
        temp = []
        found = 0
        if os.path.isdir(self.ProjectPath+'/'+self.xtalID):
            for item in glob.glob(self.ProjectPath+'/'+self.xtalID+'/*'):
                if item.startswith(self.ProjectPath+'/'+self.xtalID+'/Refine_'):
                        temp.append(int(item[item.find('_')+1:]))
                        found = 1
        if found:
            Serial = max(temp) + 1
        else:
            Serial=1
        return Serial


    def RunRefmac(self,Serial,RefmacParams):
        Serial=str(Serial)
        findTLS=''
        TLSphenix=''
        # first check if refinement is ongoing and exit if yes
        if os.path.isfile(self.ProjectPath+'/'+self.xtalID+'/REFINEMENT_IN_PROGRESS'):
            coot.info_dialog('*** REFINEMENT IN PROGRESS ***')
            return None

        RefmacParams['HKLIN']='HKLIN '+self.ProjectPath+'/'+self.xtalID+'/'+self.xtalID+'.free.mtz \\\n'
        RefmacParams['HKLOUT']='HKLOUT '+self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'/refine_'+Serial+'.mtz \\\n'
        RefmacParams['XYZIN']='XYZIN '+self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'/in.pdb \\\n'
        RefmacParams['XYZOUT']='XYZOUT '+self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'/refine_'+Serial+'.pdb \\\n'
        if os.path.isfile(self.ProjectPath+'/'+self.xtalID+'/'+self.compoundID+'.cif'):
            RefmacParams['LIBIN']='LIBIN '+self.ProjectPath+'/'+self.xtalID+'/'+self.compoundID+'.cif \\\n'
            RefmacParams['LIBOUT']='LIBOUT '+self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'/refine_'+Serial+'.cif \\\n'
        if RefmacParams['TLS'].startswith('refi'):
            findTLS='/usr/local/scripts/tobias/SParkle/helpers/FindTLSgroups.py in.pdb\n'
            RefmacParams['TLSIN']='TLSIN '+self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'/refmac.tls \\\n'
            RefmacParams['TLSOUT']='TLSOUT '+self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'/refine.tls \\\n'
            TLSphenix=' phenix.tls '

        # make new refinement folder
        os.mkdir(self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial)

        # now take protein pdb file and write it to newly create Refine_<serial> folder
        # note: the user has to make sure that the ligand file was merged into main file
        for item in coot_utils_SP.molecule_number_list():
            print coot.molecule_name(item),self.prefix+'.pdb'
            if coot.molecule_name(item).endswith(self.prefix+'.pdb'):
                coot.write_pdb_file(item,self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'/in.pdb')

        # we write 'REFINEMENT_IN_PROGRESS' immediately to avoid unncessary refiment
        os.chdir(self.ProjectPath+'/'+self.xtalID)

        os.system('/bin/rm refine.pdb refine.mtz validation_summary.txt validate_ligands.txt 2fofc.map fofc.map refine_molprobity.log')
#        for i in ['refine.pdb','refine.mtz','validation_summary.txt','validate_ligands.txt','2fofc.map','fofc.map']:
#            os.remove(i)
        os.system('touch REFINEMENT_IN_PROGRESS')

        refmacCmds = (
            '#!/bin/csh\n'
            '#PBS -joe -N REFMAC\n'
            'cd '+self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'\n'
            +findTLS+
            'refmac5 '
            +RefmacParams['HKLIN']
            +RefmacParams['HKLOUT']
            +RefmacParams['XYZIN']
            +RefmacParams['XYZOUT']
            +RefmacParams['LIBIN']
            +RefmacParams['LIBOUT']
            +RefmacParams['TLSIN']
            +RefmacParams['TLSOUT']+
            ' << EOF > refmac.log\n'
            'make -\n'
            '    hydrogen ALL -\n'
            '    hout NO -\n'
            '    peptide NO -\n'
            '    cispeptide YES -\n'
            '    ssbridge YES -\n'
            '    symmetry YES -\n'
            '    sugar YES -\n'
            '    connectivity NO -\n'
            '    link NO\n'
            +RefmacParams['NCS']+
            'refi -\n'
            '    type REST -\n'
            '    resi MLKF -\n'
            '    meth CGMAT -\n'
            +RefmacParams['BREF']
            +RefmacParams['TLS']
            +RefmacParams['TWIN']+
            'ncyc '+RefmacParams['NCYCLES']+'\n'
            'scal -\n'
            '    type SIMP -\n'
            '    LSSC -\n'
            '    ANISO -\n'
            '    EXPE\n'
            'solvent YES\n'
            'weight '+RefmacParams['MATRIX_WEIGHT']+'\n'
            'monitor MEDIUM -\n'
            '    torsion 10.0 -\n'
            '    distance 10.0 -\n'
            '    angle 10.0 -\n'
            '    plane 10.0 -\n'
            '    chiral 10.0 -\n'
            '    bfactor 10.0 -\n'
            '    bsphere 10.0 -\n'
            '    rbond 10.0 -\n'
            '    ncsr 10.0\n'
            'labin  FP=F SIGFP=SIGF FREE=FreeR_flag\n'
            'labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM\n'
            +RefmacParams['TLSADD']+'\n'
            'DNAME '+self.xtalID+'\n'
            'END\n'
            'EOF\n'
            '\n'
            'phenix.molprobity refine_%s.pdb refine_%s.mtz\n' %(Serial,Serial)+
            '/bin/mv molprobity.out refine_molprobity.log\n'
            'mmtbx.validate_ligands refine_%s.pdb refine_%s.mtz LIG > validate_ligands.txt\n' %(Serial,Serial)+
            'cd '+self.ProjectPath+'/'+self.xtalID+'\n'
            'ln -s %s/%s/Refine_%s/refine_%s.pdb refine.pdb\n' %(self.ProjectPath,self.xtalID,Serial,Serial)+
            'ln -s %s/%s/Refine_%s/refine_%s.mtz refine.mtz\n' %(self.ProjectPath,self.xtalID,Serial,Serial)+
            'ln -s Refine_%s/validate_ligands.txt .\n' %Serial+
            'ln -s Refine_%s/refine_molprobity.log .\n' %Serial+
            'mmtbx.validation_summary refine.pdb > validation_summary.txt\n'
            '\n'
            'fft hklin refine.mtz mapout 2fofc.map << EOF\n'
            'labin F1=FWT PHI=PHWT\n'
            'EOF\n'
            '\n'
            'fft hklin refine.mtz mapout fofc.map << EOF\n'
            'labin F1=DELFWT PHI=PHDELWT\n'
            'EOF\n'
             '\n'
            '/bin/rm %s/%s/REFINEMENT_IN_PROGRESS\n' %(self.ProjectPath,self.xtalID)+
            '\n'
            '#cd '+self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'\n'
            '#\n'
            '#phenix.refine in.pdb ../'+self.xtalID+'.free.mtz ../'+self.compoundID+'.cif '+TLSphenix+' '
            '# refinement.input.xray_data.labels=../'+self.xtalID+'.free.mtz:F,SIGF '
            '# optimize_xyz_weight=true '
            '# optimize_adp_weight=true ordered_solvent=False output.prefix=refine_phenix_'+Serial+'\n'
            '#\n' 
            '#phenix.molprobity refine_phenix_%s.pdb refine_phenix_%s.mtz\n' %(Serial,Serial)+
            '#/bin/mv molprobity.out refine_phenix_molprobity.log\n'
            '#mmtbx.validate_ligands refine_phenix_%s.pdb refine_phenix_%s.mtz LIG > validate_ligands_phenix.txt\n' %(Serial,Serial)+
            '#\n'
           )

        cmd = open(self.ProjectPath+'/'+self.xtalID+'/Refine_'+Serial+'/refmac.csh','w')
        cmd.write(refmacCmds)
        cmd.close()
        print '\n==> cootSP: writing input script for refinement'
        os.system('ssh %s@artemis "cd %s/%s/Refine_%s; qsub refmac.csh"' %(getpass.getuser(),self.ProjectPath,self.xtalID,Serial))
        print '\n==> cootSP: submitting job to cluster'



    def RefinementParams(self,RefmacParams):
        self.RefmacParams=RefmacParams
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        self.window.set_title("Refmac Parameters")
        self.vbox = gtk.VBox()

        self.hbox1=gtk.HBox()
        self.hbox1.add(gtk.Label('Refine'))
        self.cb = gtk.combo_box_new_text()
        self.cb.connect("changed", self.ChooseBfacRefinement)
        for item in ['isotropic','anisotropic']:
            self.cb.append_text(item)
        if 'ISOT' in self.RefmacParams['BREF']:
            self.cb.set_active(0)
        if 'ANIS' in self.RefmacParams['BREF']:
            self.cb.set_active(1)
        self.hbox1.add(self.cb)
        self.hbox1.add(gtk.Label('temperature factors'))
        self.vbox.add(self.hbox1)

        self.hbox2=gtk.HBox()
        self.hbox2.add(gtk.Label('Number of Cycles: '))
        self.Ncycles=gtk.Entry()
        self.Ncycles.add_events(gtk.gdk.KEY_RELEASE_MASK)
        self.Ncycles.connect("key-release-event", self.on_key_release_Ncycles)
        self.Ncycles.set_text(self.RefmacParams['NCYCLES'])
        self.hbox2.add(self.Ncycles)
        self.vbox.add(self.hbox2)

        self.hbox3=gtk.HBox()
        self.hbox3.add(gtk.Label('MATRIX WEIGHT: '))
        self.MATRIX_WEIGHT=gtk.Entry()
        self.MATRIX_WEIGHT.add_events(gtk.gdk.KEY_RELEASE_MASK)
        self.MATRIX_WEIGHT.connect("key-release-event", self.on_key_release_MATRIX_WEIGHT)
        self.MATRIX_WEIGHT.set_text(self.RefmacParams['MATRIX_WEIGHT'])
        self.hbox3.add(self.MATRIX_WEIGHT)
        self.vbox.add(self.hbox3)

        self.TLS = gtk.CheckButton('TLS (find TLS groups with phenix.find_tls_groups)')
        self.TLS.connect("toggled", self.TLSCallback)
        if self.RefmacParams['TLS']=='refi tlsc 10\n': self.TLS.set_active(True)
        self.vbox.pack_start(self.TLS,False)

        self.NCS = gtk.CheckButton('NCS (if applicable')
        self.NCS.connect("toggled", self.NCSCallback)
        if self.RefmacParams['NCS']=='NCSR LOCAL\n': self.NCS.set_active(True)
        self.vbox.pack_start(self.NCS,False)

        self.TWIN = gtk.CheckButton('Twin?')
        self.TWIN.connect("toggled", self.TWINCallback)
        if self.RefmacParams['TWIN']=='TWIN\n': self.TWIN.set_active(True)
        self.vbox.pack_start(self.TWIN,False)

        self.OKbutton = gtk.Button(label="OK")
        self.OKbutton.connect("clicked",self.OK)
        self.vbox.add(self.OKbutton)

        self.window.add(self.vbox)
        self.window.show_all()
        return self.RefmacParams


    def TLSCallback(self, widget):
        if widget.get_active():
            self.RefmacParams['TLS']='refi tlsc 10\n'
            self.RefmacParams['TLSIN']='refmac.tls\n'
            self.RefmacParams['TLSOUT']='out.tls\n'
            self.RefmacParams['TLSADD']='TLSO ADDU\n'
        else:
            self.RefmacParams['TLS']=''
            self.RefmacParams['TLSIN']=''
            self.RefmacParams['TLSOUT']=''
            self.RefmacParams['TLSADD']=''
        return self.RefmacParams

    def NCSCallback(self, widget):
        if widget.get_active():
            self.RefmacParams['NCS']='NCSR LOCAL\n'
        else:
            self.RefmacParams['NCS']=''
        return self.RefmacParams

    def ChooseBfacRefinement(self,widget):
        if widget.get_active_text()=='isotropic':
            self.RefmacParams['BREF']='    bref ISOT\n'
        if widget.get_active_text()=='anisotropic':
            self.RefmacParams['BREF']='    bref ANIS\n'
        return self.RefmacParams

    def on_key_release_Ncycles(self, widget, event):
        print widget.get_text()
        self.RefmacParams['NCYCLES'] = widget.get_text()
        return self.RefmacParams

    def on_key_release_MATRIX_WEIGHT(self, widget, event):
        self.RefmacParams['MATRIX_WEIGHT'] = widget.get_text()
        return self.RefmacParams

    def TWINCallback(self, widget):
        if widget.get_active():
            self.RefmacParams['TWIN']='TWIN\n'
        else:
            self.RefmacParams['TWIN']=''
        return self.RefmacParams

    def OK(self,widget):
        self.window.destroy()


    def ParamsFromPreviousCycle(self,Serial):

        RefmacParams={ 'HKLIN': '', 'HKLOUT': '',
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

        if os.path.isfile(self.ProjectPath+'/'+self.xtalID+'/Refine_'+str(Serial)+'/refmac.csh'):
            for line in open(self.ProjectPath+'/'+self.xtalID+'/Refine_'+str(Serial)+'/refmac.csh'):
                if line.startswith('refi tlsc'):
                    RefmacParams['TLS']=line
                if line.startswith('TLSO'):
                    RefmacParams['TLSADD']=line
                if line.startswith('NCSR LOCAL'):
                    RefmacParams['NCS']=line
                if line.startswith('    bref '):
                    RefmacParams['BREF']=line
                if line.startswith('ncyc'):
                    RefmacParams['Ncycles'] = line.split()[1]
                if line.startswith('weight'):
                    RefmacParams['MATRIX_WEIGHT'] = line[7:-1]
                if line.startswith('TWIN'):
                    RefmacParams['TWIN']=line

        return RefmacParams

    def GetRefinementHistory(self):
        RefinementHistory=''
        RefinementCycle = []
        found = 0
        for item in glob.glob(self.ProjectPath+'/'+self.xtalID+'/*'):
            if item.startswith(self.ProjectPath+'/'+self.xtalID+'/Refine_'):
                    RefinementCycle.append(int(item[item.find('_')+1:]))
                    found = 1
        if found:
            for cycle in sorted(RefinementCycle):
                R=0
                Rfree=0
                LigandCC=0
                try:
                    newestPDB = max(glob.iglob(self.ProjectPath+'/'+self.xtalID+'/Refine_'+str(cycle)+'/refine_'+str(cycle)+'.pdb'), key=os.path.getctime)
                    for line in open(newestPDB):
                        if line.startswith('REMARK   3   R VALUE     (WORKING + TEST SET) :'): R = line.split()[9]
                        if line.startswith('REMARK   3   FREE R VALUE                     :'): Rfree = line.split()[6]
                    if os.path.isfile(self.ProjectPath+'/'+self.xtalID+'/Refine_'+str(cycle)+'/validate_ligands.txt'):
                        for line in open(self.ProjectPath+'/'+self.xtalID+'/Refine_'+str(cycle)+'/validate_ligands.txt'):
                                if line.startswith('|  LIG'): LigandCC = line.split()[6]
                    RefinementHistory=RefinementHistory+str(cycle).rjust(10)+str(R).rjust(10)+str(Rfree).rjust(10)+str(LigandCC).rjust(10)+'\n'
                except ValueError:
                    RefinementHistory=RefinementHistory+str(cycle).rjust(10)+str(R).rjust(10)+str(Rfree).rjust(10)+str(LigandCC).rjust(10)+'\n'
            f=open(self.ProjectPath+'/'+self.xtalID+'/RefinementHistory.txt','w')
            f.write(RefinementHistory)
            f.close()
            gnuIN = (   'set term png size 300,225\n'
                        'set output "RefinementHistory.png"\n'
                        'set key right bottom\n'
                        'set xlabel "Cycle"\n'
                        'set ylabel "R/Rfree"\n'
                        'set yrange [0:0.35]\n'
                        'set y2label "ligand CC"\n'
                        'set y2range [0:1]\n'
                        'set ytics nomirror\n'
                        'set y2tics\n'
                        'set tics out\n'
                        'set xtics 1\n'
                        'set style line 1 lt 1 lw 2 linecolor rgb "red"\n'
                        'set style line 2 lt 1 lw 2 linecolor rgb "green"\n'
                        'set style line 3 lt 1 lw 2 linecolor rgb "blue"\n'
                        'plot "RefinementHistory.txt" using 1:2  title "R" axes x1y1 with line ls 1,'
                        '     "RefinementHistory.txt" using 1:3  title "Rfree" axes x1y1 with line ls 2,'
                        '     "RefinementHistory.txt" using 1:4  title "ligand CC" axes x2y2 with line ls 3\n'
            )
            f=open(self.ProjectPath+'/'+self.xtalID+'/RefinementHistory.gnu','w')
            f.write(gnuIN)
            f.close()
            os.chdir(self.ProjectPath+'/'+self.xtalID)
            os.system('gnuplot RefinementHistory.gnu')








