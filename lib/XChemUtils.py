#!/usr/bin/python

import sys
import os
import glob
import math
import subprocess
import getpass


class process:

    def __init__(self,dimple):
        self.project_directory=dimple['project_directory']
        self.delete_old=dimple['delete_old']
        self.xtalID=dimple['xtalID']
        self.compoundID=dimple['compoundID']
        self.smiles=dimple['smiles']
        self.reference=dimple['reference']

        self.mtz_free=self.xtalID+'.free.mtz'
        # need to set $CCP4_SCR, because the default directory in /home/zqr16691/tmp fills up easily 
        # and then dimple will not run
        if not os.path.isdir(self.project_directory[:self.project_directory.find('processing')+11]+'/tmp'):
            os.mkdir(self.project_directory[:self.project_directory.find('processing')+11]+'/tmp')
        self.ccp4_scratch=self.project_directory[:self.project_directory.find('processing')+11]+'/tmp'


    def get_Rfree(self):
        if not os.path.isfile(self.project_directory+'/'+self.xtalID+'/'+self.mtz_free):
            Cmds = (
                    '#!/bin/csh\n'
                    '\n'
                    'cd %s/%s\n' %(self.project_directory,self.xtalID) +
                    '\n'
                    'pointless hklin %s.mtz hklref %s.mtz hklout %s.reind.mtz << EOF > pointless.log\n' %(self.xtalID,self.reference,self.xtalID)+
                    ' tolerance 5\n'
                    'EOF\n'
                    '\n'
                    'cad hklin1 %s.reind.mtz hklin2 %s.mtz hklout cad.mtz << EOF > cad.log\n' %(self.xtalID,self.reference) +
                    '   labin file_number 1 E1=F E2=SIGF E3=IMEAN E4=SIGIMEAN\n'
                    '   labout file_number 1 E1=F E2=SIGF E3=IMEAN E4=SIGIMEAN\n'
                    '   labin file_number 2 E1=FreeR_flag\n'
                    '   labout file_number 2 E1=FreeR_flag\n'
                    '   END\n'
                    'EOF\n'
                    '\n'
                    'freerflag hklin cad.mtz hklout %s.free.mtz << EOF > freerflag.log\n' %self.xtalID +
                    '   COMPLETE FREE=FreeR_flag\n'
                    '   END\n'
                    'EOF\n'
                    )
            os.chdir('%s/%s' %(self.project_directory,self.xtalID))
            os.system(Cmds)


    def dimple(self):
        os.chdir('%s/%s' %(self.project_directory,self.xtalID))
        if os.path.isdir('%s/%s/Dimple' %(self.project_directory,self.xtalID)):
            # this 'if not' step is usually not necessary; only if process was stopped in an odd way
            if not os.path.isfile(self.project_directory+'/'+self.xtalID+'/'+self.mtz_free):
                self.get_Rfree()
            if self.delete_old==True:
                os.system('rm -fr Dimple')
                os.mkdir('%s/%s/Dimple' %(self.project_directory,self.xtalID))
            else:
                return None
        else:
            os.mkdir('%s/%s/Dimple' %(self.project_directory,self.xtalID))
            self.get_Rfree()
            os.chdir('%s/%s/Dimple' %(self.project_directory,self.xtalID))

        if self.smiles != []:
#            ligand='ligand_smiles="%s"' %self.smiles[0][0]
            ligand="acedrg -i '%s' -r LIG -o %s\n" %(self.smiles,self.compoundID)
        else:
            ligand=''
        Cmds = (
                '#!/bin/csh\n'
                '#PBS -joe -N test\n'
                '\n'
                'cd %s/%s/Dimple\n' %(self.project_directory,self.xtalID) +
                '\n'
                '#module load ccp4/6.5.007\n'
                'module load ccp4\n'
                '#ccp4/6.5-update1\n'
                '\n'
                'export CCP4_SCR='+self.ccp4_scratch+'\n'
                '\n'
                'dimple ../%s.free.mtz %s.pdb dimple' %(self.xtalID,self.reference) +
                '\n'
                'cd dimple/dimple\n'
                '\n'
                'fft hklin final.mtz mapout 2fofc.map << EOF\n'
                ' labin F1=2FOFCWT PHI=PH2FOFCWT\n'
                'EOF\n'
                '\n'
                'fft hklin final.mtz mapout fofc.map << EOF\n'
                ' labin F1=FOFCWT PHI=PH2FOFCWT\n'
                'EOF\n'
                '\n'
                '%s\n' %ligand +
                '\n'
                '#ln -s ../../%s.png .\n' %self.compoundID+
                '\n'
                'cd %s/%s\n' %(self.project_directory,self.xtalID) +
                '\n'
                '/bin/rm dimple_run_in_progress\n' 
                )
        self.Run(Cmds)



    def Run(self,Cmds):
        os.chdir(self.project_directory+'/'+self.xtalID)
        os.system('touch dimple_run_in_progress')
        os.chdir(self.project_directory+'/'+self.xtalID+'/Dimple')
        f = open('%s.sh' %self.xtalID,'w')
        f.write(Cmds)
        f.close()
        os.system('qsub %s.sh' %self.xtalID)

#def MakePng(DataPath,xtalID,compoundID,smiles):
#    os.chdir('%s/%s' %(DataPath,xtalID))
#    os.system('rm -f %s.png' %smiles)
#    Cmds = (
#            'echo "xtalID: %s"\n' %xtalID +
#            '\n'
#            'echo "======== generating 2D image of compound ======="\n'
#            '\n'
#            'module load openbabel\n'
#            '\n'
#            'obabel -:"%s" -O temp.png' %smiles+
#            '\n'
#            'convert temp.png -background white label:"%s" -gravity Center -append %s.png\n' %(compoundID,compoundID)+
#            '/bin/rm temp.png\n'
#            )
#    os.system(Cmds)

class parse:

    def GetAimlessLog(self,Logfile):
        self.Logfile=Logfile
        Aimless = { 'AutoProc': 'n/a',
                    'Run': 'n/a',
                    'SpaceGroup': 'n/a',
                    'UnitCell': 'n/a',
                    'ResolutionLow': 'n/a',
                    'ResolutionLowInnerShell': 'n/a',
                    'ResolutionHigh': 'n/a',
                    'ResolutionHighOuterShell': 'n/a',
                    'RmergeOverall': 'n/a',
                    'RmergeLow': 'n/a',
                    'RmergeHigh': 'n/a',
                    'IsigOverall': 'n/a',
                    'IsigLow': 'n/a',
                    'IsigHigh': 'n/a',
                    'CompletenessOverall': 'n/a',
                    'CompletenessLow': 'n/a',
                    'CompletenessHigh': 'n/a',
                    'MultiplicityOverall': 'n/a',
                    'MultiplicityLow': 'n/a',
                    'MultiplicityHigh': 'n/a',
                    'Alert': '#FF0000' }

        spg='n/a'
        a='n/a'
        b='n/a'
        c='n/a'
        alpha='n/a'
        beta='n/a'
        gamma='n/a'

        if 'fast_dp' in self.Logfile:
            Aimless['AutoProc']='fast_dp'
        if '3d-run' in self.Logfile:
            Aimless['AutoProc']='xia2 3d'
        if '3dii-run' in self.Logfile:
            Aimless['AutoProc']='xia2 3dii'
        if 'dials-run' in self.Logfile:
            Aimless['AutoProc']='dials'

        # get run number from logfile
        # Note: only works if file is in original directory, but not once it lifes in 'inital_model'
#        print self.Logfile.split('/')[9].split('_')[1]
#        if len(self.Logfile.split('/'))>8 and len(self.Logfile.split('/')[9].split('_'))==1:
        try:
            Aimless['Run']=self.Logfile.split('/')[9].split('_')[1]
        except IndexError:
            pass

        for line in open(self.Logfile):
            if line.startswith('Low resolution limit') and len(line.split())==6:
                Aimless['ResolutionLow'] = line.split()[3]
                Aimless['ResolutionHighOuterShell'] = line.split()[5]
            if line.startswith('High resolution limit') and len(line.split())==6:
                Aimless['ResolutionHigh'] = line.split()[3]
                Aimless['ResolutionLowInnerShell'] = line.split()[4]
            if line.startswith('Rmerge  (all I+ and I-)') and len(line.split())==8:
                Aimless['RmergeOverall'] = line.split()[5]
                Aimless['RmergeLow'] = line.split()[6]
                Aimless['RmergeHigh']  = line.split()[7]
            if line.startswith('Mean((I)/sd(I))') and len(line.split())==4:
                Aimless['IsigOverall'] = line.split()[1]
                Aimless['IsigHigh'] = line.split()[3]
                Aimless['IsigLow'] = line.split()[2]
            if line.startswith('Completeness') and len(line.split())==4:
                Aimless['CompletenessOverall'] = line.split()[1]
                Aimless['CompletenessHigh'] = line.split()[3]
                Aimless['CompletenessLow'] = line.split()[2]
            if line.startswith('Multiplicity') and len(line.split())==4:
                Aimless['MultiplicityOverall'] = line.split()[1]
                Aimless['MultiplicityHigh'] = line.split()[3]
                Aimless['MultiplicityLow'] = line.split()[3]
            if line.startswith('Average unit cell:') and len(line.split())==9:
                tmp = []
                tmp.append(line.split())
                a = int(float(tmp[0][3]))
                b = int(float(tmp[0][4]))
                c = int(float(tmp[0][5]))
                alpha = int(float(tmp[0][6]))
                beta = int(float(tmp[0][7]))
                gamma = int(float(tmp[0][8]))
            if line.startswith('Space group:'):
                Aimless['SpaceGroup']=line.replace('Space group: ','')[:-1]

        Aimless['UnitCell']=str(a)+' '+str(b)+' '+str(c)+' '+str(alpha)+' '+str(beta)+' '+str(gamma)

        # Hex Color code:
        # red:      #FF0000
        # orange:   #FF9900
        # green:    #00FF00
        # gray:     #E0E0E0

        if Aimless['ResolutionHigh']=='n/a' or Aimless['RmergeLow'] =='n/a':
            Aimless['Alert'] = '#FF0000'
        else:
            if float(Aimless['ResolutionHigh']) > 3.5 or float(Aimless['RmergeLow']) > 0.1:
                Aimless['Alert'] = '#FF0000'
            if (float(Aimless['ResolutionHigh']) <= 3.5 and float(Aimless['ResolutionHigh']) > 2.5) or \
               (float(Aimless['RmergeLow']) <= 0.1 and float(Aimless['RmergeLow']) > 0.05):
                Aimless['Alert'] = '#FF9900'
            if float(Aimless['ResolutionHigh']) <= 2.5 and float(Aimless['RmergeLow']) <= 0.05:
                Aimless['Alert'] = '#00FF00'

        return Aimless

    def PDBheader(self,pdbfile):
        PDBinfo = { 'Rcryst':         'n/a',
                    'Rfree':          'n/a',
                    'SpaceGroup':     'n/a',
                    'UnitCell':       'n/a',
                    'ResolutionHigh': 'n/a',
                    'Alert':          '#E0E0E0' }
        if os.path.isfile(pdbfile):
            for line in open(pdbfile):
                if line.startswith('REMARK   3   R VALUE     (WORKING + TEST SET) :'):
                    PDBinfo['Rcryst']=line.split()[9]
                    if float(PDBinfo['Rcryst']) > 0.4:
                        PDBinfo['Alert']='#FF0000'
                    if float(PDBinfo['Rcryst']) <= 0.4 and float(PDBinfo['Rcryst']) >= 0.3:
                        PDBinfo['Alert']='#FF9900'
                    if float(PDBinfo['Rcryst']) < 0.3:
                        PDBinfo['Alert']='#00FF00'
                if line.startswith('REMARK   3   FREE R VALUE                     :'):
                    PDBinfo['Rfree']=line.split()[6]
                    if float(PDBinfo['Rfree']) > 0.4:
                        PDBinfo['Alert']='#FF0000'
                    if float(PDBinfo['Rfree']) <= 0.4 and float(PDBinfo['Rfree']) >= 0.3:
                        PDBinfo['Alert']='#FF9900'
                    if float(PDBinfo['Rfree']) < 0.3:
                        PDBinfo['Alert']='#00FF00'
                if line.startswith('REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) :'):
                    PDBinfo['ResolutionHigh']=line.split()[7]
                if line.startswith('CRYST1'):
                    PDBinfo['UnitCell']=line.split()[1]+' '+line.split()[2]+' '+line.split()[3]+' '+ \
                                        line.split()[4]+' '+line.split()[5]+' '+line.split()[6]
                    PDBinfo['SpaceGroup']=line[55:len(line)-1].replace(' ','')
        return PDBinfo

class mtztools:
    
    def __init__(self,mtzfile):
        self.mtzfile=mtzfile
        self.space_group_dict=   {  'triclinic':    [1],
                                    'monoclinic':   [3,4,5],
                                    'orthorhombic': [16,17,18,19,20,21,22,23,24],
                                    'tetragonal':   [75,76,77,78,79,80,89,90,91,92,93,94,95,96,97,98],
                                    'hexagonal':    [143,144,145,149,150,151,152,153,154,168,169,170,
                                                     171,172,173,177,178,179,180,181,182],
                                    'rhombohedral': [146,155],
                                    'cubic':        [195,196,197,198,199]  }

    def get_bravais_lattice_from_spg_number(self,number):
        lattice=''
        for bravaislattice in self.space_group_dict:
            for spg_number in self.space_group_dict[bravaislattice]:
                if spg_number==number:
                    lattice=bravaislattice
        return lattice

    def get_unit_cell_from_mtz(self):
        unitcell=[]
        cell_line=100000
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' * Cell Dimensions :'):
                cell_line=n+2
            if n==cell_line and len(line.split())==6:
                a=line.split()[0]
                b=line.split()[1]
                c=line.split()[2]
                alpha=line.split()[3]
                beta=line.split()[4]
                gamma=line.split()[5]
        return [a,b,c,alpha,beta,gamma]

    def get_spg_number_from_mtz(self):
        spg_number=0
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' * Space group ='):
                spg_number=int(line[line.rfind(')')-3:line.rfind(')')])
        return spg_number

    def get_spg_from_mtz(self):
        spg=''
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' * Space group ='):
                spg=line[line.find("'")+1:line.rfind("'")]        
        return spg

    def get_bravais_lattice_from_mtz(self):
        bravais_lattice=''
        spg_number=self.get_spg_number_from_mtz()
        bravais_lattice=self.get_bravais_lattice_from_spg_number(spg_number)
        return bravais_lattice

    def calc_unitcell_volume_from_mtz(self):
        spg_number=self.get_spg_number_from_mtz()
        lattice=self.get_bravais_lattice_from_spg_number(spg_number)
        unitcell=self.get_unit_cell_from_mtz()
        a=float(unitcell[0])
        b=float(unitcell[1])
        c=float(unitcell[2])
        alpha=math.radians(float(unitcell[3]))
        beta= math.radians(float(unitcell[4]))
        gamma=math.radians(float(unitcell[5]))
        unitcell_volume=0
        if lattice=='triclinic':
            unitcell_volume=a*b*c* \
                            math.sqrt((1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2) \
                                      +2*(math.cos(alpha)*math.cos(beta)*math.cos(gamma)))
        if lattice=='monoclinic':
            unitcell_volume=round(a*b*c*math.sin(beta),1)
        if lattice=='orthorhombic' or lattice=='tetragonal' or lattice=='cubic':
            unitcell_volume=round(a*b*c,1)
        if lattice=='hexagonal' or lattice=='rhombohedral':
            unitcell_volume=round(a*b*c*(math.sin(math.radians(60))),1)
        return unitcell_volume

    def get_high_resolution_from_mtz(self):
        resolution_high='n/a'
        resolution_line=1000000
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' *  Resolution Range :'):
                resolution_line=n+2
            if n==resolution_line and len(line.split())==8:
                resolution_high=line.split()[5]
        return resolution_high
        
    def get_all_values_as_dict(self):
        mtz = { 'resolution_high':  'n/a',
                'unitcell':         'n/a',
                'spacegroup':       'n/a',
                'unitcell_volume':  'n/a',
                'bravais_lattice':  'n/a'   }
        a=0.0
        b=0.0
        c=0.0
        alpha_rad=0.0
        beta_rad=0.0
        gamma_rad=0.0
        resolution_line=1000000
        cell_line=1000000
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' *  Resolution Range :'):
                resolution_line=n+2
            if n==resolution_line and len(line.split())==8:
                mtz['resolution_high']=round(float(line.split()[5]),2)
            if line.startswith(' * Cell Dimensions :'):
                cell_line=n+2
            if n==cell_line and len(line.split())==6:
                a=      round(float(line.split()[0]),1)
                b=      round(float(line.split()[1]),1)
                c=      round(float(line.split()[2]),1)
                alpha=  round(float(line.split()[3]),1)
                beta=   round(float(line.split()[4]),1)
                gamma=  round(float(line.split()[5]),1)
                mtz['unitcell']=str(a)+' '+str(b)+' '+str(c)+' '+ \
                                str(alpha)+' '+str(beta)+' '+str(gamma)
                alpha_rad=math.radians(alpha)
                beta_rad= math.radians(beta)
                gamma_rad=math.radians(gamma)
            if line.startswith(' * Space group ='):
                spg_number=int(line[line.rfind(')')-3:line.rfind(')')])
                mtz['bravais_lattice']=self.get_bravais_lattice_from_spg_number(spg_number)
                mtz['spacegroup']=line[line.find("'")+1:line.rfind("'")]
        if mtz['bravais_lattice']=='triclinic':
            mtz['unitcell_volume']=a*b*c* \
                            math.sqrt((1-math.cos(alpha_rad)**2-math.cos(beta_rad)**2-math.cos(gamma_rad)**2) \
                                      +2*(math.cos(alpha_rad)*math.cos(beta_rad)*math.cos(gamma_rad)))
        elif mtz['bravais_lattice']=='monoclinic':
            mtz['unitcell_volume']=round(a*b*c*math.sin(beta_rad),1)
        elif mtz['bravais_lattice']=='orthorhombic' or mtz['bravais_lattice']=='tetragonal' or mtz['bravais_lattice']=='cubic':
            mtz['unitcell_volume']=round(a*b*c,1)
        elif mtz['bravais_lattice']=='hexagonal' or mtz['bravais_lattice']=='rhombohedral':
            mtz['unitcell_volume']=round(a*b*c*(math.sin(math.radians(60))),1)

        return mtz

class queue:

    def jobs_in_queue(self):
            jobs=[]
            qstat=subprocess.Popen(['qstat -r'],stdout=subprocess.PIPE)
            counter=0
            for n,line in enumerate(iter(qstat.stdout.readline,'')):
                jobID=''
                job_name=''
                user_name=''
                if getpass.getuser() in line:
                    counter +=1
                    jobID=line.split()[0]
                    user_name=line.split()[3]
                if "Full jobname:" in line:
                    job_name=line.split()[2]
                jobs.append([counter,jobID,True,job_name,user_name])
            return jobs

    def remove_all_jobs_from_queue(self):
            jobs_in_queue=self.jobs_in_queue()
            if jobs_in_queue != []:
                for job in jobs_in_queue:
                    os.system('qdel '+job[0])





