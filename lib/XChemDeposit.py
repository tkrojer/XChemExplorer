# last edited: 08/02/2017, 15:00

import sys
import os
import glob

from PyQt4 import QtGui, QtCore, QtWebKit

from lxml import etree

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
import XChemLog
import XChemDB
import XChemToolTips
import XChemMain
from XChemUtils import pdbtools
from XChemUtils import smilestools


def create_SF_mmcif(outDir,mtzList):
    print 'hallo'

def get_protein_sequence(database,xtalID):
    print 'hallo'

def check_depositDict(depositDict):
    # check depositDict
    for entry in depositDict:
        if 'middle_name' in depositDict[entry]:
            continue
        elif 'State_or_Province' in depositDict[entry]:
            continue
        elif depositDict[entry] == '':
            print 'ERROR'

def update_title(depositDict):

    print 'hallo'

def create_data_template_text():

    data_template_text=data_template(depositDict,sequence)

def create_Model_mmcif(outDir,pdbList):
    print 'hallo'

def update_file_locations_of_apo_structuresin_DB(database,projectDir,xce_logfile):
    Logfile=XChemLog.updateLog(xce_logfile)
    Logfile.insert('updating file information for apo structures')
    db=XChemDB.data_source(database)
    apo=db.execute_statement("select CrystalName from depositTable where StructureType is 'apo';")
    for item in apo:
        xtal=str(item[0])
        db_dict={}
        db_dict['label']=xtal+'-apo'
        db_dict['description']='apo structure for pandda.analyse'
        if os.path.isfile(os.path.join(projectDir,xtal,'dimple.pdb')):
            db_dict['PDB_file']=os.path.realpath(os.path.join(projectDir,xtal,'dimple.pdb'))
            if os.path.isfile(os.path.join(projectDir,xtal,'dimple.mtz')):
                db_dict['PDB_file']=os.path.realpath(os.path.join(projectDir,xtal,'dimple.mtz'))
                Logfile.insert('updating depositTable for apo structure '+xtal)
                db.update_depositTable(xtal,'apo',db_dict)





class templates:

    def data_template_cif(self,depositDict):

        taxonomy_dict=XChemMain.NCBI_taxonomy_ID()
        for key in taxonomy_dict:
            if taxonomy_dict[key]==depositDict['Source_organism_scientific_name']:
                pdbx_gene_src_ncbi_taxonomy_id=key
            if taxonomy_dict[key]==depositDict['Expression_system_scientific_name']:
                pdbx_host_org_ncbi_taxonomy_id=key

        structure_author_name=''
        for name in depositDict['structure_author_name'].split(';'):
            structure_author_name+='<structure_author_name=  %s>\n' %name

        primary_citation_author_name=''
        # one name must be within quotation, last name and first initial must be separated by comma and space
        for name in depositDict['primary_citation_author_name'].split(';'):
            if name.replace(' ','') == '':
                continue
            if name[name.find(',')+1:name.find(',')+2] != ' ':
                name=name.replace(',',', ')
            primary_citation_author_name+="primary '%s'\n" %name

        molecule_one_letter_sequence=';'
        counter=1
        for aa in depositDict['molecule_one_letter_sequence']:
            if counter < 70:
                molecule_one_letter_sequence+=aa
            if counter == 70:
                molecule_one_letter_sequence+='\n'+aa
                counter = 0
            counter+=1

        data_template_cif = (
            'data_UNNAMED\n'
            '#\n'
            '_pdbx_database_status.entry_id                       UNNAMED\n'
            "_pdbx_database_status.dep_release_code_coordinates   '%s'\n"                       %depositDict['Release_status_for_coordinates']+
            "_pdbx_database_status.dep_release_code_sequence      '%s'\n"                       %depositDict['Release_status_for_sequence']+
            '#\n'
            '_pdbx_deposit_group.group_id	   UNNAMED\n'
            '_pdbx_deposit_group.group_description  "%s"\n'                                     %depositDict['group_description']+
            '_pdbx_deposit_group.group_title        "%s"\n'                                     %depositDict['group_title']+
            '#\n'
            '_exptl_crystal_grow.crystal_id      1\n'
            "_exptl_crystal_grow.method          '%s'\n"                                        %depositDict['crystallization_method']+
            '_exptl_crystal_grow.pH              %s\n'                                          %depositDict['crystallization_pH']+
            '_exptl_crystal_grow.temp            %s\n'                                          %depositDict['crystallization_temperature']+
            '_exptl_crystal_grow.pdbx_details    "%s"\n'                                        %depositDict['crystallization_details']+
            '#\n'
            '_diffrn.id                     1\n'
            '_diffrn.ambient_temp           %s\n'                                               %depositDict['data_collection_temperature']+
            '_diffrn.crystal_id             1\n'
            '#\n'
            '_diffrn_source.diffrn_id                       1\n'
            '_diffrn_source.source                          %s\n'                               %depositDict['radiation_source']+
            '_diffrn_source.type                            "%s"\n'                             %depositDict['radiation_source_type']+
            '_diffrn_source.pdbx_wavelength_list            %s\n'                               %depositDict['radiation_wavelengths']+
            '#\n'
            '_diffrn_detector.detector               %s\n'                                      %depositDict['radiation_detector']+
            "_diffrn_detector.type                   '%s'\n"                                    %depositDict['radiation_detector_type']+
            '_diffrn_detector.pdbx_collection_date   %s\n'                                      %depositDict['data_collection_date']+
            '_diffrn_detector.diffrn_id              1\n'
            '#\n'
            '_diffrn_radiation.diffrn_id                        1\n'
            '_diffrn_radiation.wavelength_id                    1\n'
            "_diffrn_radiation.pdbx_diffrn_protocol             'SINGLE WAVELENGTH'\n"
            '#\n'
            '_diffrn_radiation_wavelength.id           1\n'
            '_diffrn_radiation_wavelength.wavelength   %s\n'                                    %depositDict['radiation_wavelengths']+
            '#\n'
            '#\n'
            'loop_\n'
            '_entity.id\n'
            '_entity.type\n'
            '_entity.src_method\n'
            '_entity.pdbx_description\n'
            '_entity.pdbx_mutation\n'
            '1 polymer     man "%s" %s\n'                                                          %(depositDict['Source_organism_gene'],depositDict['fragment_name_one_specific_mutation'])+
            '#\n'
            'loop_\n'
            '_entity_poly.entity_id\n'
            '_entity_poly.type\n'
            '_entity_poly.pdbx_seq_one_letter_code\n'
            '_entity_poly.pdbx_strand_id\n'
            '_entity_poly.pdbx_seq_db_id\n'
            '_entity_poly.pdbx_seq_db_name\n'
            '1 "polypeptide(L)"\n'
            +molecule_one_letter_sequence+'\n'
#        ';MDPEVTLLLQCPGGGLPQEQIQAELSPAHDRRPLPGGDEAITAIWETRLKAQPWLFDAPK\n'
#        'FRLHSATLAPIGSRGPQLLLRLGLTSYRDFLGTNWSSSAAWLRQQGATDWGDTQAYLADP\n'
#        'LGVGAALATADDFLVFLRRSRQVAEAPGLVDVPGGHPEPQALCPGGSPQHQDLAGQLVVH\n'
#        'ELFSSVLQEICDEVNLPLLTLSQPLLLGIARNETSAGRASAEFYVQCSLTSEQVRKHYLS\n'
#        'GGPEAHESTGIFFVETQNVQRLLETEMWAELCPSAKGAIILYNRVQGSPTGAALGSPALL\n'
#        'PPL\n'
            ';\n'
            '%s %s UNP\n'                                        %(depositDict['protein_chains'],depositDict['molecule_one_letter_sequence_uniprot_id'])+
            '#\n'
            'loop_\n'
            '_entity_src_gen.entity_id\n'
            '_entity_src_gen.gene_src_strain\n'
            '_entity_src_gen.pdbx_gene_src_scientific_name\n'
            '_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id\n'
            '_entity_src_gen.pdbx_host_org_scientific_name\n'
            '_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id\n'
            '1 ? "%s" %s  "%s" %s\n'                %(depositDict['Source_organism_scientific_name'],pdbx_gene_src_ncbi_taxonomy_id,depositDict['Expression_system_scientific_name'],pdbx_host_org_ncbi_taxonomy_id)+
            '#\n'
#            '#\n'
##            '_pdbx_contact_author.id                  1\n'
##            "_pdbx_contact_author.address_1           '%s'\n"                           %depositDict['contact_author_PI_address']+
##            '_pdbx_contact_author.address_2           "%s"\n'                           %depositDict['contact_author_PI_organization_name']+
##            '_pdbx_contact_author.city                %s\n'                             %depositDict['contact_author_PI_city']+
##            "_pdbx_contact_author.state_province      '%s'\n"                           %depositDict['contact_author_PI_State_or_Province']+
##            '_pdbx_contact_author.postal_code         %s\n'                             %depositDict['contact_author_PI_Zip_Code']+
##            '_pdbx_contact_author.email               %s\n'                             %depositDict['contact_author_PI_email']+
##            '_pdbx_contact_author.name_first          %s\n'                             %depositDict['contact_author_PI_first_name']+
##            '_pdbx_contact_author.name_last           %s\n'                             %depositDict['contact_author_PI_last_name']+
##            '_pdbx_contact_author.country             "%s"\n'                           %depositDict['contact_author_PI_Country']+
##            '_pdbx_contact_author.phone               %s\n'                             %depositDict['contact_author_PI_phone_number']+
##            '_pdbx_contact_author.role                "%s"\n'                           %depositDict['contact_author_PI_role']+
##            '_pdbx_contact_author.organization_type   %s\n'                             %depositDict['contact_author_PI_organization_type']+
##            '#\n'
#            '_pdbx_contact_author.id                  2\n'
#            "_pdbx_contact_author.address_1           '%s'\n"                           %depositDict['contact_author_address']+
#            '_pdbx_contact_author.address_2           "%s"\n'                           %depositDict['contact_author_organization_name']+
#            '_pdbx_contact_author.city                %s\n'                             %depositDict['contact_author_city']+
#            "_pdbx_contact_author.state_province      '%s'\n"                           %depositDict['contact_author_State_or_Province']+
#            '_pdbx_contact_author.postal_code         %s\n'                             %depositDict['contact_author_Zip_Code'].replace(' ','')+
#            '_pdbx_contact_author.email               %s\n'                             %depositDict['contact_author_email']+
#            '_pdbx_contact_author.name_first          %s\n'                             %depositDict['contact_author_first_name']+
#            '_pdbx_contact_author.name_last           %s\n'                             %depositDict['contact_author_last_name']+
#            '_pdbx_contact_author.country             "%s"\n'                           %depositDict['contact_author_Country']+
#            '_pdbx_contact_author.phone               %s\n'                             %depositDict['contact_author_phone_number']+
#            '_pdbx_contact_author.role                "%s"\n'                           %depositDict['contact_author_role']+
#            '_pdbx_contact_author.organization_type   %s\n'                             %depositDict['contact_author_organization_type']+
            'loop_\n'
            '_pdbx_contact_author.id                  \n'
            "_pdbx_contact_author.address_1           \n"
            '_pdbx_contact_author.address_2           \n'
            '_pdbx_contact_author.city                \n'
            "_pdbx_contact_author.state_province      \n"
            '_pdbx_contact_author.postal_code         \n'
            '_pdbx_contact_author.email               \n'
            '_pdbx_contact_author.name_first          \n'
            '_pdbx_contact_author.name_last           \n'
            '_pdbx_contact_author.country             \n'
            '_pdbx_contact_author.phone               \n'
            '_pdbx_contact_author.role                \n'
            '_pdbx_contact_author.organization_type   \n'
            "1 '%s' '%s' '%s' '%s' '%s' %s %s '%s' '%s' '%s' '%s' %s\n" %(depositDict['contact_author_PI_address'],depositDict['contact_author_PI_organization_name'],depositDict['contact_author_PI_city'],depositDict['contact_author_PI_State_or_Province'],depositDict['contact_author_PI_Zip_Code'],depositDict['contact_author_PI_email'],depositDict['contact_author_PI_first_name'],depositDict['contact_author_PI_last_name'],depositDict['contact_author_PI_Country'],depositDict['contact_author_PI_phone_number'],depositDict['contact_author_PI_role'],depositDict['contact_author_PI_organization_type'])+
            "2 '%s' '%s' '%s' '%s' '%s' %s %s '%s' '%s' '%s' '%s' %s\n" %(depositDict['contact_author_address'],depositDict['contact_author_organization_name'],depositDict['contact_author_city'],depositDict['contact_author_State_or_Province'],depositDict['contact_author_Zip_Code'].replace(' ',''),depositDict['contact_author_email'],depositDict['contact_author_first_name'],depositDict['contact_author_last_name'],depositDict['contact_author_Country'],depositDict['contact_author_phone_number'],depositDict['contact_author_role'],depositDict['contact_author_organization_type'])+
            '#\n'
            'loop_\n'
            '_audit_author.name\n'
            "'%s, %s.'\n" %(depositDict['contact_author_last_name'],depositDict['contact_author_first_name'][0])+
#            "'Krojer, T.'\n"
#            "'Von Delft, F.'\n"
            '#\n'
            '_citation.id                        primary\n'
            "_citation.title                     '%s'\n"   %depositDict['group_title']+
            "_citation.journal_abbrev            'To Be Published'\n"
            '#\n'
            'loop_\n'
            '_citation_author.citation_id\n'
            '_citation_author.name\n'
            +primary_citation_author_name+
#        "primary 'Krojer, T.'\n"
#        "primary 'Von Delft, F.'\n"
            '#\n'
            '_struct.entry_id                     UNNAMED\n'
            '_struct.title\n'
            ';%s\n'                                                         %depositDict['title']+
            ';\n'
            '#\n'
            '_struct_keywords.entry_id        UNNAMED\n'
            '_struct_keywords.text            "%s"\n'                       %depositDict['structure_keywords']+
            '#\n'
            '_pdbx_struct_assembly_depositor_info.id                   1\n'
            "_pdbx_struct_assembly_depositor_info.method_details       PISA\n"
            '_pdbx_struct_assembly_depositor_info.oligomeric_count     %s\n'        %depositDict['biological_assembly_chain_number']+
            '#\n'
            '#\n'
            )

        return data_template_cif



class update_depositTable(QtCore.QThread):
    def __init__(self,deposit_dict,database,xce_logfile):
        QtCore.QThread.__init__(self)
        self.deposit_dict=deposit_dict
        self.database=database
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(database)

    def run(self):
        self.Logfile.insert('all entries in the depositTable will be updated with the following values:')
        for key in self.deposit_dict:
            self.Logfile.insert(key+': '+self.deposit_dict[key])
        dbEntries=self.db.execute_statement("select CrystalName,StructureType from depositTable;")
        for item in dbEntries:
            xtal=str(item[0])
            type=str(item[1])
            db_dict=self.deposit_dict   # need to do this because individual fields might need updating for some xtals

            # try to get information about the diffraction experiment
            diffractionExperiment=self.db.execute_statement("select DataCollectionBeamline,DataCollectionDate from mainTable where CrystalName is '%s'" %xtal)
            beamline=str(diffractionExperiment[0][0])
            date=str(diffractionExperiment[0][1])
            if beamline.lower() != 'none':
                db_dict=self.tweak_deposit_dict(xtal,db_dict)
            if date.lower() != 'none':
                db_dict['data_collection_date']=date.split()[0]

            self.Logfile.insert('updating depositTable for '+xtal+' @ '+type)
            self.db.update_depositTable(xtal,type,db_dict)
        self.Logfile.insert('Note: use DBbrowser to edit individual entries')

    def tweak_deposit_dict(self,xtal,db_dict):
        dls_beamlines=['i02','i03','i04','i04-1','i23','i24']
        dls_beamline_dict = {   'i02':      ['DIAMOND BEAMLINE I02',    'DECTRIS PILATUS 6M'],
                                'i03':      ['DIAMOND BEAMLINE I03',    'DECTRIS PILATUS 6M'],
                                'i04':      ['DIAMOND BEAMLINE I04',    'DECTRIS PILATUS 6M'],
                                'i04-1':    ['DIAMOND BEAMLINE I04-1',  'DECTRIS PILATUS 6M'],
                                'i23':      ['DIAMOND BEAMLINE I23',    'DECTRIS PILATUS 12M'],
                                'i24':      ['DIAMOND BEAMLINE I24',    'DECTRIS PILATUS 6M'] ,     }

        if db_dict['radiation_source_type'] in dls_beamlines:
            db_dict['radiation_source_type']=    dls_beamline_dict[db_dict['radiation_source_type']][0]
            db_dict['radiation_detector_type']=  dls_beamline_dict[db_dict['radiation_source_type']][1]
            db_dict['radiation_detector']=       'PIXEL'
            db_dict['radiation_source']=         'SYNCHROTRON'

        return db_dict




class prepare_mmcif_files_for_deposition(QtCore.QThread):

    def __init__(self,database,xce_logfile,overwrite_existing_mmcif,projectDir,structureType):
        QtCore.QThread.__init__(self)
        self.database=database
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(database)
        self.overwrite_existing_mmcif=overwrite_existing_mmcif
        self.projectDir=projectDir
        self.data_template_dict={}
        self.data_template_apo='data_template_apo.cif'
        self.data_template_bound='data_template_bound.cif'
        self.structureType=structureType
        self.depositLog=XChemLog.depositLog('deposit.log')
        self.depositLog.text('starting preparation of mmcif files for deposition')
        self.success=0
        self.successDict={}
        self.failureDict={}
        self.n_toDeposit=0
        self.counter=1


    def run(self):

        if self.structureType=='apo':
            # remove apo entries for structures to be deposited
            self.Logfile.insert('removing apo structures from depositTable which have a ligand bound structure ready for deposition')
            toDelete=self.db.execute_statement("select CrystalName from mainTable where RefinementOutcome like '5%';")
            toDeleteList=[]
            for xtal in toDelete:
                toDeleteList.append(str(xtal[0]))
            if toDeleteList != []:
                self.db.remove_selected_apo_structures_from_depositTable(self.xce_logfile,toDeleteList)
            self.Logfile.insert('checking for apo structures in depositTable')
            toDeposit=self.db.execute_statement("select CrystalName from depositTable where StructureType is 'apo';")

        elif self.structureType=='ligand_bound':
            self.Logfile.insert('checking for models in mainTable that are ready for deposition')
            toDeposit=self.db.execute_statement("select CrystalName from mainTable where RefinementOutcome like '5%';")

        self.n_toDeposit=len(toDeposit)
        self.depositLog.text('trying to prepare mmcif files for '+str(len(toDeposit))+' structures')

        progress_step=1
        if len(toDeposit) != 0:
            progress_step=100/float(len(toDeposit))
        else:
            progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)



        for item in sorted(toDeposit):
            xtal=str(item[0])
            eventMtz=[]
            preparation_can_go_ahead=True
            self.data_template_dict={}

            self.depositLog.modelInfo(xtal,self.structureType)

            if self.structureType=='apo':
                self.out=xtal+'-apo'
            elif self.structureType=='ligand_bound':
                self.out=xtal+'-bound'

            if self.structureType=='ligand_bound':
                self.Logfile.insert(xtal+' is ready for deposition')
                self.Logfile.insert('checking refinement stage of respective PanDDA sites...')
                panddaSites=self.db.execute_statement("select CrystalName,RefinementOutcome,PANDDA_site_event_map_mtz from panddaTable where CrystalName is '%s' and PANDDA_site_ligand_placed is 'True'" %xtal)
                self.Logfile.insert('found '+str(len(panddaSites))+' ligands')
                for site in panddaSites:
                    if str(site[1]).startswith('5'):
                        self.Logfile.insert('site is ready for deposition')
                        eventMtz.append(str(site[2]))
                    else:
                        self.Logfile.insert('site is NOT ready for deposition')
                        self.updateFailureDict(xtal,'at least one PanDDA site is not ready for deposition')
                        eventMtz.append(str(site[2]))
                        preparation_can_go_ahead=False

            n_eventMtz = len(eventMtz)
            if preparation_can_go_ahead:
                self.depositLog.nEvents(xtal,n_eventMtz)
                if self.structureType=='ligand_bound':
                    ModelData=self.db.execute_statement("select RefinementPDB_latest,RefinementMTZ_latest,RefinementCIF,DataProcessingPathToLogfile,RefinementProgram,CompoundCode,CompoundSMILES,RefinementMTZfree from mainTable where CrystalName is '%s'" %xtal)
                    pdb=str(ModelData[0][0])
                    mtzFinal=str(ModelData[0][1])
                    if not os.path.isfile(mtzFinal):
                        self.Logfile.insert('cannot find refine.mtz for bound structure of '+xtal+'; skipping => ERROR')
                        self.depositLog.text('cannot find refine.mtz for bound structure of '+xtal+'; skipping => ERROR')
                        self.updateFailureDict(xtal,'cannot find refine.mtz')
                        continue
                    mtzData=str(ModelData[0][7])
                    if not os.path.isfile(mtzData):
                        self.Logfile.insert('cannot find data.mtz for bound structure of '+xtal+'; skipping => ERROR')
                        self.depositLog.text('cannot find data.mtz for bound structure of '+xtal+'; skipping => ERROR')
                        self.updateFailureDict(xtal,'cannot find data.mtz')
                        continue

                if self.structureType=='apo':
                    ModelData=self.db.execute_statement("select DimplePathToPDB,DimplePathToMTZ,RefinementCIF,DataProcessingPathToLogfile,RefinementProgram,CompoundCode,CompoundSMILES,ProjectDirectory,RefinementMTZfree from mainTable where CrystalName is '%s'" %xtal)

                    if os.path.isfile(os.path.join(str(ModelData[0][7]),xtal,str(ModelData[0][0]))):
                        pdb=os.path.join(str(ModelData[0][7]),xtal,str(ModelData[0][0]))
                    elif os.path.isfile(str(ModelData[0][0])):
                        pdb=str(ModelData[0][0])
                    else:
                        self.Logfile.insert('cannot find PDB file for apo structure of '+xtal+'; skipping => ERROR')
                        self.depositLog.text('cannot find PDB file for apo structure of '+xtal+'; skipping => ERROR')
                        self.updateFailureDict(xtal,'cannot find PDB file')
                        continue

                    if os.path.isfile(os.path.join(str(ModelData[0][7]),xtal,str(ModelData[0][1]))):
                        mtzFinal=os.path.join(str(ModelData[0][7]),xtal,str(ModelData[0][1]))
                    elif os.path.isfile(str(ModelData[0][1])):
                        mtzFinal=str(ModelData[0][1])
                    else:
                        self.Logfile.insert('cannot find dimple.mtz for apo structure of '+xtal+'; skipping => ERROR')
                        self.depositLog.text('cannot find dimple.mtz for apo structure of '+xtal+'; skipping => ERROR')
                        self.updateFailureDict(xtal,'cannot find dimple.mtz file')
                        continue

                    if os.path.isfile(os.path.join(str(ModelData[0][7]),xtal,str(ModelData[0][8]))):
                        mtzData=os.path.join(str(ModelData[0][7]),xtal,str(ModelData[0][8]))
                    elif os.path.isfile(str(ModelData[0][8])):
                        mtzData=str(ModelData[0][8])
                    else:
                        self.Logfile.insert('cannot find data.mtz for apo structure of '+xtal+'; skipping => ERROR')
                        self.depositLog.text('cannot find data.mtz for apo structure of '+xtal+'; skipping => ERROR')
                        self.updateFailureDict(xtal,'cannot find data.mtz file')
                        continue

#                    if os.path.isfile(os.path.join(str(ModelData[0][7]),xtal,str(ModelData[0][1]))):
#                        mtz=os.path.join(str(ModelData[0][7]),xtal,str(ModelData[0][1]))


                cif=str(ModelData[0][2])
                if self.structureType=='ligand_bound':
                    if not os.path.isfile(cif):
                        if not os.path.isfile(os.path.join(self.projectDir,xtal,cif)):
                            self.Logfile.insert('cannot find ligand CIF file! Please check %s and the database!' %(os.path.join(self.projectDir,xtal)))
                            self.depositLog.text('cannot find ligand CIF file for %s; skipping... => ERROR' %xtal)
                            self.Logfile.insert('cannot prepare mmcif files for %s; skipping... => ERROR' %xtal)
                            self.updateFailureDict(xtal,'cannot find CIF file for ligand')
                            continue
                log=str(ModelData[0][3])
                refSoft=str(ModelData[0][4])
                compoundID=str(ModelData[0][5])
                smiles=str(ModelData[0][6])

                # if all modelled ligands are ready for deposition, we can continue

                # remove existing mmcif files and change DB accordingly
                self.remove_existing_mmcif_files(xtal)

                # first get all meta-data for deposition, i.e. data_template file
                wavelength=self.prepare_data_template_for_xtal(xtal,compoundID,pdb)

                # make model mmcif
                self.make_model_mmcif(xtal,pdb,log,refSoft)

                # make SF mmcif
                self.make_sf_mmcif(xtal,mtzFinal,mtzData,eventMtz)

                # report any problems that came up
                mmcif=self.check_mmcif_files_and_update_db(xtal,n_eventMtz,wavelength)

                # add ligand cif file to mmcif
                if self.structureType=='ligand_bound':
                    self.add_ligand_cif_file(mmcif,cif)

                self.Logfile.insert('finished preparation of mmcif files')


            else:
                self.Logfile.insert(XChemToolTips.deposition_pandda_site_not_ready(xtal))

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.summary()

    def remove_existing_mmcif_files(self,xtal):
        if self.overwrite_existing_mmcif:
            os.chdir(os.path.join(self.projectDir,xtal))
            for cif in glob.glob('data_template*'):
                self.Logfile.insert(xtal+' -> removing exisiting file: '+cif)
                os.system('/bin/rm '+cif)
            for mmcif in glob.glob('*.mmcif'):
                self.Logfile.insert(xtal+' -> removing exisiting file: '+mmcif)
                os.system('/bin/rm '+mmcif)
                self.db.execute_statement("update depositTable set mmCIF_model_file='',mmCIF_SF_file='' where CrystalName is '%s' and StructureType is '%s'" %(xtal,self.structureType))

    def updateFailureDict(self,xtal,error):
        if xtal not in self.failureDict:
            self.failureDict[xtal]=[]
        self.failureDict[xtal].append(error)


    def summary(self):
        self.depositLog.summary(self.n_toDeposit,self.success,self.failureDict,self.structureType,self.successDict)


    def prepare_data_template_for_xtal(self,xtal,compoundID,pdb):
        # check if file exists
        if self.overwrite_existing_mmcif:
            os.chdir(os.path.join(self.projectDir,xtal))
            data_template_dict=self.db.get_deposit_dict_for_sample(xtal)

            # edit title
            data_template_dict['group_title']=data_template_dict['group_deposition_title'].replace('$ProteinName',data_template_dict['Source_organism_gene']).replace('$CompoundName',compoundID)
            self.Logfile.insert('group deposition title for '+xtal+': '+data_template_dict['group_title'])
            self.depositLog.text('group title: '+data_template_dict['group_title'])
            if self.structureType=='ligand_bound':
                title=data_template_dict['structure_title'].replace('$ProteinName',data_template_dict['Source_organism_gene']).replace('$CompoundName',compoundID)
            if self.structureType=='apo':
                title=data_template_dict['structure_title_apo'].replace('$ProteinName',data_template_dict['Source_organism_gene']).replace('$CompoundName',compoundID).replace('$n',str(self.counter))
                self.counter+=1
            data_template_dict['group_description']=data_template_dict['group_description'].replace('$ProteinName',data_template_dict['Source_organism_gene'])
            data_template_dict['title']=data_template_dict['group_title']+' -- '+title
            self.Logfile.insert('deposition title for '+xtal+': '+data_template_dict['title'])
            self.depositLog.text('title: '+data_template_dict['title'])
            if ('$ProteinName' or '$CompoundName') in data_template_dict['title']:
                self.updateFailureDict(xtal,'title not correctly formatted')

            # mutations
            mutations=data_template_dict['fragment_name_one_specific_mutation']
            if mutations.lower().replace(' ','').replace('none','').replace('null','') == '':
                data_template_dict['fragment_name_one_specific_mutation']='?'
            else:
                data_template_dict['fragment_name_one_specific_mutation']='"'+mutations.replace(' ','')+'"'

            # get protein chains
            data_template_dict['protein_chains']=''
            chains=pdbtools(pdb).GetProteinChains()
            for item in chains:
                data_template_dict['protein_chains']+=item+','
            data_template_dict['protein_chains']=data_template_dict['protein_chains'][:-1]

            if self.structureType=='ligand_bound':
                self.Logfile.insert('creating %s file for ligand bound structure of %s' %(self.data_template_bound,xtal))
                data_template=templates().data_template_cif(data_template_dict)
                site_details=self.make_site_description(xtal)
                data_template+=site_details
                f=open(os.path.join(self.projectDir,xtal,self.data_template_bound),'w')
            elif self.structureType=='apo':
                self.Logfile.insert('creating %s file for apo structure of %s' %(self.data_template_apo,xtal))
                data_template=templates().data_template_cif(data_template_dict)
                site_details=self.make_site_description(xtal)
                data_template+=site_details
                f=open(os.path.join(self.projectDir,xtal,self.data_template_apo),'w')

            self.data_template_dict=data_template_dict
            f.write(data_template)
            f.close()

            wavelength=data_template_dict['radiation_wavelengths']
            return wavelength



    def make_model_mmcif(self,xtal,pdb,log,refSoft):
        if os.path.isfile(pdb) and os.path.isfile(log):
            os.chdir(os.path.join(self.projectDir,xtal))

            if self.structureType=='ligand_bound':
                out=xtal+'-bound'
                data_template=self.data_template_bound
                refSoft='PHENIX'
            elif self.structureType=='apo':
                out=xtal+'-apo'
                data_template=self.data_template_apo
                refSoft='REFMAC'

            refSoft=pdbtools(pdb).GetRefinementProgram()

            Cmd = ( 'source '+os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/setup.sh')+'\n'
                    +os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/bin/pdb_extract')+
#                    ' -r PHENIX'
                    ' -r %s'            %refSoft+
#                    ' -iLOG initial.log'
                    ' -iPDB %s'         %pdb+
#                    ' -i %s'            %self.data_template_dict['data_integration_software']+
#                    ' -p %s'            %self.data_template_dict['phasing_software']+
                    ' -e MR'
                    ' -s AIMLESS'
                    ' -iLOG %s'         %log+
                    ' -iENT %s'         %data_template+
                    ' -o %s.mmcif > %s.mmcif.log'      %(out,out)       )

            self.Logfile.insert('running pdb_extract: '+Cmd)
            os.system(Cmd)

            # can we add here the ligand.cif?
#    '''pdb_extract -r REFMAC -iLOG initial.log -iPDB initial.pdb -e MR -s AIMLESS -iLOG aimless.log -iENT data_template.cif -o NUDT22A-x0315-model.cif'''


    def make_sf_mmcif(self,xtal,mtzFinal,mtzData,eventMtz):
        if os.path.isfile(mtzFinal):
            os.chdir(os.path.join(self.projectDir,xtal))

            mtzin = mtzFinal+' '+mtzData+' '
            if self.structureType=='ligand_bound':
                out=xtal+'-bound'
                for event in eventMtz:
                    mtzin+=event+' '

            elif self.structureType=='apo':
                out=xtal+'-apo'

            Cmd = ( 'source '+os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/setup.sh')+'\n'
                    +os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/bin/sf_convert')+
                    ' -o mmcif'
                    ' -sf %s' %mtzin+
                    ' -out %s_sf.mmcif  > %s.sf_mmcif.log' %(out,out) )

#            Cmd = ( os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/sf-convert-v1.204-prod-src/bin/sf_convert')+
#                    ' -o mmcif'
#                    ' -sf %s' %mtzin+
#                    ' -out %s_sf.mmcif > %s.sf_mmcif.log' %(xtal,xtal) )

            self.Logfile.insert('running sf_convert: '+Cmd)
            os.system(Cmd)
            os.system('/bin/rm sf_format_guess.text mtzdmp.log SF_4_validate.cif sf_information.cif')


    def check_mmcif_files_and_update_db(self,xtal,n_eventMtz,wavelength):
        os.chdir(os.path.join(self.projectDir,xtal))
        foundFiles=True

        mmcif=''
        if os.path.isfile(self.out+'.mmcif') and os.path.getsize(self.out+'.mmcif') > 20000 :
            self.Logfile.insert('found '+self.out+'.mmcif')
            mmcif=os.path.join(self.projectDir,xtal,self.out+'.mmcif')
            cross_validation='?'
            starting_model='?'
            low_reso_list=[]
            softwareLine = 100000000
            foundSoftware=False
            softwareEntry=[]
            for n,line in enumerate(open(mmcif)):
                if foundSoftware:
                    if line.split()[0] == '#':
                        softwareLine=n
                        foundSoftware=False

                    try:
                        softwareEntry.append(int(line.split()[0]))
                    except ValueError:
                        pass

                elif '_refine.pdbx_ls_cross_valid_method' in line and len(line.split()) == 2:
                    cross_validation=line.split()[1]
                elif '_refine.pdbx_starting_model' in line and len(line.split()) == 2:
                    starting_model=line.split()[1]
                elif '_reflns.d_resolution_low' in line and len(line.split()) == 2:
                    low_reso_list.append(line.split()[1])
                elif '_refine.ls_d_res_low' in line and len(line.split()) == 2:
                    low_reso_list.append(line.split()[1])
                elif '_software.language' in line:
                    foundSoftware=True

            tmpText=''
            for n,line in enumerate(open(mmcif)):
                if '_refine.pdbx_ls_cross_valid_method' in line and cross_validation == '?':
                    tmpText+='_refine.pdbx_ls_cross_valid_method               THROUGHOUT \n'

                elif '_refine.pdbx_starting_model' in line and starting_model == '?':
                    tmpText+='_refine.pdbx_starting_model                      %s \n' %self.data_template_dict['pdbx_starting_model']

                elif '_refine.pdbx_method_to_determine_struct' in line:
                    tmpText+="_refine.pdbx_method_to_determine_struct          'FOURIER SYNTHESIS'\n"

                elif '_reflns.d_resolution_low' in line and len(line.split()) == 2:
                    tmpText+='_reflns.d_resolution_low             %s\n' %min(low_reso_list)

                elif '_refine.ls_d_res_low' in line and len(line.split()) == 2:
                    tmpText+='_refine.ls_d_res_low                             %s\n' %min(low_reso_list)

                elif n == softwareLine:
                    print 'software',softwareEntry
                    tmpText+=   (   "%s %s ? ? program ? ? 'data reduction' ? ?\n"      %(str(max(softwareEntry)+1),self.data_template_dict['data_integration_software'])+
                                    '%s %s ? ? program ? ? phasing ? ?\n'                %(str(max(softwareEntry)+2),self.data_template_dict['phasing_software'])+
                                    '#\n'   )

                else:
                    tmpText+=line

            f=open(mmcif,'w')
            f.write(tmpText)
            f.close()
            if os.path.isfile(mmcif):
                self.Logfile.insert('mmcif file successfully updated')

        else:
            self.Logfile.insert('cannot find '+self.out+'.mmcif; something went wrong! => ERROR')
            self.depositLog.text('cannot find '+self.out+'.mmcif; something went wrong! => ERROR')
            self.updateFailureDict(xtal,'cannot find '+self.out+'.mmcif')
            foundFiles=False

        if os.path.isfile(self.out+'_sf.mmcif') and os.path.getsize(self.out+'_sf.mmcif') > 20000 :
            self.Logfile.insert('found '+self.out+'_sf.mmcif')
            mmcif_sf=os.path.join(self.projectDir,xtal,self.out+'_sf.mmcif')

            # now check if really all event maps are in mmcif file
            if self.structureType=='ligand_bound':
                n_eventMTZ_found=-2     # set to -2 since first two data blocks are initial/final.mtz and data.mtz
                for line in open(mmcif_sf):
                    if line.startswith('_refln.crystal_id'):
                        n_eventMTZ_found+=1
                if n_eventMTZ_found != n_eventMtz:
                    self.Logfile.insert('%s event map mtz files were specified as input, but only %s ended up in the mmcif file => ERROR' %(str(n_eventMtz),str(n_eventMTZ_found)))
                    self.depositLog.text('%s event map mtz files were specified as input, but only %s ended up in the mmcif file => ERROR' %(str(n_eventMtz),str(n_eventMTZ_found)))
                    self.updateFailureDict(xtal,'%s event mtz in input; only %s in mmcif SF file' %(str(n_eventMtz),str(n_eventMTZ_found)))
                    foundFiles=False
                else:
                    self.Logfile.insert('%s event map mtz files were specified as input, %s ended up in the mmcif file, all well so far...' %(str(n_eventMtz),str(n_eventMTZ_found)))

            self.Logfile.insert('editing wavelength information in SF mmcif file; changing wavelength to %s' %wavelength)

            apo = [     "data from inital refinement with DIMPLE, initial.mtz",
                        "data from original reflection, data.mtz"]
            bound = [   "data from final refinement with ligand, final.mtz",
                        "data from original reflection, data.mtz",
                        "data for ligand evidence map (PanDDA event map), event_map1.mtz",
                        "data for ligand evidence map (PanDDA event map), event_map2.mtz",
                        "data for ligand evidence map (PanDDA event map), event_map3.mtz",
                        "data for ligand evidence map (PanDDA event map), event_map4.mtz",
                        "data for ligand evidence map (PanDDA event map), event_map5.mtz",
                        "data for ligand evidence map (PanDDA event map), event_map6.mtz",
                        "data for ligand evidence map (PanDDA event map), event_map7.mtz"       ]

            tmpText=''
            block=-1
            cif_block_found=0
            for line in open(mmcif_sf):
                if line.startswith('_cell.length_a'):
                    block+=1
                # need to do this because the data.mtz block could have missing screw axis
                if line.startswith('_symmetry.space_group_name_H-M') and cif_block_found == 1:
                    tmpText+=symmLine
                    cif_block_found+=1
                    continue
                if line.startswith('_symmetry.space_group_name_H-M') and cif_block_found == 0:
                    symmLine=line
                    cif_block_found+=1

                if line.startswith('_cell.angle_gamma'):
                    tmpText+=line
                    if self.structureType == 'apo':
                        addLines = (    '#\n'
                                        '_diffrn.id                  1\n'
                                        '_diffrn.details             "%s"\n' %apo[block]    )
                    if self.structureType == 'ligand_bound':
                        addLines = (    '#\n'
                                        '_diffrn.id                  1\n'
                                        '_diffrn.details             "%s"\n' %bound[block]    )
                    tmpText+=addLines
                    continue
                if line.startswith('_diffrn_radiation_wavelength.wavelength'):
                    tmpText+='_diffrn_radiation_wavelength.wavelength   %s\n' %wavelength
                else:
                    tmpText+=line
            f=open(mmcif_sf,'w')
            f.write(tmpText)
            f.close()
            if os.path.isfile(mmcif_sf):
                self.Logfile.insert('mmcif SF file successfully updated')
            else:
                self.Logfile.insert('something went wrong during the update => ERROR')
                self.depositLog.text('something went wrong during the update => ERROR')
                self.updateFailureDict(xtal,'error during mmcif SF update')
                foundFiles=False

        else:
            self.Logfile.insert('cannot find '+self.out+'_sf.mmcif; something went wrong! => ERROR')
            self.depositLog.text('cannot find '+self.out+'_sf.mmcif; something went wrong! => ERROR')
            self.updateFailureDict(xtal,'cannot find '+self.out+'_sf.mmcif')
            foundFiles=False

        if foundFiles:
            self.Logfile.insert('updating database with file locations for %s.mmcif and %s_sf.mmcif' %(self.out,self.out))
            self.successDict[xtal]=[self.out+'mmcif',os.path.getsize(self.out+'.mmcif'),self.out+'_sf.mmcif',os.path.getsize(self.out+'_sf.mmcif')]
            self.success+=1
            self.db.execute_statement("update depositTable set mmCIF_model_file='%s',mmCIF_SF_file='%s' where CrystalName is '%s' and StructureType is '%s'" %(mmcif,mmcif_sf,xtal,self.structureType))
        else:
            self.Logfile.insert('could not find %s.mmcif and/or %s_sf.mmcif; removing empty files...')
            os.system('/bin/rm %s.mmcif 2> /dev/null' %self.out)
            os.system('/bin/rm %s_sf.mmcif 2> /dev/null' %self.out)

        return mmcif


    def add_ligand_cif_file(self,mmcif,ligand_cif):
        self.Logfile.insert('adding ligand cif file to %s' %mmcif)
        tmpText=''
        for line in open(mmcif):
            tmpText+=line
        for line in open(ligand_cif):
            tmpText+=line
        f=open(mmcif,'w')
        f.write(tmpText)
        f.close()

    def make_site_description(self,xtal):
        mmcif_text='_pdbx_entry_details.nonpolymer_details\n;'

        general=self.db.execute_statement("select CompoundSMILES from mainTable where CrystalName is '%s'" %xtal)
        smiles=str(general[0][0])

        if self.structureType=='apo':
            if smiles.lower() != 'none' or smiles.lower() != "null":
                mmcif_text+='smiles string of soaked compound: %s;\n' %smiles
            else:
                mmcif_text=''

        if self.structureType=='ligand_bound':

            sqlite = (  "select "
                        " PANDDA_site_index,"
                        " PANDDA_site_x,"
                        " PANDDA_site_y,"
                        " PANDDA_site_y,"
                        " PANDDA_site_name,"
                        " PANDDA_site_confidence, "
                        " PANDDA_site_comment,"
                        " PANDDA_site_occupancy,"
                        " PANDDA_site_B_average,"
                        " PANDDA_site_B_ratio_residue_surroundings,"
                        " PANDDA_site_RSCC,"
                        " PANDDA_site_RSR,"
                        " PANDDA_site_RSZD,"
                        " PANDDA_site_rmsd "
                        "from panddaTable "
                        "where CrystalName is '%s' and PANDDA_site_ligand_placed is 'True' order by PANDDA_site_index ASC" %xtal)
#        panddaSites=self.db.execute_statement("select PANDDA_site_index,PANDDA_site_x,PANDDA_site_y,PANDDA_site_y,PANDDA_site_name,PANDDA_site_confidence from panddaTable where CrystalName is '%s' and PANDDA_site_ligand_placed is 'True' order by PANDDA_site_index ASC" %xtal)
            panddaSites=self.db.execute_statement(sqlite)

            root = etree.Element(xtal)
            child1 = etree.SubElement(root, "used_for_statistical_map")
            child1.text = 'yes'
            child2 = etree.SubElement(root, "smiles_of_compound_added")
            child2.text = '%s' %smiles

            site_descpription_complete=True
            for site in panddaSites:
                SiteIndex=  str(site[0]).replace(' ','')
                x_coord=    str(site[1])
                y_coord=    str(site[2])
                z_coord=    str(site[3])
                label=      str(site[4])
                confidence= str(site[5])
                comment=    str(site[6])
                occupancy=  str(site[7])
                Baverage=   str(site[8])
                Bratio=     str(site[9])
                RSCC=       str(site[10])
                RSR=        str(site[11])
                RSZD=       str(site[12])
                RMSD=       str(site[13])

                if 'none' in (SiteIndex.lower() or x_coord.lower() or y_coord.lower() or z_coord.lower() or confidence.lower() or occupancy.lower() or Bratio.lower() or Baverage.lower() or RSZD.lower() or RMSD.lower() or RSCC.lower() or RSR.lower()):
                    site_descpription_complete=False

                child = etree.SubElement(root, "site"+SiteIndex)
                childa = etree.SubElement(child, "label")
                childa.text = '%s' %label
                childx = etree.SubElement(child, "coordinate")
                childx.text = '%s %s %s' %(x_coord,y_coord,z_coord)
                childy = etree.SubElement(child, "smiles")
                childy.text = '%s' %smiles
                childz = etree.SubElement(child, "confidence")
                childz.text = '%s' %confidence
                childb = etree.SubElement(child, "comment")
                childb.text = '%s' %comment
                childc = etree.SubElement(child, "occupancy")
                childc.text = '%s' %occupancy
                childd = etree.SubElement(child, "B_average")
                childd.text = '%s' %Baverage
                childe = etree.SubElement(child, "B_ratio")
                childe.text = '%s' %Bratio
                childf = etree.SubElement(child, "RSCC")
                childf.text = '%s' %RSCC
                childg = etree.SubElement(child, "RSR")
                childg.text = '%s' %RSR
                childh = etree.SubElement(child, "RSZD")
                childh.text = '%s' %RSZD
                childi = etree.SubElement(child, "RMSD")
                childi.text = '%s' %RMSD

#        # pretty string
            s = etree.tostring(root, pretty_print=True)
            mmcif_text+=s+';\n'

            self.depositLog.site_xml(xtal,s)

            if not site_descpription_complete:
                self.updateFailureDict(xtal,'site description incomplete')

        return mmcif_text





class prepare_for_group_deposition_upload(QtCore.QThread):

    def __init__(self,database,xce_logfile,depositDir):
        QtCore.QThread.__init__(self)
        self.database=database
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(database)
        self.depositDir=depositDir



    def run(self):

        TextIndex=''
        os.chdir(self.depositDir)

        # ligand bound structures
        toDeposit=self.db.execute_statement("select CrystalName,mmCIF_model_file,mmCIF_SF_file from depositTable where StructureType is 'ligand_bound';")
        for item in sorted(toDeposit):
            xtal=str(item[0])
            mmcif=str(item[1])
            mmcif_sf=str(item[2])
            if os.path.isfile(mmcif) and os.path.isfile(mmcif_sf):
                self.Logfile.insert('copying %s to %s' %(mmcif,self.depositDir))
                os.system('/bin/cp %s .' %mmcif)
                self.Logfile.insert('copying %s to %s' %(mmcif_sf,self.depositDir))
                os.system('/bin/cp %s .' %mmcif_sf)
            else:
                self.Logfile.insert('cannot find apo mmcif file for '+xtal+' => ERROR')

            text = (    'label: %s-ligand_bound\n'                      %xtal+
                        'description: ligand_bound structure of %s\n'   %xtal+
                        'model: %s\n'                                   %mmcif[mmcif.rfind('/')+1:]+
                        'sf: %s\n\n'                                    %mmcif_sf[mmcif_sf.rfind('/')+1:]          )
            TextIndex+=text

        f = open('index.txt','w')
        f.write(TextIndex)
        f.close()

        self.Logfile.insert('preparing tar archive...')
        os.system('tar -cvf ligand_bound_structures.tar *bound* index.txt')
        self.Logfile.insert('bzipping archive...')
        os.system('bzip2 ligand_bound_structures.tar')
        self.Logfile.insert('removing all bound mmcif files and index.txt file from '+self.depositDir)
        os.system('/bin/rm -f *bound*mmcif index.txt')
        self.Logfile.insert('done!')



        # apo structures
        TextIndex=''
        toDeposit=self.db.execute_statement("select CrystalName,mmCIF_model_file,mmCIF_SF_file from depositTable where StructureType is 'apo';")
        for item in sorted(toDeposit):
            xtal=str(item[0])
            mmcif=str(item[1])
            mmcif_sf=str(item[2])
            if 'x513' in xtal: continue
            if 'x562' in xtal: continue
            if os.path.isfile(mmcif) and os.path.isfile(mmcif_sf):
                self.Logfile.insert('copying %s to %s' %(mmcif,self.depositDir))
                os.system('/bin/cp %s .' %mmcif)
                self.Logfile.insert('copying %s to %s' %(mmcif_sf,self.depositDir))
                os.system('/bin/cp %s .' %mmcif_sf)
            else:
                self.Logfile.insert('cannot find apo mmcif file for '+xtal+' => ERROR')

            text = (    'label: %s-apo\n'                     %xtal+
                        'description: apo structure of %s\n'  %xtal+
                        'model: %s\n'                         %mmcif[mmcif.rfind('/')+1:]+
                        'sf: %s\n\n'                          %mmcif_sf[mmcif_sf.rfind('/')+1:]          )
            TextIndex+=text

        f = open('index.txt','w')
        f.write(TextIndex)
        f.close()

        self.Logfile.insert('preparing tar archive...')
        os.system('tar -cvf apo_structures.tar *apo* index.txt')
        self.Logfile.insert('bzipping archive...')
        os.system('bzip2 apo_structures.tar')
        self.Logfile.insert('removing all apo mmcif files and index.txt file from '+self.depositDir)
        os.system('/bin/rm -f *apo*mmcif index.txt')
        self.Logfile.insert('done!')




class import_PDB_IDs(QtCore.QThread):
    def __init__(self,pdbCodes,database,xce_logfile):
        QtCore.QThread.__init__(self)
        self.pdbCodes=pdbCodes
        self.database=database
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(database)

    def run(self):

        for line in self.pdbCodes.split('\n'):
            if len(line.split('/')) == 2 and '-ligand_bound' in line:
                xtal=line[:line.rfind('-ligand_bound')].replace(' ','')
                pdbID=line.split('/')[1].replace(' ','')
                self.Logfile.insert('setting PDB ID for '+xtal+' to '+pdbID)
                sqlite="UPDATE mainTable SET Deposition_PDB_ID='%s',RefinementOutcome='6 - Deposited' where CrystalName is '%s';" %(pdbID,xtal)
                self.db.execute_statement(sqlite)


class compare_smiles_in_db_with_ligand_in_pdb(QtCore.QThread):
    def __init__(self,projectDir,database,xce_logfile):
        QtCore.QThread.__init__(self)
        self.projectDir=projectDir
        self.database=database
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(database)
        self.ErrorDict={}

    def update_ErrorDict(self,xtal,message):
        if xtal not in self.ErrorDict:
            self.ErrorDict[xtal]=[]
        self.ErrorDict[xtal].append(message)

    def run(self):

        os.chdir(self.projectDir)

        progress_step=1
        if len(glob.glob('*')) != 0:
            progress_step=100/float(len(glob.glob('*')))
        else:
            progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)


        for xtal in glob.glob('*'):
            if os.path.isfile(os.path.join(xtal,'refine.pdb')):
                smiles=self.db.execute_statement("select CompoundSmiles from mainTable where CrystalName is '%s'" %xtal)
                LigandSmiles=str(smiles[0][0])
                elementDict_smiles=smilestools(LigandSmiles).ElementDict()

                pdb=pdbtools(os.path.join(xtal,'refine.pdb'))
                ligandList=pdb.ligand_details_as_list()
                for ligand in ligandList:
                    resname     = ligand[0]
                    chainID     = ligand[1]
                    resseq      = ligand[2]
                    altLoc      = ligand[3]
                    elementDict_ligand=pdb.ElementDict(resname,chainID,resseq,altLoc)
                    for element in elementDict_ligand:
                        if elementDict_ligand[element] != elementDict_smiles[element]:
                            self.update_ErrorDict(xtal, '%s %s %s %s contains different number of atoms than smiles in DB' %(resname,chainID,resseq,altLoc))
                            break


            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)


#
#
#        # need a flag which connects apo structures and bound structures
#
#        ready_to_deposit=self.db.execute_statement("select CrystalName from mainTable where RefinementOutcome like '5%'")
#
#
#
#        print 'hallo'
##        progress_step=1
##        if len(self.sample_list) != 0:
##            progress_step=100/float(len(self.sample_list))
##        progress=0
##        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
##
##            progress += progress_step
##            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
##
#
#
#    '''pdb_extract -r REFMAC -iLOG initial.log -iPDB initial.pdb -e MR -s AIMLESS -iLOG aimless.log -iENT data_template.cif -o NUDT22A-x0315-model.cif'''
#
#    '''./sf_convert -o mmcif -sf initial.mtz data.mtz -out NUDT22A-x0315-sf.cif'''

#class update_deposition_table:
#
#    def __init__(self,database):
#        self.db=XChemDB.data_source(database)
#
#    def PanDDA_models_to_deposit(self):
#        panddaModels=self.db.execute_statement('select CrystalName,PANDDA_site_index,RefinementOutcome,PANDDA_site_ligand_placed,PANDDA_site_name from panddaTable')
#        toDeposit={}
#        mismatch={}
#
#        for item in panddaModels:
#            xtalID=str(item[0])
#            site_index=str(item[1])
#            RefinementStage=str(item[2])
#            ligand_placed=str(item[3])
#            siteName=str(item[4])
#            if ligand_placed.replace(' ','').lower()=='true':
#                if xtalID in toDeposit:
#                    if not RefinementStage.startswith('3'):
#                        if xtalID not in mismatch:
#                            mismatch[xtalID]=[]
#                        mismatch[xtalID].append([xtalID,site_index])
#                        continue
#
#                if RefinementStage.startswith('3'):
#                    if xtalID not in toDeposit:
#                        toDeposit[xtalID]=[]
#                    toDeposit[xtalID].append([xtalID,site_index,siteName,RefinementStage])
#
#        return toDeposit,mismatch

class templates_old:
    def data_template_text(self,depositDict):

        structure_author_name=''
        for name in depositDict['structure_author_name'].split(';'):
            structure_author_name+='<structure_author_name=  %s>\n' %name

        primary_citation_author_name=''
        for name in depositDict['primary_citation_author_name'].split(';'):
            primary_citation_author_name+='<primary_citation_author_name=  %s>\n' %name

        molecule_one_letter_sequence=''
        counter=0
        for aa in depositDict['molecule_one_letter_sequence']:
            if counter < 70:
                molecule_one_letter_sequence+=aa
            if counter == 70:
                molecule_one_letter_sequence+='\n'+aa
                counter = 0
            counter+=1


        data_template_text = (
            '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '\n'
            '    THE DATA_TEMPLATE.TEXT FILE	FOR X-RAY STRUCTURE DEPOSITION           \n'
            '\n'
            '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '\n'
            ' 		REMINDER AND GUIDELINES FOR USING THIS FILE\n'
            '\n'
            "  1. Only strings (values) included between the 'lesser than' and 'greater than'\n"
            '     signs (<.....>) will be parsed for evaluation.\n'
            '\n'
            '  2. All the input strings CAN NOT contain the three characters (", < , >).\n'
            '     Blank spaces or carriage returns within <..> will be ignored.\n'
            '\n'
            '  3. NEVER change the data item names (first column) inside of the brackets.\n'
            '\n'
            '  4. NEVER remove the equal sign (=) after the data item in the brackets.\n'
            '\n'
            '  5. If more groups are needed, same number of data item should be added.\n'
            '\n'
            "  6. The items marked by '!' are mandatory.\n"
            '\n'
            '  7. Always check the log file (pdb_extract.log) which contains errors/warnings\n'
            '     when parsing the template file.\n'
            '\n'
            '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '++++                        START INPUT DATA BELOW                      ++++\n'
            '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '\n'
            '================CATEGORY 1:   Contact Authors=============================\n'
            'Enter information about the contact authors.(mandatory)\n'
            '\n'
            '1.  Information about the Principal investigator (PI).\n'
            '\n'
            '<contact_author_PI_id= 1 >           !(must be given 1)\n'
            '<contact_author_PI_salutation=  %s>    !(Dr./Prof./Mr./Mrs./Ms.)\n'                            %depositDict['contact_author_PI_salutation']+
            '<contact_author_PI_first_name=  %s>    !(e.g. John)\n'                                         %depositDict['contact_author_PI_first_name']+
            '<contact_author_PI_last_name=  %s>     !(e.g. Rodgers)\n'                                      %depositDict['contact_author_PI_last_name']+
            '<contact_author_PI_middle_name=  %s>\n'                                                        %depositDict['contact_author_PI_middle_name']+
            '<contact_author_PI_role= %s>  !(or responsible scientist)\n'                                   %depositDict['contact_author_PI_role']+
            '<contact_author_PI_organization_type= %s>  !(or commercial, government, other)\n'              %depositDict['contact_author_PI_organization_type']+
            '<contact_author_PI_email=  %s>              !(e.g.   name@host.domain.country)\n'              %depositDict['contact_author_PI_email']+
            '<contact_author_PI_address=  %s>            !(e.g. 610 Taylor road)\n'                         %depositDict['contact_author_PI_address']+
            '<contact_author_PI_city=  %s>               !(e.g. Piscataway)\n'                              %depositDict['contact_author_PI_city']+
            '<contact_author_PI_State_or_Province=  %s>  !(e.g.  New Jersey)\n'                             %depositDict['contact_author_PI_State_or_Province']+
            '<contact_author_PI_Zip_Code=  %s>           !(e.g.  08864)\n'                                  %depositDict['contact_author_PI_Zip_Code']+
            '<contact_author_PI_Country=  %s>        !(e.g. United States, United Kindom, . )\n'            %depositDict['contact_author_PI_Country']+
            '<contact_author_PI_fax_number=  %s>\n'                                                         %depositDict['contact_author_PI_fax_number']+
            '<contact_author_PI_phone_number=  %s>    !(e.g.  01(617) 555-1213 )\n'                         %depositDict['contact_author_PI_phone_number']+
            '\n'
            '2. Information about other contact authors (responsible scientist, investigator)\n'
            '\n'
            '<contact_author_id=  2>              (e.g. 2)\n'
            '<contact_author_salutation=  %s>\n'                                                            %depositDict['contact_author_salutation']+
            '<contact_author_first_name=  %s>\n'                                                            %depositDict['contact_author_first_name']+
            '<contact_author_last_name=  %s>\n'                                                             %depositDict['contact_author_last_name']+
            '<contact_author_middle_name=  %s>\n'                                                           %depositDict['contact_author_middle_name']+
            '<contact_author_role=  %s>\n'                                                                  %depositDict['contact_author_role']+
            '<contact_author_organization_type=  %s>\n'                                                     %depositDict['contact_author_organization_type']+
            '<contact_author_email=  %s>\n'                                                                 %depositDict['contact_author_email']+
            '<contact_author_address=  %s>\n'                                                               %depositDict['contact_author_address']+
            '<contact_author_city=  %s>\n'                                                                  %depositDict['contact_author_city']+
            '<contact_author_State_or_Province=  %s>\n'                                                     %depositDict['contact_author_State_or_Province']+
            '<contact_author_Zip_Code=  %s>\n'                                                              %depositDict['contact_author_Zip_Code']+
            '<contact_author_Country=  %s>\n'                                                               %depositDict['contact_author_Country']+
            '<contact_author_fax_number=  %s>\n'                                                            %depositDict['contact_author_fax_number']+
            '<contact_author_phone_number=  %s>\n'                                                          %depositDict['contact_author_phone_number']+
            '\n'
            '...(add more groups if needed)...\n'
            '\n'
            '================CATEGORY 2:   Release Status==============================\n'
            'Enter release status for the coordinates, structure_factor, and sequence\n'
            '\n'
            '  Status must be chosen from one of the following:\n'
            '* for coordinate & structure_factor\n'
            '  (RELEASE NOW, HOLD FOR PUBLICATION,  HOLD FOR 8 WEEKS,\n'
            '  HOLD FOR 6 MONTHS, HOLD FOR 1 YEAR)\n'
            '\n'
            '* for chemical sequence, give  RELEASE NOW  or  HOLD FOR RELEASE\n'
            '\n'
            '<Release_status_for_coordinates=  %s>      !(e.g. HOLD FOR PUBLICATION)\n'                     %depositDict['Release_status_for_coordinates']+
            '<Release_status_for_structure_factor=  %s> !(e.g. HOLD FOR PUBLICATION)\n'                     %depositDict['Release_status_for_structure_factor']+
            '<Release_status_for_sequence=  %s>         !(RELEASE NOW or  HOLD FOR RELEASE)\n'              %depositDict['Release_status_for_sequence']+
            '\n'
            '================CATEGORY 3:   Title=======================================\n'
            'Enter the title for the structure\n'
            '\n'
            '<structure_title=  %s>     !(e.g. Crystal Structure Analysis of the B-DNA)\n'                  %depositDict['structure_title']+
            '<structure_details=  %s>\n'                                                                    %depositDict['structure_details']+
            '\n'
            '================CATEGORY 4: Authors of Structure============================\n'
            'Enter authors of the deposited structures (at least one author)\n'
            '\n'
            +structure_author_name+
#        '<structure_author_name=  >  !(e.g. Surname, F.M.)\n'
#        '<structure_author_name=  >\n'
#        '<structure_author_name=  >\n'
#        '<structure_author_name=  >\n'
#        '<structure_author_name=  >\n'
#        '<structure_author_name=  >\n'
#        '<structure_author_name=  >\n'
#        '<structure_author_name=  >\n'
#        '<structure_author_name=  >\n'
            '\n'
            '...add more name if needed...\n'
            '\n'
            '================CATEGORY 5a:  Primary  Citation ============================\n'
            '\n'
            '  The primary citation is the article in which the deposited coordinates\n'
            '  were first reported. Other related citations may also be provided.\n'
            '\n'
            "  If the citation has not yet been published, give 'To be published' to the item\n"
            "  'primary_citation_journal_abbrev' and leave pages, year, volume blank.\n"
            '\n'
            'Enter the author name of primary citation\n'
            +primary_citation_author_name+
#        '<primary_citation_author_name=  >    !(e.g. Surname, F.M.)\n'
#        '<primary_citation_author_name=  >\n'
#        '<primary_citation_author_name=  >\n'
#        '<primary_citation_author_name=  >\n'
#        '<primary_citation_author_name=  >\n'
#        '<primary_citation_author_name=  >\n'
#        '<primary_citation_author_name=  >\n'
#        '<primary_citation_author_name=  >\n'
            '\n'
            '...add more name if needed...\n'
            '\n'
            'Enter journal information of the primary citation\n'
            '<primary_citation_id= %s>\n'                                                       %depositDict['primary_citation_id']+
            '<primary_citation_journal_abbrev=  %s>     (e.g. To be published)\n'               %depositDict['primary_citation_journal_abbrev']+
            '<primary_citation_title= %s>\n'                                                    %depositDict['primary_citation_title']+
            '<primary_citation_year=  %s>\n'                                                    %depositDict['primary_citation_year']+
            '<primary_citation_journal_volume=  %s>\n'                                          %depositDict['primary_citation_journal_volume']+
            '<primary_citation_page_first=  %s>\n'                                              %depositDict['primary_citation_page_first']+
            '<primary_citation_page_last=  %s>\n'                                               %depositDict['primary_citation_page_last']+
            '\n'
            '================CATEGORY 5b:  other citations (if applicable) ================\n'
            '\n'
            '1.Enter the author name of other citations\n'
            '<citation_author_id=  >    (e.g. 1, 2 ..)\n'
            '<citation_author_name=  >\n'
            '<citation_author_name=  >\n'
            '<citation_author_name=  >\n'
            '<citation_author_name=  >\n'
            '<citation_author_name=  >\n'
            '<citation_author_name=  >\n'
            '<citation_author_name=  >\n'
            '<citation_author_name=  >\n'
            '\n'
            '...add more name if needed...\n'
            '\n'
            '1. Enter journal information of the primary citation\n'
            '<citation_id= 1 >               (e.g. 1, 2, 3 ...)\n'
            '<citation_journal_abbrev=  >\n'
            '<citation_title=  >\n'
            '<citation_year=  >\n'
            '<citation_journal_volume=  >\n'
            '<citation_page_first=  >\n'
            '<citation_page_last=  >\n'
            '\n'
            '...(add more other citations if needed)...\n'
            '\n'
            '================CATEGORY 6:   Molecule Information========================\n'
            'Enter the names of the molecules (entities) that are in the asymmetric unit\n'
            '\n'
            'NOTE: The name of molecule should be obtained from the appropriate\n'
            '      sequence database reference, if available. Otherwise the gene name or\n'
            '      other common name of the entity may be used.\n'
            '      e.g. HIV-1 integrase for protein , RNA Hammerhead Ribozyme\n'
            '\n'
            '1. For entity 1\n'
            '<molecule_id= 1 >        (e.g. 1 )\n'
            '<molecule_name=  %s>       (e.g.  RNA Hammerhead Ribozyme )\n'                     %depositDict['molecule_name']+
            '<molecule_type= polymer >    (e.g. polymer , non-polymer, macrolide  )\n'
            '<molecule_source_method= man>   (e.g. man , nat, syn)\n'
            '\n'
            '2. For entity 2\n'
            '<molecule_id=  >     (e.g. 2 )\n'
            '<molecule_name=  >\n'
            '<molecule_type=  >\n'
            '<molecule_source_method= >\n'
            '\n'
            '...(add more group if needed)...\n'
            '\n'
            '================CATEGORY 7:   Molecule Details============================\n'
            'Enter additional information about each entity, if known. (optional)\n'
            '\n'
            '1. For entity 1\n'
            '<Molecular_entity_id= 1 >       (e.g. 1 )\n'
            '<Fragment_name=  %s>             (e.g. ligand binding domain, hairpin)\n'          %depositDict['Fragment_name']+
            '<Specific_mutation=  %s>         (e.g. C280S)\n'                                   %depositDict['Specific_mutation']+
            '<Enzyme_Comission_number=  %s>   (if known: e.g. 2.7.7.7)\n'                       %depositDict['Enzyme_Comission_number']+
            '\n'
            '2. For entity 2\n'
            '<Molecular_entity_id=  >       (e.g.  2 )\n'
            '<Fragment_name=  >\n'
            '<Specific_mutation=  >\n'
            '<Enzyme_Comission_number=  >\n'
            '\n'
            '...(add more group if needed)...\n'
            '\n'
            '================CATEGORY 8:   Genetically Manipulated Source=============\n'
            'Enter data in the genetically manipulated source category\n'
            '\n'
            '   If the biomolecule has been genetically manipulated, describe its\n'
            '   source and expression system here.\n'
            '\n'
            '1. For entity 1\n'
            '<Manipulated_entity_id= 1 >               !(e.g. 1 )\n'
            '<Source_organism_scientific_name=  %s>      !(e.g. Homo sapiens)\n'                %depositDict['Source_organism_scientific_name']+
            '<Source_organism_gene=  %s>                 (e.g. RPOD, ALKA...)\n'                %depositDict['Source_organism_gene']+
            '<Source_organism_strain=  %s>               (e.g. BH10 ISOLATE, K-12...)\n'        %depositDict['Source_organism_strain']+
            '<Expression_system_scientific_name=  %s>    (e.g. Escherichia coli)\n'             %depositDict['Expression_system_scientific_name']+
            '<Expression_system_strain= %s >	          (e.g. BL21(DE3))\n'                   %depositDict['Expression_system_strain']+
            '<Expression_system_vector_type=  %s>	  (e.g. plasmid)\n'                         %depositDict['Expression_system_vector_type']+
            '<Expression_system_plasmid_name=  %s>       (e.g. pET26)\n'                        %depositDict['Expression_system_plasmid_name']+
            '<Manipulated_source_details=  %s>           (any other relevant information)\n'    %depositDict['Manipulated_source_details']+
            '\n'
            '2. For entity 2\n'
            '\n'
            '...(add more group if needed)...\n'
            '\n'
            '================CATEGORY 9:   Natural Source (optional) ===================\n'
            'Enter data in the natural source category  (if applicable)\n'
            '\n'
            '   If the biomolecule was derived from a natural source, describe it here.\n'
            '\n'
            '1. For entity 1\n'
            '<natural_source_entity_id=  >          (e.g. 1, 2..)\n'
            '<natural_source_scientific_name=  >    (e.g. Homo sapiens)\n'
            '<natural_source_organism_strain=  >    (e.g. DH5a , BMH 71-18)\n'
            '<natural_source_details=  >            (e.g. organ, tissue, cell ..)\n'
            '\n'
            '2. For entity 2\n'
            '\n'
            '...(add more group if needed)...\n'
            '\n'
            '================CATEGORY 10:  Synthetic Source (optional)==================\n'
            'If the biomolecule has not been genetically manipulated or synthesized,\n'
            'describe its source here.\n'
            '\n'
            '1. For entity 1\n'
            '<synthetic_source_entity_id=  >          (e.g. 1,2. )\n'
            '<synthetic_source_description=  >      (if known)\n'
            '\n'
            '\n'
            '2. For entity 2\n'
            '\n'
            '...(add more group if needed)...\n'
            '\n'
            '\n'
            '================CATEGORY 11:   Keywords===================================\n'
            'Enter a list of keywords that describe important features of the deposited\n'
            'structure.\n'
            '\n'
            '   Example: beta barrel, protein-DNA complex, double helix, hydrolase, etc.\n'
            '\n'
            '<structure_keywords=  %s>  !(e.g. beta barrel)\n'                                      %depositDict['structure_keywords']+
            '\n'
            '================CATEGORY 12:   Biological Assembly ======================\n'
            'Enter data in the biological assembly category (if applicable)\n'
            '\n'
            'Enter the number of polymer chains that form the assembly in solution\n'
            '\n'
            '<biological_assembly_chain_number=  %s>  !(e.g.  1 for monomer, 2 for dimer ..)\n'     %depositDict['biological_assembly_chain_number']+
            '\n'
            '================CATEGORY 13:   Methods and Conditions=====================\n'
            'Enter the crystallization conditions for each crystal\n'
            '\n'
            '1. For crystal 1:\n'
            '<crystal_number= 1 >	            (e.g. 1, )\n'
            '<crystallization_method=  %s>      (e.g. BATCH MODE, EVAPORATION)\n'                   %depositDict['crystallization_method']+
            '<crystallization_pH=  %s>          (e.g. 7.5 ...)\n'                                   %depositDict['crystallization_pH']+
            '<crystallization_temperature=  %s> (e.g. 298) (in Kelvin)\n'                           %depositDict['crystallization_temperature']+
            '<crystallization_details=  %s>     (e.g. PEG 4000, NaCl etc.)\n'                       %depositDict['crystallization_details']+
            '\n'
            '...(add more crystal groups if needed)...\n'
            '\n'
            '================CATEGORY 14:   Radiation Source (experiment)============\n'
            'Enter the details of the source of radiation, the X-ray generator,\n'
            'and the wavelength for each diffraction.\n'
            '\n'
            '1. For experiment 1:\n'
            '<radiation_experiment= 1 >      !(e.g. 1, 2, ...)\n'
            '<radiation_source=  %s>           !(e.g. SYNCHROTRON, ROTATING ANODE ..)\n'            %depositDict['radiation_source']+
            '<radiation_source_type=  %s>      !(e.g. NSLS BEAMLINE X8C ..)\n'                      %depositDict['radiation_source_type']+
            '<radiation_wavelengths=  %s>       !(e.g. 1.502, or a list 0.987,0.988 ..)\n'          %depositDict['radiation_wavelengths']+
            '<radiation_detector=  %s>         !(e.g. CCD, AREA DETECTOR, IMAGE PLATE ..)\n'        %depositDict['radiation_detector']+
            '<radiation_detector_type=  %s>     !(e.g. ADSC QUANTUM 1,  ..)\n'                      %depositDict['radiation_detector_type']+
            '<radiation_detector_details=  >    (e.g. mirrors...)\n'
            '<data_collection_date=  %s>             !(e.g. 2004-01-07)\n'                          %depositDict['data_collection_date']+
            '<data_collection_temperature=  %s>      !(e.g. 100 for crystal  1:)\n'                 %depositDict['data_collection_temperature']+
            '<data_collection_protocol=  %s>          !(e.g. SINGLE WAVELENGTH, MAD, ...)\n'        %depositDict['data_collection_protocol']+
            '<data_collection_monochromator=  >     (e.g. GRAPHITE, Ni FILTER ...)\n'
            '<data_collection_monochromatic_or_laue=  M >  !(default M, give L if Laue diffr.)\n'
            '\n'
            '\n'
            '....(add more experiment group if needed)....\n'
            '\n'
            '================CATEGORY 15:   refinement details (optional)============\n'
            'Enter the details of the structure refinement. (if applicable)\n'
            '\n'
            '<refinement_detail=   >\n'
            '<refinement_start_model=   >    (e.g. pdbid 100D)\n'
            '\n'
            '================CATEGORY 16:   database (optional)======================\n'
            'Enter the database name for each molecule (entity), (IF KNOWN).\n'
            '\n'
            '1. For entity 1\n'
            '<database_entity_id= 1  >  (e.g. 1 )\n'
            '<database_name=  >  (e.g. BMCD, BMRB, EMDB, PDB, NDB, TargetTrack )\n'
            '<database_code=  >  (e.g. 1ABC, 100D, TNKS2_HUMAN )\n'
            '<database_accession=   >  (e.g. 100D, Q9H2K2  )\n'
            '\n'
            '\n'
            '...(add more group if needed)...\n'
            '\n'
            '================CATEGORY 17:   Ligand binding (optional)==================\n'
            '\n'
            '<binding_assay_id=   >     !(A unique identifier such as 1,2,..)\n'
            '<binding_assay_target_sequence_one_letter_code=   >  (Chemical sequence if known)\n'
            '<binding_assay_ligand_descriptor_type=   >  (e.g. SMILES, SMILES_CANONICAL, InChI,InChIKey)\n'
            '<binding_assay_ligand_descriptor=   >   (e.g. Cc1cccc(c1)C1COc2cc(O)c(O)cc2C1 )\n'
            "<binding_assay_assay_type=   >        (Type of binding assay. e.g. 'competitive binding')\n"
            '<binding_assay_assay_value_type=   > (e.g. IC50, EC50,Ki,Kd)\n'
            '<binding_assay_assay_value=   >     (The value measured. e.g. 8300.0)\n'
            '<binding_assay_assay_pH=   >       (pH value at which the assay was performed. e.g. 6.4)\n'
            '<binding_assay_assay_temperature=   > (temperature (K) at which the assay was performed. e.g.273)\n'
            '<binding_assay_details=   >  (details of the measurement).\n'
            '\n'
            '\n'
            '\n'
            '================CATEGORY 18:   Structure Genomic (optional)==============\n'
            "If it is the structure genome's project, give the information\n"
            '\n'
            '<SG_project_id=  1>\n'
            '<SG_project_name=  %s>        (e.g. PSI, Protein Structure Initiative)\n'              %depositDict['SG_project_name']+
            '<full_name_of_SG_center=  %s>   (e.g. Berkeley Structural Genomic Center)\n'           %depositDict['full_name_of_SG_center']+
            '\n'
            '\n'
            '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '+++ Categories below were extracted from the coordinate. Please check. +++\n'
            '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '\n'
            '================CATEGORY 19:   Crystallographic Data=======================\n'
            'Enter crystallographic data\n'
            '\n'
            '<space_group = P 43 21 2> (use International Table conventions)\n'
            '<unit_cell_a     =   71.467  >\n'
            '<unit_cell_b     =   71.467  >\n'
            '<unit_cell_c     =  150.264  >\n'
            '<unit_cell_alpha =  90.00  >\n'
            '<unit_cell_beta  =  90.00  >\n'
            '<unit_cell_gamma =  90.00  >\n'
            '\n'
            '\n'
            '================CATEGORY 20:   Sequence Information =======================\n'
            'Enter one letter sequence for each polymeric entity in asymmetric unit\n'
            '\n'
            '--------------------------------------------------------------------------\n'
            '			  SOME DEFINITIONS\n'
            '     An ENTITY is defined as any unique molecule present in the asymmetric\n'
            '     unit. Each unique biological polymer (protein or nucleic acids) in the\n'
            '     structure is considered an entity. Thus, if there are five copies of\n'
            '     a single protein in the asymmetric unit, the molecular entity is still\n'
            '     only one. Water and non-polymers like ions, ligands and sugars are\n'
            '     also entities.\n'
            '\n'
            '     Here we only consider the sequences of polymeric entities (protein or\n'
            '     nucleic acid).\n'
            '\n'
            '	        GUIDELINES FOR COMPLETING THIS CATEGORY\n'
            '\n'
            '      * The unique chemical sequences are extract from the coordinate file.\n'
            '      Never give a TER cards in the pdb file inside of a complete entity.\n'
            '\n'
            '      * In a PDB or mmCIF format file, all residues of a single polymeric\n'
            '      entity should have one chain ID. Multiple copies of the same entity\n'
            '      should each be assigned a unique chain ID. The multiple chain IDs\n'
            "      should be separated by commas as 'A,B,C,...'.\n"
            '\n'
            "      * Four question marks '????' are used to denote a polymer chain breaks,\n"
            '      which may be caused by missing residues due to unobsered electron density\n'
            '      or linkage problem due to poor geometry. Replace these question marks\n'
            '      with the sequence of residues missing from the coordinates. Delete these\n'
            '      question marks if the break is due to the poor geometry. Also add any\n'
            '      residues missing from the N- and/or C-termini here.\n'
            '\n'
            '      * If there are non-standard residues in the coordinates, this program\n'
            '      lists them according to the three letter code used in the coordinate\n'
            '      file as (ABC).\n'
            '\n'
            '      * If any residue was modeled as Ala or Gly due to lack of the side-chain\n'
            '      density, the sequence extracted here will represent them as A or G\n'
            '      respectively. Correct this to the original sequence that was present in\n'
            '      the crystal.\n'
            '\n'
            '\n'
            '----------------------------------------------------------------------------\n'
            '\n'
            '   Below is the one letter chemical sequence extracted from your PDB coordinate\n'
            '   file. The molecular entities are grouped and listed together.\n'
            '\n'
            '   PLEASE CHECK THE SEQUENCE of each entity carefully and modify it, as necessary.\n'
            '   Make sure that you REVIEW THE FOLLOWING:\n'
            '    * chain breaks due to missing residues,\n'
            '    * missing residues in the N- and/or C-termini,\n'
            '    * non-standard residues and\n'
            '    * cases of residues modeled as Ala or Gly due to missing side-chain density.\n'
            '\n'
            '<molecule_entity_id=1 >\n'
            '<molecule_entity_type=polypeptide(L) >\n'
            '<molecule_one_letter_sequence=\n'
            +molecule_one_letter_sequence+'\n'
#        'AQNPNCNIMIFHPTKEEFNDFDKYIAYMESQGAHRAGLAKIIPPKEWKARETYDNISEILIATPLQQVAS\n'
#        'GRAGVFTQYHKKKKAMTVGEYRHLANSKKYQTPPHQNFEDLERKYWKNRIYNSPIYGADISGSLFDENTK\n'
#        'QWNLGHLGTIQDLLEKECGVVIEGVNTPYLYFGMWKTTFAWHTEDMDLYSINYLHLGEPKTWYVVPPEHG\n'
#        'QRLERLARELFPGSSRGCGAFLRHKVALISPTVLKENGIPFNRITQEAGEFMVTFPYGYHAGFNHGFNCA\n'
#        'EAINFATPRWIDYGKMASQCSCGEARVTFSMDAFVRILQPERYDLWKRGQD >\n'
            '< molecule_chain_id=A >\n'
            '< target_DB_id=  > (if known)\n'
            '< sequence_database_id=  > (if known)\n'
            '< sequence_database_name=  > (if known)\n'
            '\n'
            '<molecule_entity_id=  >\n'
            '<molecule_entity_type=  >\n'
            '<molecule_one_letter_sequence= >\n'
            '<molecule_chain_id=  >\n'
            '<target_DB_id=  >  (if known)\n'
            '<sequence_database_id=   > (if known)\n'
            '<sequence_database_name=   > (if known)\n'
            '\n'
            '\n'
            '\n'
            '\n'
            '=====================================END==================================\n'
        )

        return data_template_text



#def find_apo_structures(panddaDir,projectDir,database):
#
#    # first check if structure is already present in DB and if so if all the
#    # information concur
#
#    # need to update pandda directory for every exported structure so that
#    # we know where to look for the pandda.log file that contains the relevant information
#
#    # update CrystalName_of_pandda_input in DB
#
#    # in DB: update StructureType field accordingly
#
#    # newer pandda versions seem to have severl copies of pandda.log with names like
#    # pandda-2016-09-01-2139.log
#    panddaLog=glob.glob(os.path.join(panddaDir,'pandda*log'))
#    panddaLog.sort(key=os.path.getmtime)
#
#    panddaVersion='unknown'
#    readindApoStructures = False
#    apoStructures = []
#    apoStructureDict = {}
#    for files in panddaLog:
#        for line in open(files):
#            if line.startswith('-  Pandda Version'):
#                if len(line.split()) >= 4:
#                    panddaVersion=line.split()[3]
#            if 'No Statistical Maps Found:' in line:
#                readindApoStructures=True
#            if readindApoStructures:
#                if 'Pickling Object: processed_datasets' in line:
#                    if line.split() >= 2:
#                        # e.g. line.split() -> ['Pickling', 'Object:', 'processed_datasets/NUDT22A-x0055/pickles/dataset.pickle']
#                        xtal=line.split()[2].split('/')[1]
#                        if os.path.isfile(os.path.join(panddaDir,'processed_datasets',xtal,xtal+'-pandda-input.pdb')):
#                            apoStructures.append(xtal)
#            if 'Pre-existing statistical maps (from previous runs) have been found and will be reused:' in line:
#                readindApoStructures=False
#        apoStructureDict[panddaDir]=apoStructures
#
#    return apoStructureDict


class prepare_bound_models_for_deposition(QtCore.QThread):

    def __init__(self,database,xce_logfile,overwrite_existing_mmcif,projectDir):
        QtCore.QThread.__init__(self)
        self.database=database
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(database)
        self.overwrite_existing_mmcif=overwrite_existing_mmcif
        self.projectDir=projectDir
        self.data_template='data_template.cif'

    def run(self):
        self.Logfile.insert('checking for models in mainTable that are ready for deposition')
        toDeposit=self.db.execute_statement("select CrystalName from mainTable where RefinementOutcome like '5%';")
        for item in toDeposit:
            xtal=str(item[0])
            self.Logfile.insert(xtal+' is ready for deposition')
            preparation_can_go_ahead=True
            self.Logfile.insert('checking refinement stage of respective PanDDA sites...')
            panddaSites=self.db.execute_statement("select CrystalName,RefinementOutcome,PANDDA_site_event_map_mtz from panddaTable where CrystalName is '%s' and PANDDA_site_ligand_placed is 'True'" %xtal)
            self.Logfile.insert('found '+str(len(panddaSites))+' ligands')
            eventMtz=[]
            for site in panddaSites:
                if str(site[1]).startswith('5'):
                    self.Logfile.insert('site is ready for deposition')
                    eventMtz.append(str(site[2]))
                else:
                    self.Logfile.insert('site is NOT ready for deposition')
                    preparation_can_go_ahead=False

            n_eventMtz = len(eventMtz)
            if preparation_can_go_ahead:
                ModelData=self.db.execute_statement("select RefinementPDB_latest,RefinementMTZ_latest,RefinementCIF,DataProcessingPathToLogfile,RefinementProgram,CompoundCode from mainTable where CrystalName is '%s'" %xtal)
                pdb=str(ModelData[0][0])
                mtz=str(ModelData[0][1])
                cif=str(ModelData[0][2])
                if not os.path.isfile(cif):
                    if not os.path.isfile(os.path.join(self.projectDir,xtal,cif)):
                        self.Logfile.insert('cannot find ligand CIF file! Please check %s and the database!' %(os.path.join(self.projectDir,xtal)))
                        self.Logfile.insert('cannot prepare mmcif files for %s; skipping... => ERROR' %xtal)
                        continue
                log=str(ModelData[0][3])
                refSoft=str(ModelData[0][4])
                compoundID=str(ModelData[0][5])

                # if all modelled ligands are ready for deposition, we can continue

                # first get all meta-data for deposition, i.e. data_template file
                wavelength=self.prepare_data_template_for_xtal(xtal,compoundID,pdb)

                # make model mmcif
                self.make_model_mmcif(xtal,pdb,log,refSoft)

                # make SF mmcif
                self.make_sf_mmcif(xtal,mtz,eventMtz)

                # report any problems that came up
                mmcif=self.check_mmcif_files_and_update_db(xtal,n_eventMtz,wavelength)

                # add ligand cif file to mmcif
                self.add_ligand_cif_file(mmcif,cif)

                self.Logfile.insert('finished preparation of mmcif files')


            else:
                self.Logfile.insert(XChemToolTips.deposition_pandda_site_not_ready(xtal))

    def prepare_data_template_for_xtal(self,xtal,compoundID,pdb):
        # check if file exists
        if self.overwrite_existing_mmcif:
            os.chdir(os.path.join(self.projectDir,xtal))
            data_template_dict=self.db.get_deposit_dict_for_sample(xtal)

            # edit title
            data_template_dict['title']=data_template_dict['structure_title'].replace('$ProteinName',data_template_dict['Source_organism_gene']).replace('$CompoundName',compoundID)
            self.Logfile.insert('deposition title for '+xtal+': '+data_template_dict['title'])

            # get protein chains
            data_template_dict['protein_chains']=''
            chains=pdbtools(pdb).GetProteinChains()
            for item in chains:
                data_template_dict['protein_chains']+=item+','
            data_template_dict['protein_chains']=data_template_dict['protein_chains'][:-1]

            self.Logfile.insert('creating %s file for %s' %(self.data_template,xtal))
            if self.data_template.endswith('cif'):
                data_template=templates().data_template_cif(data_template_dict)
            else:
                data_template=templates().data_template_text(data_template_dict)
            f=open(os.path.join(self.projectDir,xtal,self.data_template),'w')
            f.write(data_template)
            f.close()

            wavelength=data_template_dict['radiation_wavelengths']
            return wavelength



    def make_model_mmcif(self,xtal,pdb,log,refSoft):
        if os.path.isfile(pdb) and os.path.isfile(log):
            os.chdir(os.path.join(self.projectDir,xtal))
            Cmd = ( 'source '+os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/setup.sh')+'\n'
                    +os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/bin/pdb_extract')+
                    ' -r PHENIX'
#                    ' -r %s'            %refSoft+
#                    ' -iLOG initial.log'
                    ' -iPDB %s'         %pdb+
                    ' -e MR'
                    ' -s AIMLESS'
                    ' -iLOG %s'         %log+
                    ' -iENT %s'         %self.data_template+
                    ' -o %s.mmcif > %s.mmcif.log'      %(xtal,xtal)       )

            self.Logfile.insert('running pdb_extract: '+Cmd)
            os.system(Cmd)

            # can we add here the ligand.cif?
#    '''pdb_extract -r REFMAC -iLOG initial.log -iPDB initial.pdb -e MR -s AIMLESS -iLOG aimless.log -iENT data_template.cif -o NUDT22A-x0315-model.cif'''


    def make_sf_mmcif(self,xtal,mtz,eventMtz):
        if os.path.isfile(mtz):
            os.chdir(os.path.join(self.projectDir,xtal))
            mtzin = mtz+' '
            for event in eventMtz:
                mtzin+=event+' '

            Cmd = ( 'source '+os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/setup.sh')+'\n'
                    +os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/bin/sf_convert')+
                    ' -o mmcif'
                    ' -sf %s' %mtzin+
                    ' -out %s_sf.mmcif  > %s.sf_mmcif.log' %(xtal,xtal) )

#            Cmd = ( os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/sf-convert-v1.204-prod-src/bin/sf_convert')+
#                    ' -o mmcif'
#                    ' -sf %s' %mtzin+
#                    ' -out %s_sf.mmcif > %s.sf_mmcif.log' %(xtal,xtal) )

            self.Logfile.insert('running sf_convert: '+Cmd)
            os.system(Cmd)
            os.system('/bin/rm sf_format_guess.text mtzdmp.log SF_4_validate.cif sf_information.cif')


    def check_mmcif_files_and_update_db(self,xtal,n_eventMtz,wavelength):
        os.chdir(os.path.join(self.projectDir,xtal))
        foundFiles=True

        mmcif=''
        if os.path.isfile(xtal+'.mmcif') and os.path.getsize(xtal+'.mmcif') > 20 :
            self.Logfile.insert('found '+xtal+'.mmcif')
            mmcif=os.path.join(self.projectDir,xtal,xtal+'.mmcif')
        else:
            self.Logfile.insert('cannot find '+xtal+'.mmcif; something went wrong! => ERROR')
            foundFiles=False

        if os.path.isfile(xtal+'_sf.mmcif') and os.path.getsize(xtal+'_sf.mmcif') > 20 :
            self.Logfile.insert('found '+xtal+'_sf.mmcif')
            mmcif_sf=os.path.join(self.projectDir,xtal,xtal+'_sf.mmcif')

            # now check if really all event maps are in mmcif file
            n_eventMTZ_found=-1     # set to -1 since first data block should be measured data
            for line in open(mmcif_sf):
                if line.startswith('_refln.crystal_id'):
                    n_eventMTZ_found+=1
            if n_eventMTZ_found != n_eventMtz:
                self.Logfile.insert('%s event map mtz files were specified as input, but only %s ended up in the mmcif file => ERROR' %(str(n_eventMtz),str(n_eventMTZ_found)))
                foundFiles=False
            else:
                self.Logfile.insert('%s event map mtz files were specified as input, %s ended up in the mmcif file, all well so far...' %(str(n_eventMtz),str(n_eventMTZ_found)))

            self.Logfile.insert('editing wavelength information in SF mmcif file; changing wavelength to %s' %wavelength)
            tmpText=''
            for line in open(mmcif_sf):
                if line.startswith('_diffrn_radiation_wavelength.wavelength'):
                    tmpText+='_diffrn_radiation_wavelength.wavelength   %s\n' %wavelength
                else:
                    tmpText+=line
            f=open(mmcif_sf,'w')
            f.write(tmpText)
            f.close()
            if os.path.isfile(mmcif_sf):
                self.Logfile.insert('mmcif SF file successfully update')
            else:
                self.Logfile.insert('something went wrong during the update => ERROR')
                foundFiles=False

        else:
            self.Logfile.insert('cannot find '+xtal+'_sf.mmcif; something went wrong! => ERROR')
            foundFiles=False

        if foundFiles:
            self.Logfile.insert('updating database with file locations for %s.mmcif and %s_sf.mmcif' %(xtal,xtal))
            self.db.execute_statement("update depositTable set mmCIF_model_file='%s',mmCIF_SF_file='%s' where CrystalName is '%s'" %(mmcif,mmcif_sf,xtal))
        else:
            self.Logfile.insert('could not find %s.mmcif and/or %s_sf.mmcif; removing empty files...')
            os.system('/bin/rm %s.mmcif 2> /dev/null' %xtal)
            os.system('/bin/rm %s_sf.mmcif 2> /dev/null' %xtal)

        return mmcif


    def add_ligand_cif_file(self,mmcif,ligand_cif):
        self.Logfile.insert('adding ligand cif file to %s' %mmcif)
        tmpText=''
        for line in open(mmcif):
            tmpText+=line
        for line in open(ligand_cif):
            tmpText+=line
        f=open(mmcif,'w')
        f.write(tmpText)
        f.close()

    def make_site_description(self):
        print 'hallo'
