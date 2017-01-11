# last edited: 10/01/2017, 15:00

import sys
import os
import glob

from PyQt4 import QtGui, QtCore, QtWebKit

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
import XChemLog
import XChemDB

deposition = {}


def data_template(depositDict):

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
        +molecule_one_letter_sequence+
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

def data_template_cif(depositDict):

    structure_author_name=''
    for name in depositDict['structure_author_name'].split(';'):
        structure_author_name+='<structure_author_name=  %s>\n' %name

    primary_citation_author_name=''
    for name in depositDict['primary_citation_author_name'].split(';'):
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
        '1 polymer     man "NUDT22"\n'
        '#\n'
        'loop_\n'
        '_entity_poly.entity_id\n'
        '_entity_poly.type\n'
        '_entity_poly.pdbx_seq_one_letter_code\n'
        '_entity_poly.pdbx_strand_id\n'
        '_entity_poly.pdbx_seq_db_id\n'
        '_entity_poly.pdbx_seq_db_name\n'
        '1 "polypeptide(L)"\n'
        +molecule_one_letter_sequence+
#        ';MDPEVTLLLQCPGGGLPQEQIQAELSPAHDRRPLPGGDEAITAIWETRLKAQPWLFDAPK\n'
#        'FRLHSATLAPIGSRGPQLLLRLGLTSYRDFLGTNWSSSAAWLRQQGATDWGDTQAYLADP\n'
#        'LGVGAALATADDFLVFLRRSRQVAEAPGLVDVPGGHPEPQALCPGGSPQHQDLAGQLVVH\n'
#        'ELFSSVLQEICDEVNLPLLTLSQPLLLGIARNETSAGRASAEFYVQCSLTSEQVRKHYLS\n'
#        'GGPEAHESTGIFFVETQNVQRLLETEMWAELCPSAKGAIILYNRVQGSPTGAALGSPALL\n'
#        'PPL\n'
        ';\n'
        'A Q9BRQ3 UNP\n'
        '#\n'
        'loop_\n'
        '_entity_src_gen.entity_id\n'
        '_entity_src_gen.gene_src_strain\n'
        '_entity_src_gen.pdbx_gene_src_scientific_name\n'
        '_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id\n'
        '_entity_src_gen.pdbx_host_org_scientific_name\n'
        '_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id\n'
        '1 ? "%s" %s  "%s" %s\n'                %(depositDict['Source_organism_scientific_name'],'?',depositDict['Expression_system_scientific_name'],'?')+
        '#\n'
        '#\n'
        '_pdbx_contact_author.id                  1\n'
        "_pdbx_contact_author.address_1           '%s'\n"                           %depositDict['contact_author_PI_address']+
#        '_pdbx_contact_author.address_2           "XXX Institute/Company"\n'        XXX missing XXX
        '_pdbx_contact_author.city                %s\n'                             %depositDict['contact_author_PI_city']+
        "_pdbx_contact_author.state_province      '%s'\n"                           %depositDict['contact_author_PI_State_or_Province']+
        '_pdbx_contact_author.postal_code         %s\n'                             %depositDict['contact_author_PI_Zip_Code']+
        '_pdbx_contact_author.email               %s\n'                             %depositDict['contact_author_PI_email']+
        '_pdbx_contact_author.name_first          %s\n'                             %depositDict['contact_author_PI_first_name']+
        '_pdbx_contact_author.name_last           %s\n'                             %depositDict['contact_author_PI_last_name']+
        '_pdbx_contact_author.country             "%s"\n'                           %depositDict['contact_author_PI_Country']+
        '_pdbx_contact_author.phone               %s\n'                             %depositDict['contact_author_PI_phone_number']+
        '_pdbx_contact_author.role                "%s"\n'                           %depositDict['contact_author_PI_role']+
        '_pdbx_contact_author.organization_type   %s\n'                             %depositDict['contact_author_PI_organization_type']+
        '#\n'
        '_pdbx_contact_author.id                  2\n'
        "_pdbx_contact_author.address_1           '%s'\n"                           %depositDict['contact_author_address']+
#        '_pdbx_contact_author.address_2           "XXX Institute/Company"\n'        XXX missing XXX
        '_pdbx_contact_author.city                %s\n'                             %depositDict['contact_author_city']+
        "_pdbx_contact_author.state_province      '%s'\n"                           %depositDict['contact_author_State_or_Province']+
        '_pdbx_contact_author.postal_code         %s\n'                             %depositDict['contact_author_Zip_Code']+
        '_pdbx_contact_author.email               %s\n'                             %depositDict['contact_author_email']+
        '_pdbx_contact_author.name_first          %s\n'                             %depositDict['contact_author_first_name']+
        '_pdbx_contact_author.name_last           %s\n'                             %depositDict['contact_author_last_name']+
        '_pdbx_contact_author.country             "%s"\n'                           %depositDict['contact_author_Country']+
        '_pdbx_contact_author.phone               %s\n'                             %depositDict['contact_author_phone_number']+
        '_pdbx_contact_author.role                "%s"\n'                           %depositDict['contact_author_role']+
        '_pdbx_contact_author.organization_type   %s\n'                             %depositDict['contact_author_organization_type']+
        '#\n'
        'loop_\n'
        '_audit_author.name\n'
        "'Krojer, T.'\n"
        "'Von Delft, F.'\n"
        '#\n'
        '_citation.id                        primary\n'
        '_citation.title                     "your citation title"\n'
        '_citation.journal_abbrev            "To Be Published"\n'
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
        ';Your entry title\n'
        ';\n'
        '#\n'
        '_struct_keywords.entry_id        UNNAMED\n'
        '_struct_keywords.text            "your keyword"\n'
        '#\n'
        '_pdbx_struct_assembly_depositor_info.id                   1\n'
        "_pdbx_struct_assembly_depositor_info.method_details       'gel filtration'\n"
        '_pdbx_struct_assembly_depositor_info.oligomeric_count     3\n'
        '#\n'
        '#\n'
        )

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

def find_apo_structures(panddaDir):

    # first check if structure is already present in DB and if so if all the
    # information concur

    # need to update pandda directory for every exported structure so that
    # we know where to look for the pandda.log file that contains the relevant information

    # update CrystalName_of_pandda_input in DB

    # in DB: update StructureType field accordingly

    # newer pandda versions seem to have severl copies of pandda.log with names like
    # pandda-2016-09-01-2139.log
    panddaLog=glob.glob(os.path.join(panddaDir,'pandda*log'))
    panddaLog.sort(key=os.path.getmtime)

    panddaVersion='unknown'
    readindApoStructures = False
    apoStructures = []
    apoStructureDict = {}
    for files in panddaLog:
        for line in open(files):
            if line.startswith('-  Pandda Version'):
                if len(line.split()) >= 4:
                    panddaVersion=line.split()[3]
            if 'No Statistical Maps Found:' in line:
                readindApoStructures=True
            if readindApoStructures:
                if 'Pickling Object: processed_datasets' in line:
                    if line.split() >= 2:
                        # e.g. line.split() -> ['Pickling', 'Object:', 'processed_datasets/NUDT22A-x0055/pickles/dataset.pickle']
                        xtal=line.split()[2].split('/')[1]
                        if os.path.isfile(os.path.join(panddaDir,'processed_datasets',xtal,xtal+'-x0785-pandda-input.pdb')):
                            apoStructures.append(xtal)
            if 'Pre-existing statistical maps (from previous runs) have been found and will be reused:' in line:
                readindApoStructures=False
        apoStructureDict[panddaDir]=apoStructures

    return apoStructureDict


#class update_depositTable(QtCore.QThread):
#    def __init__(self,deposit_dict,database,xce_logfile):
#        QtCore.QThread.__init__(self)
#        self.deposit_dict=deposit_dict
#        self.database=database
#        self.Logfile=XChemLog.updateLog(xce_logfile)
#        self.db=XChemDB.data_source(database)
#        self.header,self.data=self.db.load_samples_from_data_source()
#        self.xtal_db_dict={}
#        self.update_xtal_db_dict()
#
#    def update_xtal_db_dict(self):
#        sampleID_column=0
#        for n,entry in enumerate(self.header):
#            if entry=='CrystalName':
#                sampleID_column=n
#                break
#        for line in self.data:
#            if str(line[sampleID_column]) != '':
#                db_dict={}
#                for n,entry in enumerate(line):
#                    if n != sampleID_column:
#                        db_dict[str(self.header[n])]=str(entry)
#                self.xtal_db_dict[str(line[sampleID_column])]=db_dict
#        self.Logfile.insert('read all samples in mainTable')
#
#    def run(self):
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

class update_deposition_table:

    def __init__(self,database):
        self.db=XChemDB.data_source(database)

    def PanDDA_models_to_deposit(self):
        panddaModels=self.db.execute_statement('select CrystalName,PANDDA_site_index,RefinementOutcome,PANDDA_site_ligand_placed,PANDDA_site_name from panddaTable')
        toDeposit={}
        mismatch={}

        for item in panddaModels:
            xtalID=str(item[0])
            site_index=str(item[1])
            RefinementStage=str(item[2])
            ligand_placed=str(item[3])
            siteName=str(item[4])
            if ligand_placed.replace(' ','').lower()=='true':
                if xtalID in toDeposit:
                    if not RefinementStage.startswith('3'):
                        if xtalID not in mismatch:
                            mismatch[xtalID]=[]
                        mismatch[xtalID].append([xtalID,site_index])
                        continue

                if RefinementStage.startswith('3'):
                    if xtalID not in toDeposit:
                        toDeposit[xtalID]=[]
                    toDeposit[xtalID].append([xtalID,site_index,siteName,RefinementStage])

        return toDeposit,mismatch


