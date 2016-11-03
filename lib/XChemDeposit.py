# last edited: 03/11/2016, 17:00

import sys
import os
import glob

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
import XChemLog


deposition = {}


def data_template():

    data_tempplate = (
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
        '<contact_author_PI_salutation=  >    !(Dr./Prof./Mr./Mrs./Ms.)\n'
        '<contact_author_PI_first_name=  >    !(e.g. John)\n'
        '<contact_author_PI_last_name=  >     !(e.g. Rodgers)\n'
        '<contact_author_PI_middle_name=  >\n'
        '<contact_author_PI_role= principal investigator/group leader>  !(or responsible scientist)\n'
        '<contact_author_PI_organization_type= academic>  !(or commercial, government, other)\n'
        '<contact_author_PI_email=  >              !(e.g.   name@host.domain.country)\n'
        '<contact_author_PI_address=  >            !(e.g. 610 Taylor road)\n'
        '<contact_author_PI_city=  >               !(e.g. Piscataway)\n'
        '<contact_author_PI_State_or_Province=  >  !(e.g.  New Jersey)\n'
        '<contact_author_PI_Zip_Code=  >           !(e.g.  08864)\n'
        '<contact_author_PI_Country=  >        !(e.g. United States, United Kindom, . )\n'
        '<contact_author_PI_fax_number=  >\n'
        '<contact_author_PI_phone_number=  >    !(e.g.  01(617) 555-1213 )\n'
        '\n'
        '2. Information about other contact authors (responsible scientist, investigator)\n'
        '\n'
        '<contact_author_id=  >              (e.g. 2)\n'
        '<contact_author_salutation=  >\n'
        '<contact_author_first_name=  >\n'
        '<contact_author_last_name=  >\n'
        '<contact_author_middle_name=  >\n'
        '<contact_author_role=  >\n'
        '<contact_author_organization_type=  >\n'
        '<contact_author_email=  >\n'
        '<contact_author_address=  >\n'
        '<contact_author_city=  >\n'
        '<contact_author_State_or_Province=  >\n'
        '<contact_author_Zip_Code=  >\n'
        '<contact_author_Country=  >\n'
        '<contact_author_fax_number=  >\n'
        '<contact_author_phone_number=  >\n'
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
        '<Release_status_for_coordinates=  >      !(e.g. HOLD FOR PUBLICATION)\n'
        '<Release_status_for_structure_factor=  > !(e.g. HOLD FOR PUBLICATION)\n'
        '<Release_status_for_sequence=  >         !(RELEASE NOW or  HOLD FOR RELEASE)\n'
        '\n'
        '================CATEGORY 3:   Title=======================================\n'
        'Enter the title for the structure\n'
        '\n'
        '<structure_title=  >     !(e.g. Crystal Structure Analysis of the B-DNA)\n'
        '<structure_details=  >\n'
        '\n'
        '================CATEGORY 4: Authors of Structure============================\n'
        'Enter authors of the deposited structures (at least one author)\n'
        '\n'
        '<structure_author_name=  >  !(e.g. Surname, F.M.)\n'
        '<structure_author_name=  >\n'
        '<structure_author_name=  >\n'
        '<structure_author_name=  >\n'
        '<structure_author_name=  >\n'
        '<structure_author_name=  >\n'
        '<structure_author_name=  >\n'
        '<structure_author_name=  >\n'
        '<structure_author_name=  >\n'
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
        '<primary_citation_author_name=  >    !(e.g. Surname, F.M.)\n'
        '<primary_citation_author_name=  >\n'
        '<primary_citation_author_name=  >\n'
        '<primary_citation_author_name=  >\n'
        '<primary_citation_author_name=  >\n'
        '<primary_citation_author_name=  >\n'
        '<primary_citation_author_name=  >\n'
        '<primary_citation_author_name=  >\n'
        '\n'
        '...add more name if needed...\n'
        '\n'
        'Enter journal information of the primary citation\n'
        '<primary_citation_id= primary>\n'
        '<primary_citation_journal_abbrev=  >     (e.g. To be published)\n'
        '<primary_citation_title=  >\n'
        '<primary_citation_year=  >\n'
        '<primary_citation_journal_volume=  >\n'
        '<primary_citation_page_first=  >\n'
        '<primary_citation_page_last=  >\n'
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
        '<molecule_name=  >       (e.g.  RNA Hammerhead Ribozyme )\n'
        '<molecule_type= polymer >    (e.g. polymer , non-polymer, macrolide  )\n'
        '<molecule_source_method= >   (e.g. man , nat, syn)\n'
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
        '<Fragment_name=  >             (e.g. ligand binding domain, hairpin)\n'
        '<Specific_mutation=  >         (e.g. C280S)\n'
        '<Enzyme_Comission_number=  >   (if known: e.g. 2.7.7.7)\n'
        '\n'
2. For entity 2
<Molecular_entity_id=  >       (e.g.  2 )
<Fragment_name=  >
<Specific_mutation=  >
<Enzyme_Comission_number=  >

...(add more group if needed)...

================CATEGORY 8:   Genetically Manipulated Source=============
Enter data in the genetically manipulated source category

  If the biomolecule has been genetically manipulated, describe its
  source and expression system here.

1. For entity 1
<Manipulated_entity_id= 1 >               !(e.g. 1 )
<Source_organism_scientific_name=  >      !(e.g. Homo sapiens)
<Source_organism_gene=  >                 (e.g. RPOD, ALKA...)
<Source_organism_strain=  >               (e.g. BH10 ISOLATE, K-12...)
<Expression_system_scientific_name=  >    (e.g. Escherichia coli)
<Expression_system_strain=  >	          (e.g. BL21(DE3))
<Expression_system_vector_type=  >	  (e.g. plasmid)
<Expression_system_plasmid_name=  >       (e.g. pET26)
<Manipulated_source_details=  >           (any other relevant information)

2. For entity 2

...(add more group if needed)...

================CATEGORY 9:   Natural Source (optional) ===================
Enter data in the natural source category  (if applicable)

  If the biomolecule was derived from a natural source, describe it here.

1. For entity 1
<natural_source_entity_id=  >          (e.g. 1, 2..)
<natural_source_scientific_name=  >    (e.g. Homo sapiens)
<natural_source_organism_strain=  >    (e.g. DH5a , BMH 71-18)
<natural_source_details=  >            (e.g. organ, tissue, cell ..)

2. For entity 2

...(add more group if needed)...

================CATEGORY 10:  Synthetic Source (optional)==================
If the biomolecule has not been genetically manipulated or synthesized,
describe its source here.

1. For entity 1
<synthetic_source_entity_id=  >          (e.g. 1,2. )
<synthetic_source_description=  >      (if known)


2. For entity 2

...(add more group if needed)...


================CATEGORY 11:   Keywords===================================
Enter a list of keywords that describe important features of the deposited
structure.

  Example: beta barrel, protein-DNA complex, double helix, hydrolase, etc.

<structure_keywords=  >  !(e.g. beta barrel)

================CATEGORY 12:   Biological Assembly ======================
Enter data in the biological assembly category (if applicable)

Enter the number of polymer chains that form the assembly in solution

<biological_assembly_chain_number=  >  !(e.g.  1 for monomer, 2 for dimer ..)

================CATEGORY 13:   Methods and Conditions=====================
Enter the crystallization conditions for each crystal

1. For crystal 1:
<crystal_number= 1 >	            (e.g. 1, )
<crystallization_method=  >      (e.g. BATCH MODE, EVAPORATION)
<crystallization_pH=  >          (e.g. 7.5 ...)
<crystallization_temperature=  > (e.g. 298) (in Kelvin)
<crystallization_details=  >     (e.g. PEG 4000, NaCl etc.)

...(add more crystal groups if needed)...

================CATEGORY 14:   Radiation Source (experiment)============
Enter the details of the source of radiation, the X-ray generator,
and the wavelength for each diffraction.

1. For experiment 1:
<radiation_experiment= 1 >      !(e.g. 1, 2, ...)
<radiation_source=  >           !(e.g. SYNCHROTRON, ROTATING ANODE ..)
<radiation_source_type=  >      !(e.g. NSLS BEAMLINE X8C ..)
<radiation_wavelengths=  >       !(e.g. 1.502, or a list 0.987,0.988 ..)
<radiation_detector=  >         !(e.g. CCD, AREA DETECTOR, IMAGE PLATE ..)
<radiation_detector_type=  >     !(e.g. ADSC QUANTUM 1,  ..)
<radiation_detector_details=  >    (e.g. mirrors...)
<data_collection_date=  >             !(e.g. 2004-01-07)
<data_collection_temperature=  >      !(e.g. 100 for crystal  1:)
<data_collection_protocol=  >          !(e.g. SINGLE WAVELENGTH, MAD, ...)
<data_collection_monochromator=  >     (e.g. GRAPHITE, Ni FILTER ...)
<data_collection_monochromatic_or_laue=  M >  !(default M, give L if Laue diffr.)


....(add more experiment group if needed)....

================CATEGORY 15:   refinement details (optional)============
Enter the details of the structure refinement. (if applicable)

<refinement_detail=   >
<refinement_start_model=   >    (e.g. pdbid 100D)

================CATEGORY 16:   database (optional)======================
Enter the database name for each molecule (entity), (IF KNOWN).

1. For entity 1
<database_entity_id= 1  >  (e.g. 1 )
<database_name=  >  (e.g. BMCD, BMRB, EMDB, PDB, NDB, TargetTrack )
<database_code=  >  (e.g. 1ABC, 100D, TNKS2_HUMAN )
<database_accession=   >  (e.g. 100D, Q9H2K2  )


...(add more group if needed)...

================CATEGORY 17:   Ligand binding (optional)==================

<binding_assay_id=   >     !(A unique identifier such as 1,2,..)
<binding_assay_target_sequence_one_letter_code=   >  (Chemical sequence if known)
<binding_assay_ligand_descriptor_type=   >  (e.g. SMILES, SMILES_CANONICAL, InChI,InChIKey)
<binding_assay_ligand_descriptor=   >   (e.g. Cc1cccc(c1)C1COc2cc(O)c(O)cc2C1 )
<binding_assay_assay_type=   >        (Type of binding assay. e.g. 'competitive binding')
<binding_assay_assay_value_type=   > (e.g. IC50, EC50,Ki,Kd)
<binding_assay_assay_value=   >     (The value measured. e.g. 8300.0)
<binding_assay_assay_pH=   >       (pH value at which the assay was performed. e.g. 6.4)
<binding_assay_assay_temperature=   > (temperature (K) at which the assay was performed. e.g.273)
<binding_assay_details=   >  (details of the measurement).



================CATEGORY 18:   Structure Genomic (optional)==============
If it is the structure genome's project, give the information

<SG_project_id=  1>
<SG_project_name=  >        (e.g. PSI, Protein Structure Initiative)
<full_name_of_SG_center=  >   (e.g. Berkeley Structural Genomic Center)


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++ Categories below were extracted from the coordinate. Please check. +++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

================CATEGORY 19:   Crystallographic Data=======================
Enter crystallographic data

<space_group = P 43 21 2> (use International Table conventions)
<unit_cell_a     =   71.467  >
<unit_cell_b     =   71.467  >
<unit_cell_c     =  150.264  >
<unit_cell_alpha =  90.00  >
<unit_cell_beta  =  90.00  >
<unit_cell_gamma =  90.00  >


================CATEGORY 20:   Sequence Information =======================
Enter one letter sequence for each polymeric entity in asymmetric unit

--------------------------------------------------------------------------
			  SOME DEFINITIONS
     An ENTITY is defined as any unique molecule present in the asymmetric
     unit. Each unique biological polymer (protein or nucleic acids) in the
     structure is considered an entity. Thus, if there are five copies of
     a single protein in the asymmetric unit, the molecular entity is still
     only one. Water and non-polymers like ions, ligands and sugars are
     also entities.

     Here we only consider the sequences of polymeric entities (protein or
     nucleic acid).

	        GUIDELINES FOR COMPLETING THIS CATEGORY

      * The unique chemical sequences are extract from the coordinate file.
      Never give a TER cards in the pdb file inside of a complete entity.

      * In a PDB or mmCIF format file, all residues of a single polymeric
      entity should have one chain ID. Multiple copies of the same entity
      should each be assigned a unique chain ID. The multiple chain IDs
      should be separated by commas as 'A,B,C,...'.

      * Four question marks '????' are used to denote a polymer chain breaks,
      which may be caused by missing residues due to unobsered electron density
      or linkage problem due to poor geometry. Replace these question marks
      with the sequence of residues missing from the coordinates. Delete these
      question marks if the break is due to the poor geometry. Also add any
      residues missing from the N- and/or C-termini here.

      * If there are non-standard residues in the coordinates, this program
      lists them according to the three letter code used in the coordinate
      file as (ABC).

      * If any residue was modeled as Ala or Gly due to lack of the side-chain
      density, the sequence extracted here will represent them as A or G
      respectively. Correct this to the original sequence that was present in
      the crystal.


----------------------------------------------------------------------------

  Below is the one letter chemical sequence extracted from your PDB coordinate
  file. The molecular entities are grouped and listed together.

  PLEASE CHECK THE SEQUENCE of each entity carefully and modify it, as necessary.
  Make sure that you REVIEW THE FOLLOWING:
   * chain breaks due to missing residues,
   * missing residues in the N- and/or C-termini,
   * non-standard residues and
   * cases of residues modeled as Ala or Gly due to missing side-chain density.

<molecule_entity_id=1 >
<molecule_entity_type=polypeptide(L) >
<molecule_one_letter_sequence=
AQNPNCNIMIFHPTKEEFNDFDKYIAYMESQGAHRAGLAKIIPPKEWKARETYDNISEILIATPLQQVAS
GRAGVFTQYHKKKKAMTVGEYRHLANSKKYQTPPHQNFEDLERKYWKNRIYNSPIYGADISGSLFDENTK
QWNLGHLGTIQDLLEKECGVVIEGVNTPYLYFGMWKTTFAWHTEDMDLYSINYLHLGEPKTWYVVPPEHG
QRLERLARELFPGSSRGCGAFLRHKVALISPTVLKENGIPFNRITQEAGEFMVTFPYGYHAGFNHGFNCA
EAINFATPRWIDYGKMASQCSCGEARVTFSMDAFVRILQPERYDLWKRGQD >
< molecule_chain_id=A >
< target_DB_id=  > (if known)
< sequence_database_id=  > (if known)
< sequence_database_name=  > (if known)

<molecule_entity_id=  >
<molecule_entity_type=  >
<molecule_one_letter_sequence= >
<molecule_chain_id=  >
<target_DB_id=  >  (if known)
<sequence_database_id=   > (if known)
<sequence_database_name=   > (if known)




=====================================END==================================
