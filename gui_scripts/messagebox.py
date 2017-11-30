import sys, os
from PyQt4 import QtGui, QtCore, QtWebKit

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'gui_scripts'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))

import layout
import XChemToolTips
import XChemMain


def clearWidgetsAndLayout(layout):
#    i=0
#    while layout.count():
    for i in reversed(range(layout.count())):
        print i
        widgetToRemove = layout.itemAt( i ).widget()
        # remove it from the layout list
        layout.removeWidget( widgetToRemove )
        # remove it from the gui
        widgetToRemove.setParent( None )
#        child = layout.takeAt(0)
#        if child.widget() is not None:
#            print i,child.widget()
##            child.widget().deleteLater()
#        elif child.layout() is not None:
#            print i,child.layout()
##            clearLayout(child.layout())
#        i += 1
    print 'here'

class DepositionDetails():
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
        ################################################################################################################
        #                                                                                                              #
        #                                                DEPOSITION Details                                            #
        #                                                                                                              #
        ################################################################################################################

#        clearLayout(xce_object.messageBoxLayout)

        vbox = QtGui.QVBoxLayout()

        deposit_tab_widget = QtGui.QTabWidget()
        deposit_tab_list = ['Contact',
                            'General',
                            'Authors',
                            'Citation',
                            'Molecule',
                            'Misc',
                            'Methods',
                            'Software']

        deposit_tab_dict = {}
        for page in deposit_tab_list:
            tab = QtGui.QWidget()
            vb = QtGui.QVBoxLayout(tab)
            deposit_tab_widget.addTab(tab, page)
            deposit_tab_dict[page] = [tab, vb]

        ## PI and scientist info
        vb = QtGui.QVBoxLayout()
        hbox = QtGui.QHBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Principal Investigator'), 0, 0)

        grid.addWidget(QtGui.QLabel('Salutation'), 1, 0)
        xce_object.contact_author_PI_salutation = QtGui.QLineEdit()
        xce_object.contact_author_PI_salutation.setText('Dr.')
        xce_object.contact_author_PI_salutation.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_salutation, 1, 1)

        grid.addWidget(QtGui.QLabel('First name'), 2, 0)
        xce_object.contact_author_PI_first_name = QtGui.QLineEdit()
        xce_object.contact_author_PI_first_name.setText('')
        xce_object.contact_author_PI_first_name.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_first_name, 2, 1)

        grid.addWidget(QtGui.QLabel('Last name'), 3, 0)
        xce_object.contact_author_PI_last_name = QtGui.QLineEdit()
        xce_object.contact_author_PI_last_name.setText('')
        xce_object.contact_author_PI_last_name.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_last_name, 3, 1)

        grid.addWidget(QtGui.QLabel('Middle name'), 4, 0)
        xce_object.contact_author_PI_middle_name = QtGui.QLineEdit()
        xce_object.contact_author_PI_middle_name.setText('')
        xce_object.contact_author_PI_middle_name.setFixedWidth(200)
        xce_object.contact_author_PI_middle_name.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.contact_author_PI_middle_name, 4, 1)

        grid.addWidget(QtGui.QLabel('PI role'), 5, 0)
        xce_object.contact_author_PI_role = QtGui.QComboBox()
        PIroles = ['group leader', 'principal investigator', 'investigator']
        for item in PIroles: xce_object.contact_author_PI_role.addItem(item)
        grid.addWidget(xce_object.contact_author_PI_role, 5, 1)

        grid.addWidget(QtGui.QLabel('Organization type'), 6, 0)
        xce_object.contact_author_PI_organization_type = QtGui.QComboBox()
        Organizations = ['academic', 'commercial', 'government']
        for item in Organizations: xce_object.contact_author_PI_organization_type.addItem(item)
        grid.addWidget(xce_object.contact_author_PI_organization_type, 6, 1)

        grid.addWidget(QtGui.QLabel('Organization Name'), 7, 0)
        xce_object.contact_author_PI_organization_name = QtGui.QLineEdit()
        xce_object.contact_author_PI_organization_name.setText('')
        xce_object.contact_author_PI_organization_name.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_organization_name, 7, 1)

        grid.addWidget(QtGui.QLabel('Email'), 8, 0)
        xce_object.contact_author_PI_email = QtGui.QLineEdit()
        xce_object.contact_author_PI_email.setText('')
        xce_object.contact_author_PI_email.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_email, 8, 1)

        grid.addWidget(QtGui.QLabel('Street'), 9, 0)
        xce_object.contact_author_PI_address = QtGui.QLineEdit()
        xce_object.contact_author_PI_address.setText('')
        xce_object.contact_author_PI_address.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_address, 9, 1)

        grid.addWidget(QtGui.QLabel('City'), 10, 0)
        xce_object.contact_author_PI_city = QtGui.QLineEdit()
        xce_object.contact_author_PI_city.setText('')
        xce_object.contact_author_PI_city.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_city, 10, 1)

        grid.addWidget(QtGui.QLabel('State'), 11, 0)
        xce_object.contact_author_PI_State_or_Province = QtGui.QLineEdit()
        xce_object.contact_author_PI_State_or_Province.setText('')
        xce_object.contact_author_PI_State_or_Province.setFixedWidth(200)
        xce_object.contact_author_PI_State_or_Province.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.contact_author_PI_State_or_Province, 11, 1)

        grid.addWidget(QtGui.QLabel('ZIP code'), 12, 0)
        xce_object.contact_author_PI_Zip_Code = QtGui.QLineEdit()
        xce_object.contact_author_PI_Zip_Code.setText('')
        xce_object.contact_author_PI_Zip_Code.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_Zip_Code, 12, 1)

        grid.addWidget(QtGui.QLabel('Country'), 13, 0)
        xce_object.contact_author_PI_Country = QtGui.QLineEdit()
        xce_object.contact_author_PI_Country.setText('')
        xce_object.contact_author_PI_Country.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_Country, 13, 1)

        grid.addWidget(QtGui.QLabel('Phone'), 14, 0)
        xce_object.contact_author_PI_phone_number = QtGui.QLineEdit()
        xce_object.contact_author_PI_phone_number.setText('')
        xce_object.contact_author_PI_phone_number.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_PI_phone_number, 14, 1)

        frame.setLayout(grid)
        hbox.addWidget(frame)

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Responsible Scientist'), 0, 0)

        grid.addWidget(QtGui.QLabel('Salutation'), 1, 0)
        xce_object.contact_author_salutation = QtGui.QLineEdit()
        xce_object.contact_author_salutation.setText('Dr.')
        xce_object.contact_author_salutation.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_salutation, 1, 1)

        grid.addWidget(QtGui.QLabel('First name'), 2, 0)
        xce_object.contact_author_first_name = QtGui.QLineEdit()
        xce_object.contact_author_first_name.setText('')
        xce_object.contact_author_first_name.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_first_name, 2, 1)

        grid.addWidget(QtGui.QLabel('Last name'), 3, 0)
        xce_object.contact_author_last_name = QtGui.QLineEdit()
        xce_object.contact_author_last_name.setText('')
        xce_object.contact_author_last_name.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_last_name, 3, 1)

        grid.addWidget(QtGui.QLabel('Middle name'), 4, 0)
        xce_object.contact_author_middle_name = QtGui.QLineEdit()
        xce_object.contact_author_middle_name.setText('')
        xce_object.contact_author_middle_name.setFixedWidth(200)
        xce_object.contact_author_middle_name.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.contact_author_middle_name, 4, 1)

        grid.addWidget(QtGui.QLabel('Role'), 5, 0)

        xce_object.contact_author_role = QtGui.QComboBox()
        ScientistRoles = ['responsible scientist', 'investigator']
        for item in ScientistRoles: xce_object.contact_author_role.addItem(item)
        grid.addWidget(xce_object.contact_author_role, 5, 1)

        grid.addWidget(QtGui.QLabel('Organization type'), 6, 0)

        xce_object.contact_author_organization_type = QtGui.QComboBox()
        for item in Organizations: xce_object.contact_author_organization_type.addItem(item)
        grid.addWidget(xce_object.contact_author_organization_type, 6, 1)

        grid.addWidget(QtGui.QLabel('Organization Name'), 7, 0)
        xce_object.contact_author_organization_name = QtGui.QLineEdit()
        xce_object.contact_author_organization_name.setText('')
        xce_object.contact_author_organization_name.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_organization_name, 7, 1)

        grid.addWidget(QtGui.QLabel('Email'), 8, 0)
        xce_object.contact_author_email = QtGui.QLineEdit()
        xce_object.contact_author_email.setText('')
        xce_object.contact_author_email.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_email, 8, 1)

        grid.addWidget(QtGui.QLabel('Street'), 9, 0)
        xce_object.contact_author_address = QtGui.QLineEdit()
        xce_object.contact_author_address.setText('')
        xce_object.contact_author_address.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_address, 9, 1)

        grid.addWidget(QtGui.QLabel('City'), 10, 0)
        xce_object.contact_author_city = QtGui.QLineEdit()
        xce_object.contact_author_city.setText('')
        xce_object.contact_author_city.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_city, 10, 1)

        grid.addWidget(QtGui.QLabel('State'), 11, 0)
        xce_object.contact_author_State_or_Province = QtGui.QLineEdit()
        xce_object.contact_author_State_or_Province.setText('')
        xce_object.contact_author_State_or_Province.setFixedWidth(200)
        xce_object.contact_author_State_or_Province.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.contact_author_State_or_Province, 11, 1)

        grid.addWidget(QtGui.QLabel('ZIP code'), 12, 0)
        xce_object.contact_author_Zip_Code = QtGui.QLineEdit()
        xce_object.contact_author_Zip_Code.setText('')
        xce_object.contact_author_Zip_Code.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_Zip_Code, 12, 1)

        grid.addWidget(QtGui.QLabel('Country'), 13, 0)
        xce_object.contact_author_Country = QtGui.QLineEdit()
        xce_object.contact_author_Country.setText('')
        xce_object.contact_author_Country.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_Country, 13, 1)

        grid.addWidget(QtGui.QLabel('Phone'), 14, 0)
        xce_object.contact_author_phone_number = QtGui.QLineEdit()
        xce_object.contact_author_phone_number.setText('')
        xce_object.contact_author_phone_number.setFixedWidth(200)
        grid.addWidget(xce_object.contact_author_phone_number, 14, 1)

        frame.setLayout(grid)
        hbox.addWidget(frame)

        vb.addLayout(hbox)
        vb.addWidget(QtGui.QLabel(XChemToolTips.deposition_interface_note()))
        vb.addStretch(1)

        deposit_tab_dict['Contact'][1].addLayout(vb)

        ## release status
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Release status'), 0, 0)

        grid.addWidget(QtGui.QLabel('Release Status for sequence'), 4, 0)

        xce_object.Release_status_for_sequence = QtGui.QComboBox()
        codeStatus = ['RELEASE NOW', 'HOLD FOR RELEASE']
        for item in codeStatus: xce_object.Release_status_for_sequence.addItem(item)
        grid.addWidget(xce_object.Release_status_for_sequence, 4, 1)

        grid.addWidget(QtGui.QLabel('Release Status for coordinates/ SF'), 8, 0)
        xce_object.Release_status_for_coordinates = QtGui.QComboBox()
        coordStatus = ['RELEASE NOW', 'HOLD FOR PUBLICATION', 'HOLD FOR 4 WEEKS', 'HOLD FOR 6 MONTHS',
                       'HOLD FOR 1 YEAR']
        for item in coordStatus: xce_object.Release_status_for_coordinates.addItem(item)
        grid.addWidget(xce_object.Release_status_for_coordinates, 8, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Title & Details'), 0, 0)
        note = (
            'Note: supported wildcards: $ProteinName,$CompoundName; e.g. "Crystal Structure of human JMJD2D in complex with N2317a"')
        grid.addWidget(QtGui.QLabel(note), 1, 0)

        grid.addWidget(QtGui.QLabel('Group deposition title'), 2, 0)
        xce_object.group_deposition_title = QtGui.QLineEdit()
        xce_object.group_deposition_title.setText('PanDDA analysis group deposition')
        xce_object.group_deposition_title.setFixedWidth(600)
        #        xce_object.group_deposition_title.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.group_deposition_title, 2, 1)

        grid.addWidget(QtGui.QLabel('Description'), 3, 0)
        xce_object.group_description = QtGui.QLineEdit()
        xce_object.group_description.setText(
            'XDomainX of XOrganismX $ProteinName screened against the XXX Fragment Library by X-ray Crystallography at the XChem facility of Diamond Light Source beamline I04-1')
        xce_object.group_description.setFixedWidth(600)
        grid.addWidget(xce_object.group_description, 3, 1)

        grid.addWidget(QtGui.QLabel('Structure Title (ligand bound)'), 4, 0)
        xce_object.structure_title = QtGui.QLineEdit()
        xce_object.structure_title.setText('Crystal Structure of $ProteinName in complex with $CompoundName')
        xce_object.structure_title.setFixedWidth(600)
        grid.addWidget(xce_object.structure_title, 4, 1)

        note = ('\n\nApo Structure:\nonly use if you want to deposit PanDDA models!')
        grid.addWidget(QtGui.QLabel(note), 6, 0)

        grid.addWidget(QtGui.QLabel('Structure Title (apo)'), 7, 0)
        xce_object.structure_title_apo = QtGui.QLineEdit()
        xce_object.structure_title_apo.setText(
            'Crystal Structure of $ProteinName after initial refinement with no ligand modelled (structure $n)')
        xce_object.structure_title_apo.setFixedWidth(600)
        grid.addWidget(xce_object.structure_title_apo, 7, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['General'][1].addLayout(vb)

        ## authors
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Deposition authors (e.g. Surname, F.M.)'), 0, 0)

        xce_object.structure_author_name_List = []

        for column in range(0, 2):
            for row in range(1, 15):
                structure_author_name = QtGui.QLineEdit()
                structure_author_name.setText('')
                structure_author_name.setFixedWidth(300)
                grid.addWidget(structure_author_name, row, column)
                xce_object.structure_author_name_List.append(structure_author_name)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Authors'][1].addLayout(vb)

        ## primary citation
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Primary Citation'), 0, 0)

        grid.addWidget(QtGui.QLabel('ID'), 1, 0)
        xce_object.primary_citation_id = QtGui.QLineEdit()
        xce_object.primary_citation_id.setText('primary')
        xce_object.primary_citation_id.setFixedWidth(500)
        grid.addWidget(xce_object.primary_citation_id, 1, 1)

        grid.addWidget(QtGui.QLabel('Journal'), 2, 0)
        xce_object.primary_citation_journal_abbrev = QtGui.QLineEdit()
        xce_object.primary_citation_journal_abbrev.setText('To be published')
        xce_object.primary_citation_journal_abbrev.setFixedWidth(500)
        grid.addWidget(xce_object.primary_citation_journal_abbrev, 2, 1)

        grid.addWidget(QtGui.QLabel('Title'), 3, 0)
        xce_object.primary_citation_title = QtGui.QLineEdit()
        xce_object.primary_citation_title.setText('')
        xce_object.primary_citation_title.setFixedWidth(500)
        xce_object.primary_citation_title.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.primary_citation_title, 3, 1)

        grid.addWidget(QtGui.QLabel('Year'), 4, 0)
        xce_object.primary_citation_year = QtGui.QLineEdit()
        xce_object.primary_citation_year.setText('')
        xce_object.primary_citation_year.setFixedWidth(500)
        xce_object.primary_citation_year.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.primary_citation_year, 4, 1)

        grid.addWidget(QtGui.QLabel('Volume'), 5, 0)
        xce_object.primary_citation_journal_volume = QtGui.QLineEdit()
        xce_object.primary_citation_journal_volume.setText('')
        xce_object.primary_citation_journal_volume.setFixedWidth(500)
        xce_object.primary_citation_journal_volume.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.primary_citation_journal_volume, 5, 1)

        grid.addWidget(QtGui.QLabel('Page, first'), 6, 0)
        xce_object.primary_citation_page_first = QtGui.QLineEdit()
        xce_object.primary_citation_page_first.setText('')
        xce_object.primary_citation_page_first.setFixedWidth(500)
        xce_object.primary_citation_page_first.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.primary_citation_page_first, 6, 1)

        grid.addWidget(QtGui.QLabel('Page, last'), 7, 0)
        xce_object.primary_citation_page_last = QtGui.QLineEdit()
        xce_object.primary_citation_page_last.setText('')
        xce_object.primary_citation_page_last.setFixedWidth(500)
        xce_object.primary_citation_page_last.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.primary_citation_page_last, 7, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        ## citation authors
        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        xce_object.set_primary_citation_authors = QtGui.QCheckBox('same as deposition authors')
        xce_object.layout_funcs.add_checkbox(xce_object, xce_object.set_primary_citation_authors,
                                       'xce_object.set_primary_citation_as_structure_authors')
        grid.addWidget(xce_object.set_primary_citation_authors, 0, 0)

        xce_object.primary_citation_author_name_List = []

        for column in range(0, 2):
            for row in range(1, 15):
                primary_citation_author_name = QtGui.QLineEdit()
                primary_citation_author_name.setText('')
                primary_citation_author_name.setFixedWidth(300)
                grid.addWidget(primary_citation_author_name, row, column)
                xce_object.primary_citation_author_name_List.append(primary_citation_author_name)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Citation'][1].addLayout(vb)

        ## molecule info
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('Entity 1'), 1, 0)

        grid.addWidget(QtGui.QLabel('Molecule Name'), 2, 0)
        xce_object.molecule_name = QtGui.QLineEdit()
        xce_object.molecule_name.setText('')
        xce_object.molecule_name.setFixedWidth(300)
        xce_object.molecule_name.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.molecule_name, 2, 1)
        grid.addWidget(QtGui.QLabel('(e.g. RNA Hammerhead Ribozyme)'), 2, 2)

        grid.addWidget(QtGui.QLabel('Fragment Name'), 3, 0)
        xce_object.fragment_name_one = QtGui.QLineEdit()
        xce_object.fragment_name_one.setText('')
        xce_object.fragment_name_one.setFixedWidth(300)
        xce_object.fragment_name_one.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.fragment_name_one, 3, 1)
        grid.addWidget(QtGui.QLabel('(e.g. ligand binding domain, hairpin)'), 3, 2)

        grid.addWidget(QtGui.QLabel('Specific Mutation'), 4, 0)
        xce_object.fragment_name_one_specific_mutation = QtGui.QLineEdit()
        xce_object.fragment_name_one_specific_mutation.setText('')
        xce_object.fragment_name_one_specific_mutation.setFixedWidth(300)
        xce_object.fragment_name_one_specific_mutation.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.fragment_name_one_specific_mutation, 4, 1)
        grid.addWidget(QtGui.QLabel('(e.g. C280S)'), 4, 2)

        grid.addWidget(QtGui.QLabel('Enzyme Comission Number'), 5, 0)
        xce_object.fragment_name_one_enzyme_comission_number = QtGui.QLineEdit()
        xce_object.fragment_name_one_enzyme_comission_number.setText('')
        xce_object.fragment_name_one_enzyme_comission_number.setFixedWidth(300)
        xce_object.fragment_name_one_enzyme_comission_number.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.fragment_name_one_enzyme_comission_number, 5, 1)
        grid.addWidget(QtGui.QLabel('(if known: e.g. 2.7.7.7)'), 5, 2)

        grid.addWidget(QtGui.QLabel('Genetically Manipulated Source'), 6, 0)

        grid.addWidget(QtGui.QLabel('Source organism scientific name'), 7, 0)

        xce_object.Source_organism_scientific_name = QtGui.QComboBox()
        taxonomy_dict = XChemMain.NCBI_taxonomy_ID()
        for item in taxonomy_dict:
            xce_object.Source_organism_scientific_name.addItem(taxonomy_dict[item])
        grid.addWidget(xce_object.Source_organism_scientific_name, 7, 1)

        grid.addWidget(QtGui.QLabel('Source organism gene'), 8, 0)
        xce_object.Source_organism_gene = QtGui.QLineEdit()
        xce_object.Source_organism_gene.setText('')
        xce_object.Source_organism_gene.setFixedWidth(300)
        grid.addWidget(xce_object.Source_organism_gene, 8, 1)
        grid.addWidget(QtGui.QLabel('(e.g. RPOD, ALKA...)'), 8, 2)

        grid.addWidget(QtGui.QLabel('Source organism strain'), 9, 0)
        xce_object.Source_organism_strain = QtGui.QLineEdit()
        xce_object.Source_organism_strain.setText('')
        xce_object.Source_organism_strain.setFixedWidth(300)
        xce_object.Source_organism_strain.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.Source_organism_strain, 9, 1)
        grid.addWidget(QtGui.QLabel('(e.g. BH10 ISOLATE, K-12...)'), 9, 2)

        grid.addWidget(QtGui.QLabel('Expression system scientific name'), 10, 0)

        xce_object.Expression_system_scientific_name = QtGui.QComboBox()
        for item in taxonomy_dict:
            xce_object.Expression_system_scientific_name.addItem(taxonomy_dict[item])
        grid.addWidget(xce_object.Expression_system_scientific_name, 10, 1)

        grid.addWidget(QtGui.QLabel('Expression system strain'), 11, 0)
        xce_object.Expression_system_strain = QtGui.QLineEdit()
        xce_object.Expression_system_strain.setText('')
        xce_object.Expression_system_strain.setFixedWidth(300)
        xce_object.Expression_system_strain.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.Expression_system_strain, 11, 1)
        grid.addWidget(QtGui.QLabel('(e.g. BL21(DE3))'), 11, 2)

        grid.addWidget(QtGui.QLabel('Expression system vector type'), 12, 0)
        xce_object.Expression_system_vector_type = QtGui.QLineEdit()
        xce_object.Expression_system_vector_type.setText('')
        xce_object.Expression_system_vector_type.setFixedWidth(300)
        xce_object.Expression_system_vector_type.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.Expression_system_vector_type, 12, 1)
        grid.addWidget(QtGui.QLabel('(e.g. plasmid)'), 12, 2)

        grid.addWidget(QtGui.QLabel('Expression_system_plasmid_name'), 13, 0)
        xce_object.Expression_system_plasmid_name = QtGui.QLineEdit()
        xce_object.Expression_system_plasmid_name.setText('')
        xce_object.Expression_system_plasmid_name.setFixedWidth(300)
        xce_object.Expression_system_plasmid_name.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.Expression_system_plasmid_name, 13, 1)
        grid.addWidget(QtGui.QLabel('(e.g. pET26)'), 13, 2)

        grid.addWidget(QtGui.QLabel('Manipulated_source_details'), 14, 0)
        xce_object.Manipulated_source_details = QtGui.QLineEdit()
        xce_object.Manipulated_source_details.setText('')
        xce_object.Manipulated_source_details.setFixedWidth(300)
        xce_object.Manipulated_source_details.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(xce_object.Manipulated_source_details, 14, 1)
        grid.addWidget(QtGui.QLabel('(any other relevant information)'), 14, 2)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Molecule'][1].addLayout(vb)

        ## misc
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('Keywords'), 1, 0)
        xce_object.structure_keywords = QtGui.QLineEdit()
        xce_object.structure_keywords.setText('SGC - Diamond I04-1 fragment screening, PanDDA, XChemExplorer')
        xce_object.structure_keywords.setFixedWidth(300)
        grid.addWidget(xce_object.structure_keywords, 1, 1)
        grid.addWidget(QtGui.QLabel('(e.g. beta barrel, protein-DNA complex)'), 1, 2)

        grid.addWidget(QtGui.QLabel('Biological Assembly'), 2, 0)
        xce_object.biological_assembly_chain_number = QtGui.QLineEdit()
        xce_object.biological_assembly_chain_number.setText('')
        xce_object.biological_assembly_chain_number.setFixedWidth(300)
        grid.addWidget(xce_object.biological_assembly_chain_number, 2, 1)
        grid.addWidget(QtGui.QLabel('(e.g.  1 for monomer, 2 for dimer ..)'), 2, 2)

        grid.addWidget(QtGui.QLabel('Sequence UNIPROT ID'), 3, 0)
        xce_object.molecule_one_letter_sequence_uniprot_id = QtGui.QLineEdit()
        xce_object.molecule_one_letter_sequence_uniprot_id.setText('')
        xce_object.molecule_one_letter_sequence_uniprot_id.setFixedWidth(300)
        grid.addWidget(xce_object.molecule_one_letter_sequence_uniprot_id, 3, 1)
        grid.addWidget(QtGui.QLabel('(e.g.  Q6B0I6)'), 3, 2)

        grid.addWidget(QtGui.QLabel('Sequence'), 4, 0)
        xce_object.molecule_one_letter_sequence = QtGui.QTextEdit()
        xce_object.molecule_one_letter_sequence.setText('')
        xce_object.molecule_one_letter_sequence.setFixedWidth(300)
        grid.addWidget(xce_object.molecule_one_letter_sequence, 4, 1, 7, 2)

        grid.addWidget(QtGui.QLabel('Structural Genomic (optional)'), 8, 0)

        grid.addWidget(QtGui.QLabel('Project Name'), 9, 0)
        xce_object.SG_project_name = QtGui.QLineEdit()
        xce_object.SG_project_name.setText('')
        xce_object.SG_project_name.setFixedWidth(300)
        grid.addWidget(xce_object.SG_project_name, 9, 1)
        grid.addWidget(QtGui.QLabel('(e.g. PSI, Protein Structure Initiative)'), 9, 2)

        grid.addWidget(QtGui.QLabel('Full Name'), 10, 0)
        xce_object.full_name_of_SG_center = QtGui.QLineEdit()
        xce_object.full_name_of_SG_center.setText('')
        xce_object.full_name_of_SG_center.setFixedWidth(300)
        grid.addWidget(xce_object.full_name_of_SG_center, 10, 1)
        grid.addWidget(QtGui.QLabel('(e.g. Berkeley Structural Genomic Center)'), 10, 2)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Misc'][1].addLayout(vb)

        ## methods
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('Crystallization'), 1, 0)

        grid.addWidget(QtGui.QLabel('Method'), 2, 0)

        xce_object.crystallization_method = QtGui.QComboBox()
        for item in XChemMain.crystal_growth_methods(): xce_object.crystallization_method.addItem(item)
        grid.addWidget(xce_object.crystallization_method, 2, 1)

        grid.addWidget(QtGui.QLabel('pH'), 3, 0)
        xce_object.crystallization_pH = QtGui.QLineEdit()
        xce_object.crystallization_pH.setText('')
        xce_object.crystallization_pH.setFixedWidth(300)
        grid.addWidget(xce_object.crystallization_pH, 3, 1)
        grid.addWidget(QtGui.QLabel('(e.g. 7.5 ...)'), 3, 2)

        grid.addWidget(QtGui.QLabel('Temperature'), 4, 0)
        xce_object.crystallization_temperature = QtGui.QLineEdit()
        xce_object.crystallization_temperature.setText('')
        xce_object.crystallization_temperature.setFixedWidth(300)
        grid.addWidget(xce_object.crystallization_temperature, 4, 1)
        grid.addWidget(QtGui.QLabel('(e.g. 298) (in Kelvin)'), 4, 2)

        grid.addWidget(QtGui.QLabel('Condition'), 5, 0)
        xce_object.crystallization_details = QtGui.QLineEdit()
        xce_object.crystallization_details.setText('')
        xce_object.crystallization_details.setFixedWidth(300)
        grid.addWidget(xce_object.crystallization_details, 5, 1)
        grid.addWidget(QtGui.QLabel('(e.g. PEG 4000, NaCl etc.)'), 5, 2)

        grid.addWidget(QtGui.QLabel('Diffraction Experiment'), 6, 0)
        note = ('Note: this information will only be used if it is\n'
                'not already available in the mainTable!\n'
                'Ignore if data were collected at DLS')
        grid.addWidget(QtGui.QLabel(note), 7, 0)

        grid.addWidget(QtGui.QLabel('Source'), 8, 0)

        xce_object.radiation_source = QtGui.QComboBox()
        for item in XChemMain.radiationSource(): xce_object.radiation_source.addItem(item)
        grid.addWidget(xce_object.radiation_source, 8, 1)

        grid.addWidget(QtGui.QLabel('Source Type'), 9, 0)

        xce_object.radiation_source_type = QtGui.QComboBox()
        for item in XChemMain.wwBeamlines(): xce_object.radiation_source_type.addItem(item)
        grid.addWidget(xce_object.radiation_source_type, 9, 1)

        grid.addWidget(QtGui.QLabel('Wavelength'), 10, 0)
        xce_object.radiation_wavelengths = QtGui.QLineEdit()
        xce_object.radiation_wavelengths.setText('')
        xce_object.radiation_wavelengths.setFixedWidth(300)
        grid.addWidget(xce_object.radiation_wavelengths, 10, 1)
        grid.addWidget(QtGui.QLabel('(e.g. 1.502)'), 10, 2)

        grid.addWidget(QtGui.QLabel('Detector'), 11, 0)

        xce_object.radiation_detector = QtGui.QComboBox()
        for item in XChemMain.detector(): xce_object.radiation_detector.addItem(item)
        grid.addWidget(xce_object.radiation_detector, 11, 1)

        grid.addWidget(QtGui.QLabel('Detector Type'), 12, 0)

        xce_object.radiation_detector_type = QtGui.QComboBox()
        for item in XChemMain.detectorType(): xce_object.radiation_detector_type.addItem(item)
        grid.addWidget(xce_object.radiation_detector_type, 12, 1)

        grid.addWidget(QtGui.QLabel('Date'), 13, 0)
        xce_object.data_collection_date = QtGui.QLineEdit()
        xce_object.data_collection_date.setText('')
        xce_object.data_collection_date.setFixedWidth(300)
        grid.addWidget(xce_object.data_collection_date, 13, 1)
        grid.addWidget(QtGui.QLabel('(e.g. 2004-01-07)'), 13, 2)

        grid.addWidget(QtGui.QLabel('Temperature'), 14, 0)
        xce_object.data_collection_temperature = QtGui.QLineEdit()
        xce_object.data_collection_temperature.setText('')
        xce_object.data_collection_temperature.setFixedWidth(300)
        grid.addWidget(xce_object.data_collection_temperature, 14, 1)
        grid.addWidget(QtGui.QLabel('(e.g. 100) (in Kelvin)'), 14, 2)

        grid.addWidget(QtGui.QLabel('Protocol'), 15, 0)
        xce_object.data_collection_protocol = QtGui.QLineEdit()
        xce_object.data_collection_protocol.setText('SINGLE WAVELENGTH')
        xce_object.data_collection_protocol.setFixedWidth(300)
        grid.addWidget(xce_object.data_collection_protocol, 15, 1)
        grid.addWidget(QtGui.QLabel('(e.g. SINGLE WAVELENGTH, MAD, ...)'), 15, 2)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Methods'][1].addLayout(vb)

        ## software
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('PDB starting model'), 1, 0)
        xce_object.pdbx_starting_model = QtGui.QLineEdit()
        xce_object.pdbx_starting_model.setText('')
        xce_object.pdbx_starting_model.setFixedWidth(300)
        grid.addWidget(xce_object.pdbx_starting_model, 1, 1)
        grid.addWidget(QtGui.QLabel('(e.g. 7.5 ...)'), 1, 2)

        grid.addWidget(QtGui.QLabel('Data reduction'), 2, 0)
        xce_object.data_integration_software = QtGui.QComboBox()
        for item in XChemMain.data_integration_software(): xce_object.data_integration_software.addItem(item)
        grid.addWidget(xce_object.data_integration_software, 2, 1)

        grid.addWidget(QtGui.QLabel('Phasing'), 3, 0)
        xce_object.phasing_software = QtGui.QComboBox()
        for item in XChemMain.phasing_software(): xce_object.phasing_software.addItem(item)
        grid.addWidget(xce_object.phasing_software, 3, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)
        vb.addStretch(1)

        deposit_tab_dict['Software'][1].addLayout(vb)

        vbox.addWidget(deposit_tab_widget)

        hbox = QtGui.QHBoxLayout()
        button = QtGui.QPushButton('Load\nFile')
        button.clicked.connect(xce_object.load_deposit_config_file)
        hbox.addWidget(button)
        button = QtGui.QPushButton('Save\nFile')
        button.clicked.connect(xce_object.save_deposit_config_file)
        hbox.addWidget(button)
        button = QtGui.QPushButton('Load from\nDatabase')
        button.clicked.connect(xce_object.load_deposit_from_database)
        button.setEnabled(False)
        hbox.addWidget(button)
        button = QtGui.QPushButton('Save to\nDatabase')
        button.clicked.connect(xce_object.save_deposit_to_database)
        hbox.addWidget(button)

        vbox.addLayout(hbox)

        xce_object.messageBoxLayout.addLayout(vbox, 0, 0)
        xce_object.messageBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)


class CheckAutoProcessing():
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def query(self, xce_object):
        clearWidgetsAndLayout(xce_object.messageBoxLayout)
#        print 'nn',xce_object.messageBoxLayout.count()
#        for i in range(xce_object.messageBoxLayout.count()):    print i
        warning = ('*** WARNING ***\n'
                   'You did not select a target!\n'
                   'In this case we will only parse the project directory!\n'
                   'Please note that this option is usually only useful in case you reprocessed your data.\n'
                   'Do you want to continue?')
        xce_object.messageBox.setText(warning)
        xce_object.messageBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
        xce_object.messageBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)

        reply = xce_object.messageBox.exec_();
        if reply == 0:
            start_thread = True
        else:
            start_thread = False
        return start_thread