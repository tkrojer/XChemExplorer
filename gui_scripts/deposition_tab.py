import sys, os
from PyQt4 import QtGui, QtCore, QtWebKit

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'gui_scripts'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))

import layout
import XChemToolTips
import XChemMain


class DepositionTab():
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
        ################################################################################################################
        #                                                                                                              #
        #                                                DEPOSITION TAB                                                #
        #                                                                                                              #
        ################################################################################################################
        xce_object.deposition_vbox = QtGui.QVBoxLayout()

        scroll = QtGui.QScrollArea()
        xce_object.deposition_vbox.addWidget(scroll)
        scrollContent = QtGui.QWidget(scroll)
        scrollLayout = QtGui.QVBoxLayout(scrollContent)
        scrollContent.setLayout(scrollLayout)

        # deposition page heading
        deposition_page_heading = self.layout_funcs.add_depo_heading('HTML export & ZENODO upload')
        deposition_page_heading.setStyleSheet("font: 30pt Arial Bold")

        # introduction text
        introduction_text = self.layout_funcs.add_depo_text(XChemToolTips.html_summary_introduction())

        # zenodo example image
        zenodo_image = QtGui.QLabel()
        zenodo_pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'html_summary_page.png'))
        zenodo_image.setPixmap(zenodo_pixmap)

        # export html directory heading
        export_html_heading = self.layout_funcs.add_depo_heading('1. Specify HTML export directory in the settings tab')

        # export html directory text
        export_html_text = self.layout_funcs.add_depo_text(XChemToolTips.html_export_directory_background())

        # default for dls
        dls_dir_text = self.layout_funcs.add_depo_text(
            'Note: default for labxchem project at DLS is <labxchem_directory>/processing/html.')

        # user specific directory
        usr_dir_text = self.layout_funcs.add_depo_text('In your case: ' + xce_object.html_export_directory)

        # html prep heading
        html_prep_heading = self.layout_funcs.add_depo_heading("2. Prepare files for HTML export")

        # html prep test
        html_prep_text = self.layout_funcs.add_depo_text(XChemToolTips.html_export_step())

        # html export button
        html_export_button = QtGui.QPushButton('Export to HTML')
        html_export_button.clicked.connect(xce_object.export_to_html)
        html_export_button.setMaximumWidth(200)

        # prepare ICM heading
        ICM_heading = self.layout_funcs.add_depo_heading("3. Prepare ICM files")

        # ICM background text
        ICM_bgr_text = self.layout_funcs.add_depo_text(XChemToolTips.icb_file_background())

        # ICM prep text
        ICM_prep_text = self.layout_funcs.add_depo_text(XChemToolTips.prepare_ICB_files())

        # Open ICM button
        ICM_button = QtGui.QPushButton('Open ICM-pro')
        ICM_button.clicked.connect(xce_object.open_icm)
        ICM_button.setMaximumWidth(200)

        # ICM loading image
        ICM_image = QtGui.QLabel()
        ICM_pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'drag_and_drop_icb_file.png'))
        ICM_image.setPixmap(ICM_pixmap)

        # ICM run image
        ICM_run_image = QtGui.QLabel()
        ICM_run_pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'run_icm_script.png'))
        ICM_run_image.setPixmap(ICM_run_pixmap)

        # zenodo upload heading 1
        zenodo_upload_heading = self.layout_funcs.add_depo_heading("4. Prepare files for ZENODO upload")

        # zenodo upload text 1
        zenodo_upload_text = self.layout_funcs.add_depo_text(
            XChemToolTips.zenodo_upload_start(xce_object.html_export_directory))

        # prepare files button
        prep_files_button = QtGui.QPushButton('prepare files')
        prep_files_button.clicked.connect(xce_object.prepare_files_for_zenodo_upload)
        prep_files_button.setMaximumWidth(200)

        # zenodo upload heading 2
        zenodo_upload_heading2 = self.layout_funcs.add_depo_heading("5. ZENODO")

        # zenodo upload text 2
        zenodo_upload_text2 = self.layout_funcs.add_depo_text(
            XChemToolTips.zenodo_upload_part_one(xce_object.html_export_directory))

        # zenodo upload image
        zenodo_upload_image = QtGui.QLabel()
        zenodo_upload_pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image',
                                                          'new_zenodo_upload.png'))
        zenodo_upload_image.setPixmap(zenodo_upload_pixmap)

        # zenodo upload ID heading
        zenodo_upload_ID_heading = self.layout_funcs.add_depo_heading("6. ZENODO upload ID")

        # zenodo upload ID text
        zenodo_upload_ID_text = self.layout_funcs.add_depo_text(XChemToolTips.zenodo_upload_part_two())

        # zenodo upload ID image
        zenodo_upload_image2 = QtGui.QLabel()
        zenodo_upload_pixmap2 = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image',
                                                           'zenodo_upload_id.png'))
        zenodo_upload_image2.setPixmap(zenodo_upload_pixmap2)

        # zenodo upload ID text 2
        zenodo_upload_ID_text2 = self.layout_funcs.add_depo_text(XChemToolTips.zenodo_upload_part_three())

        # zenodo upload ID entry box
        hbox_zenodo_upload_id = QtGui.QHBoxLayout()
        hbox_zenodo_upload_id.addWidget(QtGui.QLabel('upload ID:'))
        xce_object.zenodo_upload_id_entry = QtGui.QLineEdit()
        xce_object.zenodo_upload_id_entry.setFixedWidth(200)
        hbox_zenodo_upload_id.addWidget(xce_object.zenodo_upload_id_entry)
        hbox_zenodo_upload_id.addStretch(1)

        # update html button
        update_html_button = QtGui.QPushButton('update html files with upload ID')
        update_html_button.clicked.connect(xce_object.update_html_for_zenodo_upload)
        update_html_button.setMaximumWidth(300)

        # zenodo upload html files heading
        upload_html_heading = self.layout_funcs.add_depo_heading("7. ZENODO upload HTML files")

        # zenodo upload html text
        upload_html_text = self.layout_funcs.add_depo_text(XChemToolTips.zenodo_upload_part_four(xce_object.
                                                                                                 html_export_directory))

        deposition_widget_list = [deposition_page_heading, QtGui.QLabel('  '), introduction_text, QtGui.QLabel('  '),
                                  zenodo_image, QtGui.QLabel('  '), export_html_heading, export_html_text, dls_dir_text,
                                  usr_dir_text, QtGui.QLabel('  '), html_prep_heading, html_prep_text,
                                  html_export_button, QtGui.QLabel('  '), ICM_heading, ICM_bgr_text, ICM_prep_text,
                                  ICM_button, QtGui.QLabel('  '), ICM_image, QtGui.QLabel('  '), ICM_run_image,
                                  QtGui.QLabel('  '), zenodo_upload_heading, zenodo_upload_text, prep_files_button,
                                  QtGui.QLabel('  '), zenodo_upload_heading2, zenodo_upload_text2, zenodo_upload_image,
                                  QtGui.QLabel('  '), zenodo_upload_ID_heading, zenodo_upload_ID_text,
                                  zenodo_upload_image2, zenodo_upload_ID_text2]

        deposition_widget_list2 = [update_html_button, QtGui.QLabel('  '),
                                   upload_html_heading, upload_html_text, QtGui.QLabel('  ')]

        self.layout_funcs.add_to_box(scrollLayout, deposition_widget_list)
        scrollLayout.addLayout(hbox_zenodo_upload_id)
        self.layout_funcs.add_to_box(scrollLayout, deposition_widget_list2)

        # container settings
        scrollLayout.addStretch(1)
        scroll.setWidget(scrollContent)


