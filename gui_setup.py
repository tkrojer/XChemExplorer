class tables():
    def __init__(self):
        pass

    def setup(self, object):
        # Table column settings

        # functions that use tables.overview_datasource_table_columns:
        #
        # - select_datasource_columns_to_display() - dropdown in datasource top menu (select columns to show)
        # - populate_data_source_table()           - appears to be completely unused, so commented out
        # - populate_and_update_datasource_table() - used within select_datasource_columns_to_display and
        #                                            update_all_tables()

        object.overview_datasource_table_columns = ['Sample ID',
                                                  'Compound ID',
                                                  'Smiles',
                                                  'Visit',
                                                  'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                                  'Refinement\nRfree',
                                                  'Data Collection\nDate',
                                                  'Puck',
                                                  'PuckPosition',
                                                  'Ligand\nConfidence']

        # functions that use tables.datasets_summary_table_columns:
        #
        # - populate_datasets_summary_table() - appears in create_widgets_for_autoprocessing_results_only()
        # - user_update_selected_autoproc_datasets_summary_table()
        #                                            - appears in create_widgets_for_autoprocessing_results_only()

        object.datasets_summary_table_columns = ['Sample ID',
                                               'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                               'DataProcessing\nSpaceGroup',
                                               'DataProcessing\nRfree',
                                               'SoakDB\nBarcode',
                                               'GDA\nBarcode',
                                               'Rmerge\nLow',
                                               'auto-assigned',
                                               'DataCollection\nOutcome',
                                               'img1',
                                               'img2',
                                               'img3',
                                               'img4',
                                               'img5',
                                               'Show\nDetails',
                                               'Show Diffraction\nImage'
                                               ]

        # functions that use tables.datasets_reprocess_columns:
        #
        # - update_datasets_reprocess_table() - appears in search_for_datasets()

        object.datasets_reprocess_columns = ['Dataset ID',
                                           'Sample ID',
                                           'Run\nxia2',
                                           'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                           'Rmerge\nLow',
                                           'Dimple\nRfree',
                                           'DataProcessing\nSpaceGroup',
                                           'DataProcessing\nUnitCell',
                                           'DataProcessing\nStatus']

        # functions that use tables.maps_table_columns:
        #
        # - update_datasets_reprocess_table() - appears in create_maps_table()

        object.maps_table_columns = ['Sample ID',
                                   'Select',
                                   'Compound ID',
                                   'Smiles',
                                   'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                   'Dimple\nRcryst',
                                   'Dimple\nRfree',
                                   'DataProcessing\nSpaceGroup',
                                   'Reference\nSpaceGroup',
                                   'Difference\nUC Volume (%)',
                                   'Reference File',
                                   'DataProcessing\nUnitCell',
                                   'Dimple\nStatus',
                                   'Compound\nStatus',
                                   'LastUpdated']

        # functions that use tables.pandda_table_columns:
        #
        # - populate_pandda_analyse_input_table() - appears in update_all_tables()

        object.pandda_table_columns = ['Sample ID',
                                     'Refinement\nSpace Group',
                                     'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                     'Dimple\nRcryst',
                                     'Dimple\nRfree',
                                     'Crystal Form\nName', ]

        # functions that use tables.refinement_table_columns:
        #
        # - populate_and_update_refinement_table() - appears in update_all_tables

        object.refinement_table_columns = ['Sample ID',
                                         'Compound ID',
                                         'Refinement\nSpace Group',
                                         'Refinement\nResolution',
                                         'Refinement\nRcryst',
                                         'Refinement\nRfree',
                                         'Refinement\nOutcome',
                                         'PanDDA site details',
                                         'Refinement\nStatus']