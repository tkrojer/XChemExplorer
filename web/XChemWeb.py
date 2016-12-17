# last edited: 16/12/2016, 15:00

import os,sys

def create_ICM_input_file(html_export_directory,database):
    icm_in = (
        '#!/usr/local/bin/icm\n'
        '# Author: Brian Marsden\n'
        'panddaDir="%s"\n' %html_export_directory+
        '\n'
        'set directory panddaDir\n'
        'connect molcart filename="%s"\n' %database+

#        'query molcart "select p.ID,p.CrystalName,p.PANDDA_site_event_index,p.PANDDA_site_confidence,p.CrystalName || '_event'|| p.PANDDA_site_event_index as ModelName,m.CompoundCode,m.CompoundSMILES,p.PANDDA_site_name,p.PANDDA_site_confidence as LigandConfidence,p.RefinementOutcome as ModelStatus,p.PANDDA_site_comment,p.PANDDA_site_x,p.PANDDA_site_y,p.PANDDA_site_z, p.PANDDA_site_spider_plot,m.DataProcessingResolutionHigh,m.DataProcessingSpaceGroup,m.DataProcessingUnitCell,m.RefinementPDB_latest,m.RefinementMTZ_latest,p.PANDDA_site_event_map from panddaTable as p, mainTable as m where p.CrystalName=m.CrystalName and p.PANDDA_site_ligand_placed='True' and (LigandConfidence like '1%' or LigandConfidence like '2%' or LigandConfidence like '3%' or LigandConfidence like '4%') order by p.CrystalName,ModelStatus desc,PANDDA_site_event_index" name="T"\n'


        'query molcart "select p.ID,p.CrystalName,p.PANDDA_site_event_index,p.PANDDA_site_confidence,p.CrystalName || '
        "'_event'|| p.PANDDA_site_event_index as ModelName,m.CompoundCode,m.CompoundSMILES,p.PANDDA_site_name,p.PANDDA_site_confidence "
        "as LigandConfidence,p.RefinementOutcome as ModelStatus,p.PANDDA_site_comment,p.PANDDA_site_x,p.PANDDA_site_y,p.PANDDA_site_z, "
        "p.PANDDA_site_spider_plot,m.DataProcessingResolutionHigh,m.DataProcessingSpaceGroup,m.DataProcessingUnitCell,"
        "m.RefinementPDB_latest,m.RefinementMTZ_latest,p.PANDDA_site_event_map from panddaTable as p, mainTable as m "
        "where p.CrystalName=m.CrystalName and p.PANDDA_site_ligand_placed='True' and (LigandConfidence like '1%' "
        "or LigandConfidence like '2%' or LigandConfidence like '3%' or LigandConfidence like '4%') "
        'order by p.CrystalName,ModelStatus desc,PANDDA_site_event_index" name="T"\n'


        'numToDo=Nof(T)\n'
        '\n'
        'macro generateICB s_eventID s_pdb s_ligandSMILES s_eventMap R_coords\n'
        'l_confirm=no\n'
        'l_info=no\n'
        'l_commands=no\n'
        'l_warn=no\n'
        'print s_eventID s_pdb s_ligandSMILES s_eventMap R_coords\n'
        '# Read pdb\n'
        'read pdb s_pdb\n'
        'if Nof(a_1.?lig)==0 then\n'
        '  print "No ligand!"\n'
        '  delete a_*.*\n'
        '  return\n'
        'endif\n'
        '# Fix ligand topology\n'
        'set bond topology a_1.?lig s_ligandSMILES\n'
        '# Read event map\n'
        'read map s_eventMap\n'
        'contourEDS Name( map )[1] {2.0} {"cyan"} a_1.?lig | Sphere(a_1.?lig !a_1.?lig 7.5) yes yes\n'
        'assign sstructure\n'
        'cool a_ no\n'
        'color background refresh rgb={0,0,0}\n'
        'display xstick Res(Sphere(a_1.?lig !a_1.?lig 7.5))\n'
        'read binary s_icmhome+"shapes" name="star"\n'
        'display star\n'
        'translate star R_coords\n'
        'center star margin=0.0\n'
        'undisplay star\n'
        'center a_1.?lig\n'
        'write png "mapImages/"+s_eventID+"_large.png" window={800,600} GRAPHICS.quality=Max(image graphic)\n'
        'write png "mapImages/"+s_eventID+"_small.png" window={150,150} GRAPHICS.quality=Max(image graphic)\n'
        'delete star\n'
        'writeProject "icbs/"+s_eventID+".icb" no\n'
        'delete maps\n'
        'delete a_*.\n'
        'delete grob\n'
        'endmacro\n'
        '\n'
        'for i=1,numToDo\n'
        '  ligandSMILES=T.CompoundSMILES[i]\n'
        '  # PDB file path\n'
        '  pdbID="pdbs/"+T.ModelName[i]+".pdb"\n'
        '  if (T.PANDDA_site_name[i]!="" & T.PANDDA_site_confidence[i]!="None") then\n'
        '  # Get event ID\n'
        '  eventID=T.ModelName[i]+"_"+T.CompoundCode[i]\n'
        '  # Event map path\n'
        '  eventMap="maps/"+T.ModelName[i]+".ccp4"\n'
        '  # Ligand centre\n'
        '  ligR3={$T.PANDDA_site_x[i],$T.PANDDA_site_y[i],$T.PANDDA_site_z[i]}\n'
        '  generateICB eventID,pdbID,ligandSMILES,eventMap,ligR3\n'
        '  endif\n'
        'endfor\n'
        )

    f=open('%s/dsEvent_sqlite.icm' %html_export_directory,'w')
    f.write(icm_in)
    f.close()