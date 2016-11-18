# last edited: 18/11/2016, 15:00

import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

import XChemDB


def update_data_source(db_file,crystalString,status):
    db=XChemDB.data_source(db_file)

    print "update mainTable set PANDDAStatus = '%s' where CrystalName in (%s)" %(status,"'"+crystalString.replace(",","','")+"'")
    db.execute_statement("update mainTable set PANDDAStatus = '%s' where CrystalName in (%s)" %(status,"'"+crystalString.replace(",","','")+"'"))


if __name__=='__main__':
    db_file=sys.argv[1]
    crystalString=sys.argv[2]
    status=sys.argv[3]

    update_data_source(db_file,crystalString,status)
