import re
import libtbx.phil
import sys

from iotbx.pdb import hierarchy
from scitbx.array_family import flex
############################################################################

PROGRAM = 'wrapper_merge_confs'

DESCRIPTION = """
    A tool to alter the merged model from (giant.merge_conformations) to 
     edit occupancy of residues that are not in occupancy groups, and thus 
     will not be refined. To be used in xce pandda refinement pipeline, 
     before giant.quick_refine
"""

############################################################################
blank_arg_prepend = {'.pdb': 'input.pdb='}

master_phil = libtbx.phil.parse("""
input {
    refmac_params_file = None
        .help = 'The params file that contains the occupancy groups to be compared to' 
        .type = str
    pdb = 'multi-state-model.pdb'
        .help = 'Input pdb file that conatins residues with occupancy over 1.0'
        .type = str
}
output {
    pdb = 'multi-state-model.pdb'
        .help = 'output pdb file'
        .type = str
}
""", process_includes=True)

def over_one_occupancies(pdb):
    """ Find all residues that have occupancy summing over one in their altlocs.
    
    Code is altered from xce-dev-tk:XChemUtils.pdbtools.check_occupancies()
    
    """
    error_text = ''
    warning_text = ''
    residue_dict = {}
    all_over_one_occ = set()
    for line in open(pdb):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atomname = str(line[12:16]).replace(' ', '')
            resname = str(line[17:20]).replace(' ', '')
            chainID = str(line[21:22]).replace(' ', '')
            resseq = str(line[22:26]).replace(' ', '')
            alt_loc = str(line[16:17]).replace(' ', '')
            occupancy = str(line[56:60]).replace(' ', '')
            if resname + '-' + chainID + '-' + resseq not in residue_dict:
                residue_dict[resname + '-' + chainID + '-' + resseq] = []
            if alt_loc == '':
                alt_loc = '0'
            residue_dict[resname + '-' + chainID + '-' + resseq].append(
                [alt_loc, occupancy, atomname, resname, chainID, resseq])
    for item in residue_dict:
        alt_loc_dict = {}
        for atom in residue_dict[item]:
            alt_loc = atom[0]
            occupancy = atom[1]
            atomname = atom[2]
            resname = atom[3]
            chainID = atom[4]
            resseq = atom[5]
            if alt_loc not in alt_loc_dict:
                alt_loc_dict[alt_loc] = []
            alt_loc_dict[alt_loc].append([occupancy, atomname, resname, chainID, resseq])
            if float(occupancy) > 1:
                error_text += 'ERROR: %s %s %s (%s) %s: occupancy is %s\n' % (
                chainID, resname, resseq, alt_loc, atomname, occupancy)
        occupancySumList = []
        for alt_loc in alt_loc_dict:
            occupancySum = 0.0
            nAtom = float(len(alt_loc_dict[alt_loc]))
            for n, atom in enumerate(alt_loc_dict[alt_loc]):
                occupancy = atom[0]
                occupancySum += float(occupancy)
                atomname = atom[1]
                resname = atom[2]
                chainID = atom[3]
                resseq = atom[4]
                if n == 0:
                    occStart = occupancy
                else:
                    if occupancy != occStart:
                        warning_text += '%s %s %s (%s) %s: occupancy differs for altLoc -> %s\n' % (
                        chainID, resname, resseq, alt_loc, atomname, occupancy)
            occupancySumList.append(occupancySum / nAtom)
        occAdd = 0.0
        for entry in occupancySumList:
            occAdd += entry
        if occAdd > 1:
            error_text += 'ERROR: ' + item + ' -> summarised occupanies of alternative conformations are > 1.0 (' + str(
                occupancySumList) + ')\n'
            all_over_one_occ.update([(chainID,resseq)])

    return all_over_one_occ

def get_residues_in_occ_groups(refmac_params_file):
    """ Read the refmac params file and return residue/chain in an occupancy group """

    # open refmac params file
    occupancy_lines = [line for line in open(refmac_params_file) if line.startswith('occupancy')]

    all_occ_group_residues = set()
    for occ_line in occupancy_lines:
        # Get the residue in the occupancy group lines
        try:
            resi = re.search('resi(.+?)alte',occ_line).group(1)
            resi = "".join(resi.split())
        except AttributeError:
            ### can't find regular expression
            resi = ''

        # Get the chain in occupancy group lines
        try:
            chain = re.search('chain(.+?)resi',occ_line).group(1)
            chain = "".join(chain.split())
        except AttributeError:
            ### Can't find reqular expression
            chain = ''

        if resi and chain != "":
            all_occ_group_residues.update([(chain,resi)])

    return all_occ_group_residues

def alter_pdb_occupancy_not_in_occ_group(pdb,all_over_one_occ,all_occ_group_residues,output_pdb):
    """ Generate PDB file with occupancies reduced to 1/altlocs for residues not refined
        
        Take input pdb file path, output pdb file path, and two sets of 
        residue/chain pairs,and alter occupancy of residues in difference 
        between the sets """

    if all_over_one_occ - all_occ_group_residues:
        pdb_in = hierarchy.input(file_name = pdb)

        sel_cache = pdb_in.hierarchy.atom_selection_cache()
        for residue in all_over_one_occ - all_occ_group_residues:
            altlocs = []
            over_one_not_occ_group = sel_cache.selection("resseq {} and chain {}".format(residue[1],residue[0]))
            residue_hier = pdb_in.hierarchy.select(over_one_not_occ_group)


            for ag in residue_hier.atom_groups():
                altlocs.append(ag.altloc)

            new_occ = 1.0/len(altlocs)

            for ag in pdb_in.hierarchy.atom_groups():
                if ag.parent().resseq == residue[1]:
                    ag.atoms().set_occ(flex.double(ag.atoms().size(), new_occ))

        # Write output pdb file
        pdb_in.hierarchy.write_pdb_file(file_name=output_pdb,
                                        crystal_symmetry=pdb_in.input.crystal_symmetry())

def replace_occupancies(pdb, refmac_params_file, output_pdb):
    """"Generate PDB file with occupancies reduced to 1/altlocs for residues not refined"""

    over_one_occ = over_one_occupancies(pdb)
    all_occ_group_residues = get_residues_in_occ_groups(refmac_params_file)
    alter_pdb_occupancy_not_in_occ_group(pdb, over_one_occ, all_occ_group_residues, output_pdb)

#######################################################################

def run(params):
    replace_occupancies(pdb = params.input.pdb,
                        refmac_params_file = params.input.refmac_params_file,
                        output_pdb = params.output.pdb)
#######################################################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)