import re
abbr_AA_dict = {'GLY':'G', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I',
                'PRO':'P', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 'SER':'S',
                'THR':'T', 'CYS':'C', 'MET':'M', 'ASN':'N', 'GLN':'Q',
                'ASP':'D', 'GLU':'E', 'LYS':'K', 'ARG':'R', 'HIS':'H'
                }
class PDBFileParser():
    def __init__(self, pdbfile):
        self.pdbfile = pdbfile
    def get_sequence(self, chain):
        seq = ""
        with open(self.pdbfile, 'r', encoding='utf-8') as fr:
            early_break = False
            for line in fr:
                if line.startswith('SEQRES'):
                    ls = line.replace("\n", "")
                    ls = ls.strip()
                    ls = re.sub(r'\s+', ' ', ls)
                    ls = ls.split(" ")
                    if ls[2] == chain:
                        early_break = True
                        for aa in ls[4:]:
                            seq += abbr_AA_dict[aa]
                elif early_break:
                    break

        return seq

    def get_atoms(self, chain):
        atoms, atom = [], []
        amino_acid = ''
        with open(self.pdbfile, 'r', encoding='utf-8') as fr:
            early_break = False
            firstAtom = True
            for line in fr:
                if line.startswith('ATOM'):
                    ls = line.replace("\n", "")
                    ls = ls.strip()
                    ls = re.sub(r'\s+', ' ', ls)
                    ls = ls.split(" ")
                    if ls[4] == chain:
                        early_break = True

                        if firstAtom:
                            amino_acid = ls[3]
                            firstAtom = False
                            atom.append([float(v) for v in ls[6:9]])
                        else:
                            if amino_acid == ls[3]:
                                atom.append([float(v) for v in ls[6:9]])
                            else:
                                atoms.append(atom)
                                atom = []
                                amino_acid = ls[3]
                                atom.append([float(v) for v in ls[6:9]])

                elif early_break:
                    atoms.append(atom)
                    break

            return atoms

if __name__ == "__main__":
    parser = PDBFileParser(r'e:\downloads\4ZCF.pdb')
    atoms = parser.get_atoms('C')