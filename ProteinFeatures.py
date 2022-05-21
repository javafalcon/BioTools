from alphafold import AlphaFoldDB, PDBFileParser
import numpy as np
from sklearn.preprocessing import MinMaxScaler
def DSSPFeature(parser):
    """
    DSSP codes for secondary structure:
    H: Alpha helix; B: Isolated beta-bridge residue; E: Strand; G: 3-10 helix
    I: Pi helix; T: Turn; S: Bend; -: None
    """
    dssp_code={}
    for i,s in enumerate(['H','B','E','G','I','T','S','-']):
        x = [0 for _ in range(8)]
        x[i] = 1
        dssp_code[s] = x

    dsspfeature = []
    # ss = parser.secondaryStructure()
    ss = parser.dsspData()
    for i in range(len(ss)):
        t =[]
        t = t + dssp_code[ss[i][0]]
        for a in ss[i][1:]:
            t.append(a)
        dsspfeature.append(t)
    scaler = MinMaxScaler()
    dsspfeature = np.array(dsspfeature)
    dsspfeature = scaler.fit_transform(dsspfeature)
    # if len(seq) > len(parser.seq): # 匹对的序列比原始序列短，
    #     n = seq.index(parser.seq)
    #     if n > 0: # 原始序列左端被裁剪
    #         for _ in range(n):
    #             dsspfeature.insert(0, dssp_code['-'])
    #     else:# 原始序列右端被裁剪
    #         for _ in range(n):
    #             dsspfeature.append(dssp_code['-'])

    return dsspfeature

def OneHot(seqs):
    text = 'ACDEFGHIKLMNPQRSTVWY'
    aminoacid_idx = {text[i]:i for i in range(20)}
    onehot = []
    for seq in seqs:
        v = [[0 for _ in range(21)] for _ in range(len(seq))]
        for i,a in enumerate(seq):
            k = aminoacid_idx.get(a,20)
            v[i][k] = 1
        onehot.append(v)
    return onehot
if __name__ == "__main__":
    # seqs= ['HMINKKSLLQNLLSKCKTTFQQSFTNANITLKDEKWLKNVRTAYFVCDHDGSVELAYLPNVLPKELVEEFTEKFESIQTGRKKDTGYSGILDNSMPFNYVTADLSQELGQYLSEIVNPQINYYISKLLTCVSSRTINYLVSLNDSYYALNNCLYPSTAFNSLKPSNDGHRIRKPHKDNLDITPSSLFYFGNFQNTEGYLELTDKNCKVFVQPGDVLFFKGNEYKHVVANITSGWRIGLVYFAHKGSKTKPYYEDTQKNSLKIHKETK']
    # alphafold = AlphaFoldDB()
    # pdbfiles = alphafold.getHomoProtPDBFiles(seqs[0])
    # pdbfile, start, end = pdbfiles[0]
    pdbfile = r'E:\Repoes\pytorch\BioTools\protpdb\AF-B1AXP6-F1-model_v2.pdb'
    pdbparser = PDBFileParser(pdbfile)
    dsspfeature = DSSPFeature(pdbparser)