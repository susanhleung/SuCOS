import argparse, os
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers

def get_FeatureMapScore(small_m, large_m):
    from rdkit.Chem.FeatMaps import FeatMaps
    from rdkit import RDConfig
    fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
    #    keep = ('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic')

    fmParams = {}
    for k in fdef.GetFeatureFamilies():
        fparams = FeatMaps.FeatMapParams()
        fmParams[k] = fparams

    keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
            'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
    featLists = []
    for m in [small_m, large_m]:
        rawFeats = fdef.GetFeaturesForMol(m)
        # filter that list down to only include the ones we're intereted in
        featLists.append([f for f in rawFeats if f.GetFamily() in keep])
    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
    fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
    return fm_score

def main(ref_file, prb_file):
    if type(ref_file) == str:
        if os.path.splitext(ref_file)[-1] == '.sdf':
            reflig = Chem.MolFromMolFile(ref_file, sanitize=True)


    if type(prb_file) == str:
        if os.path.splitext(prb_file)[-1] == '.sdf':
            prb_mols = Chem.SDMolSupplier(prb_file, sanitize=True)
    
    try: reflig
    except NameError:
        raise ValueError("Incorrect file format for ref lig" )
    try: prb_mols
    except NameError:
        raise ValueError("Incorrect file format for prb lig" )

    w = Chem.SDWriter("%s_SuCOS_score.sdf" % os.path.splitext(prb_file)[0])
    for prb_mol in prb_mols:
        ##############################################
        ####### Feature map
        ##############################################
        fm_score = get_FeatureMapScore(reflig, prb_mol)
        ##############################################

        #tversky_ind = rdShapeHelpers.ShapeTverskyIndex(reflig, prb_mol, 1.0, 0.0)
        #SuCOS_score = 0.5*fm_score + 0.5*tversky_ind
        
        protrude_dist = rdShapeHelpers.ShapeProtrudeDist(reflig, prb_mol,
                allowReordering=False)
        SuCOS_score = 0.5*fm_score + 0.5*(1 - protrude_dist)

        print "********************************"
        print "SuCOS score:\t%f" % SuCOS_score
        print "********************************"
        prb_mol.SetProp("SuCOS_score", str(SuCOS_score))
        w.write(prb_mol)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="run SuCOS")
    parser.add_argument('--lig1', help='first ligand')
    parser.add_argument('--lig2', help='second ligand')

    args = parser.parse_args()

    ref_file = args.lig1
    prb_file = args.lig2

    main(ref_file, prb_file)
