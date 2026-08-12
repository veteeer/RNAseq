import pandas as pd
import argparse
import os
import PAM50_module.pam50_classifier as clf

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Assign called variants')
    parser.add_argument('--input', type=str, help='Path to quant.genes.sf', required=True)
    parser.add_argument('--prefix', type=str, help='Output dir/prefix', required=True)
    parser.add_argument('--refdir', type=str, help='Path to reference folder', required=True)

    args = parser.parse_args()
    input_file = args.input
    prefix = args.prefix
    refdir = args.refdir

    df = pd.read_csv(input_file, sep="\t", index_col=0)
    df.index.name = 'Gene id'

    mapping = pd.read_csv(os.path.join(refdir, 'ensembl_tx2gene.tsv'), sep='\t', index_col=0)
    df = df.join(mapping)

    sample = df[['Gene name','TPM']]
    sample.set_index('Gene name', inplace=True)
    sample.to_csv(prefix + '_tpm.tsv', sep='\t')

    model = clf.Classifier_PAM50()
    subtype = model.predict(sample['TPM'])
    with open(prefix + "_pred_subtype.txt", "w", encoding="utf-8") as file:
        file.write(subtype)

    proba_pred = model.predict_proba(sample['TPM'])
    proba_pred.name = 'Probability'
    proba_pred.to_csv(prefix +'_pred_proba.tsv', sep='\t')

    

    