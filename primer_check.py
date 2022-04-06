from optparse import Values
from Bio import SeqIO
from itertools import product, combinations_with_replacement
from scipy import stats
from math import log10
from matplotlib import rcParams
import seaborn as sns, matplotlib.pyplot as plt
import numpy as np, argparse, sys, pandas as pd,  plotly.graph_objects as go

#Arguments
parser = argparse.ArgumentParser(description='Estimate self-pairing tendency for a sequence set.')
parser.add_argument('-f', '--fasta',  metavar='\b', help='Path to fasta file.', default=None)
parser.add_argument('-v', '--html', metavar='\b', help='Path to csv file to create html heatmaps', default=None)
parser.add_argument('-s', '--split-images', metavar='\b', help='Path to a folder where to save splitting result', default=None)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args=parser.parse_args()


if args.fasta:
#Read sequences from fasta/file and 1-hot encode each base
#Extract 2-mers from sequences & 1-hot encode each base
#Map sequence id to the corresponding set of kmers
    seqs = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    encoding_dict = {
        "A":np.asarray([1,0,0,0]),
        "C":np.asarray([0,1,0,0]),
        "G":np.asarray([0,0,1,0]),
        "T":np.asarray([0,0,0,1])
    }

    decoding_dict = {
        "[1, 0, 0, 0]":"A",
        "[0, 1, 0, 0]":"C",
        "[0, 0, 1, 0]":"G",
        "[0, 0, 0, 1]":"T"
    }

    for id in seqs.keys():
        ohe_seq = np.asarray([encoding_dict[nt] for nt in seqs[id].seq])
        if len(ohe_seq)%2 == 1:
            two_mers = np.asarray([ohe_seq[i:i+2] for i in range(0,len(ohe_seq),2)], dtype=object)
        else:
            two_mers = np.asarray([ohe_seq[i:i+2] if i != len(ohe_seq)-1 else ohe_seq[i] for i in range(0,len(ohe_seq),2)], dtype=object)
        seqs[id] = [seqs[id],ohe_seq,two_mers]

    #Generate all pairs of sequence ids (including self) - first sequence is used as ref in kmer sliding
    id_list = list(seqs.keys())
    id_dict = {id:{} for id in id_list}
    if len(id_list) < 99:
        combo_list = list(combinations_with_replacement(id_list, 2))
    else:
        combo_list = list(product(id_list, repeat=2))
        
    #For each sequence combination - slide 2-mers from the split sequence along the reference sequence and compute sums of one-hot encodings for each aligned base (use index sliding)
    def slide_twomer(reference,twomer):
        null_limit = 1.0000000000000001e-34
        twomer_sum = np.asarray([0,0,0,0], dtype='object')
        for i,_ in enumerate(reference):
            if i == 0:
                twomer_sum += twomer[-1]+reference[i]
            elif i == len(reference)-1:
                twomer_sum += 2*twomer[0]+reference[i-1]+twomer[-1]+2*reference[i]
            else:
                twomer_sum += twomer[0]+reference[i]+twomer[-1]+reference[i-1]
            aa_agr = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([2,null_limit,null_limit,null_limit])))
            ac_agr = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([1,1,null_limit,null_limit])))
            ag_agr = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([1,null_limit,1,null_limit])))
            at_agr = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([1,null_limit,null_limit,1])))
            cc_agr = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([null_limit,2,null_limit,null_limit])))
            gc_arg = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([null_limit,1,1,null_limit])))
            ct_agr = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([null_limit,1,null_limit,1])))
            gt_agr = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([null_limit,null_limit,1,1])))
            gg_agr = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([null_limit,null_limit,2,null_limit])))
            tt_agr = np.sum(np.square(twomer_sum/np.amax(twomer_sum) - np.asarray([null_limit,null_limit,null_limit,2])))
        return {'aa':aa_agr,'ac':ac_agr,'ag':ag_agr,'at':at_agr,"cc":cc_agr,'gc':gc_arg,'ct':ct_agr,'gt':gt_agr, 'tt':tt_agr, 'gg':gg_agr}

    #Compute phred scored p values for all sequence combinations
    sqdiff = pd.DataFrame()
    for combo in combo_list:
        ref_id = combo[0]
        twomer_id = combo[1]
        reference = seqs[ref_id][1]
        twomers = seqs[twomer_id][2]
        aggregate_sum = {'aa':0,'ac':0,'ag':0,'at':0,"cc":0,'gc':0,'ct':0,'gt':0, 'tt':0, 'gg':0}

        for twomer in twomers: #Sum for all kmers & aggregate 
            twomer_dict = slide_twomer(reference, twomer)
            for key in aggregate_sum.keys():
                aggregate_sum[key] += twomer_dict[key]
        
        aggregate_sum = {key:[-10*log10(stats.chi2.pdf(value/len(twomers),9))] for key,value in aggregate_sum.items()}
        aggregate_sum["ref_id"] = [ref_id]
        aggregate_sum["split_id"] = [twomer_id]
        print(f"Finished: {aggregate_sum}")
        sqdiff = sqdiff.append(pd.DataFrame.from_dict(aggregate_sum))
    sqdiff.to_csv(f"{str(args.fasta).replace('.fasta','')}_report.csv",header=True,index=False)


    # Visualize p values for all sequences in the original file
    if len(sqdiff["ref_id"].unique()) > 99:
        pairs = ['aa', 'ac', 'ag', 'at', "cc", 'gc', 'ct', 'gt',  'tt','gg']
        for id in sqdiff["ref_id"].unique():
            df = sqdiff.loc[sqdiff["ref_id"] == id]
            df["name_pairs"] = df["ref_id"]+"-"+df["split_id"]
            df.set_index(["ref_id",'split_id'], inplace=True)
            df.reset_index(drop=True,inplace=True)
            df.set_index("name_pairs", inplace=True)
            fig = go.Figure()
            fig.add_trace(go.Surface(x=df.columns.tolist(), y=df.index.tolist(), z=df.values.tolist(), colorscale="Viridis", showscale = False))

            # Update plot sizing
            fig.update_layout(width=800,height=900,autosize=False,margin=dict(t=0, b=0, l=0, r=0),template="plotly_white",)

            # Update 3D scene options
            fig.update_scenes(aspectratio=dict(x=1, y=1, z=0.7),aspectmode="manual")

            # Add dropdown
            fig.update_layout(updatemenus=[dict(buttons=list([dict(args=["type", "surface"],label="3D Surface",method="restyle"),dict(args=["type", "heatmap"],label="Heatmap",method="restyle")]),
                    direction="down", pad={"r": 10, "t": 10}, showactive=True, x=0.1, xanchor="left", y=1.1,yanchor="top"),])
            fig.write_html(f'heatmap_{id}.html',full_html=False,include_plotlyjs='cdn')
    else:
        sqdiff["name_pairs"] = sqdiff["ref_id"]+"-"+sqdiff["split_id"]
        sqdiff.set_index(["ref_id",'split_id'], inplace=True)
        sqdiff.reset_index(drop=True,inplace=True)
        sqdiff.set_index("name_pairs", inplace=True)
        fig = go.Figure()
        fig.add_trace(go.Surface(x=sqdiff.columns.tolist(), y=sqdiff.index.tolist(), z=sqdiff.values.tolist(), colorscale="Viridis", showscale = False))

        # Update plot sizing
        fig.update_layout(width=800,height=900,autosize=False,margin=dict(t=0, b=0, l=0, r=0),template="plotly_white",)

        # Update 3D scene options
        fig.update_scenes(aspectratio=dict(x=1, y=1, z=0.7),aspectmode="manual")

        # Add dropdown
        fig.update_layout(updatemenus=[dict(buttons=list([dict(args=["type", "surface"],label="3D Surface",method="restyle"),dict(args=["type", "heatmap"],label="Heatmap",method="restyle")]),
                direction="down", pad={"r": 10, "t": 10}, showactive=True, x=0.1, xanchor="left", y=1.1,yanchor="top"),])
        fig.write_html(f'heatmap.html',full_html=False,include_plotlyjs='cdn')
        fig.show()
# create figures
    if args.split_images:
        sqdiff = pd.read_csv(f"{str(args.fasta).replace('.fasta','')}_report.csv")
        rcParams['figure.figsize'] = 30,0.46*len(sqdiff['ref_id'].unique())
        for id in sqdiff["ref_id"].unique():
            df = sqdiff.loc[sqdiff["ref_id"] == id]
            df = df.set_index(["ref_id",'split_id'])
            heatmap = sns.heatmap(data=df, annot=True, cmap="afmhot", cbar=False)
            plt.xticks(rotation=90)
            plt.yticks(rotation=0)
            plt.savefig(f"{args.split_images}/heatmap_{id}.svg")
            plt.clf()

if args.html:
    sqdiff = pd.read_csv(args.html)
    if len(sqdiff["ref_id"].unique()) > 99:
        pairs = ['aa', 'ac', 'ag', 'at', "cc", 'gc', 'ct', 'gt',  'tt','gg']
        for id in sqdiff["ref_id"].unique():
            df = sqdiff.loc[sqdiff["ref_id"] == id]
            df["name_pairs"] = df["ref_id"]+"-"+sqdiff["split_id"]
            df.set_index(["ref_id",'split_id'], inplace=True)
            df.reset_index(drop=True,inplace=True)
            df.set_index("name_pairs", inplace=True)
            fig = go.Figure()
            fig.add_trace(go.Surface(x=df.columns.tolist(), y=df.index.tolist(), z=df.values.tolist(), colorscale="Viridis"))

            # Update plot sizing
            fig.update_layout(width=800,height=900,autosize=False,margin=dict(t=0, b=0, l=0, r=0),template="plotly_white",)

            # Update 3D scene options
            fig.update_scenes(aspectratio=dict(x=1, y=1, z=0.7),aspectmode="manual")

            # Add dropdown
            fig.update_layout(updatemenus=[dict(buttons=list([dict(args=["type", "surface"],label="3D Surface",method="restyle"),dict(args=["type", "heatmap"],label="Heatmap",method="restyle")]),
                    direction="down", pad={"r": 10, "t": 10}, showactive=True, x=0.1, xanchor="left", y=1.1,yanchor="top"),])
            fig.write_html(f'heatmap_{id}.html',full_html=False,include_plotlyjs='cdn')
    else:
        sqdiff["name_pairs"] = sqdiff["ref_id"]+"-"+sqdiff["split_id"]
        sqdiff.set_index(["ref_id",'split_id'], inplace=True)
        sqdiff.reset_index(drop=True,inplace=True)
        sqdiff.set_index("name_pairs", inplace=True)
        fig = go.Figure()
        fig.add_trace(go.Surface(x=sqdiff.columns.tolist(), y=sqdiff.index.tolist(), z=sqdiff.values.tolist(), colorscale="Viridis", showscale = False))

        # Update plot sizing
        fig.update_layout(width=800,height=900,autosize=False,margin=dict(t=0, b=0, l=0, r=0),template="plotly_white",)

        # Update 3D scene options
        fig.update_scenes(aspectratio=dict(x=1, y=1, z=0.7),aspectmode="manual")

        # Add dropdown
        fig.update_layout(updatemenus=[dict(buttons=list([dict(args=["type", "surface"],label="3D Surface",method="restyle"),dict(args=["type", "heatmap"],label="Heatmap",method="restyle")]),
                direction="down", pad={"r": 10, "t": 10}, showactive=True, x=0.1, xanchor="left", y=1.1,yanchor="top"),])
        fig.write_html(f'heatmap.html',full_html=False,include_plotlyjs='cdn')
        # fig.write_json(f'heatmap.json')
        fig.show()

    if args.split_images:
        sqdiff = pd.read_csv(args.html)
        rcParams['figure.figsize'] = 30,0.46*len(sqdiff['ref_id'].unique())
        for id in sqdiff["ref_id"].unique():
            df = sqdiff.loc[sqdiff["ref_id"] == id]
            df = df.set_index(['ref_id','split_id'])
            heatmap = sns.heatmap(data=df, annot=True, cmap="afmhot", cbar=False)
            plt.xticks(rotation=90)
            plt.yticks(rotation=0)
            plt.savefig(f"{args.split_images}/heatmap_{id}.svg")
            plt.clf()
