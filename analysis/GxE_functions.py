# Benjamin Adam Catching
# 2022-10-24
# NIH-NIAID-LVD-QVEU
# EV-D68 GxE Project

# Misc functions

# Define TCID50 titer function
def calc_titer(temp_dataframe, base_dilution=0):
    titers = []
    for line_num in range(len(temp_dataframe)):
        #print(temp_dataframe.iloc[line_num])
        dils = [12] + [int(x) / 12 for x in temp_dataframe.iloc[line_num][3:]]
        effective_dil = 0
        #print(dils)
        for i in range(len(dils)):
            if dils[i] >= .5:
                if dils[i+1] <= .5:
                    prev_dil = dils[i]
                    next_dil = dils[i+1]
                    effect_dil = (prev_dil - 0.5) / (prev_dil - next_dil)
                    #print(i, prev_dil, next_dil, effect_dil)
                    titers.append(10 ** (effect_dil + i + base_dilution))
                    break
            else:
                #print(i)
                titers.append(0)#np.NaN)
                break
    return titers

def mut_count(list_of_list_of_muts):
    temp_dict = {}
    for rep in list_of_list_of_muts:
        for mut in rep:
            if mut in temp_dict.keys():
                temp_dict[mut] += 1
            else:
                temp_dict[mut] = 1
    return temp_dict

def co_frequencies(df_1, df_2, mut_pair_1, mut_pair_2):
    co_frequency = []
    for passage in range(6):
        passage += 1
        #print(passage)
        #pass_df = passage_df[(passage_df['strain']=='Fermon') & (passage_df['passage'] == passage)]
        for cell in ['RD', 'A549', 'SH-SY5Y']:
            #cell_df = pass_df[pass_df['cell'] == cell]
            for temperature in [33, 37]:
                for replicate in 'ABCDEFGH':
                    mut1_df = df_1[(df_1['temperature'] == temperature) & 
                                   (df_1['cell'] == cell) &
                                   (df_1['passage'] == passage) &
                                   (df_1['protein'] == '3D') &
                                   (df_1['protein position'] == int(mut_pair_1[0])) &
                                   (df_1['mutant aa'] == mut_pair_1[1]) &
                                   (df_1['replicate'] == replicate)]['nucleotide percent'].to_numpy()
                    
                    mut2_df = df_2[(df_2['temperature'] == temperature) & 
                                   (df_2['cell'] == cell) &
                                   (df_2['passage'] == passage) &
                                   (df_2['protein'] == '3D') &
                                   (df_2['protein position'] == int(mut_pair_2[0])) &
                                   (df_2['mutant aa'] == mut_pair_2[1]) &
                                   (df_2['replicate'] == replicate)]['nucleotide percent'].to_numpy()
                    
                    if len(mut1_df) != 0 and len(mut2_df) != 0:
                        freq_1 = float(mut1_df[0])
                        freq_2 = float(mut2_df[0])
                    elif len(mut1_df) == 0 and len(mut2_df) != 0:
                        freq_1 = 0
                        freq_2 = float(mut2_df[0])
                    elif len(mut1_df) != 0 and len(mut2_df) == 0:
                        freq_1 = float(mut1_df[0])
                        freq_2 = 0
                    else:
                        freq_1 = 0
                        freq_2 = 0
                    co_frequency.append([freq_1, freq_2])
    return np.array(co_frequency)


def computeEnrichment(t,d,P="b",ids="ext_gene",tailed=2):
    """
    param t: term
    param d: input DataFrame
    param P: parameter to test (beta value)
    param ids: mapping between the sleuth and database

    return list of values: p-value, mean of the input parameters, mean variance of the input population,
                            number of overlapping genes in the term list and the input file
    
    """
    
    # Set the threshold for number of minimum number of times the mean reactome beta is larger
    # than the null beta
    threshold=21
    
    # Import gene list for the given reactome
    DBgenes=MSigDB.loc[MSigDB.term==t].genes

    # Mean of all the beta values from Sleuth Wald test
    totalMean=d[P].values.mean()
    # 
    targets=d.loc[d[ids].isin(DBgenes.values)]
    #print(targets)
    nrowS=targets.shape[0]
    if nrowS>2:#If more that two genes in term set... #nrowS: number of overlapping genes in pathway and data. 
        print(t)
        
        #print(meanSelect)
        nrowT=d.shape[0]#nrowT: number of rows in dataframe (number of genes)
        meanVar=targets[P].abs().mean()
        meanSelect=totalMean-targets[P].values.mean()
        positives=0
        denom=0
        if tailed==1:
            if meanSelect < 0:#choose which tail to test...
                while (positives<threshold)&(denom<MaxReps):
                    denom+=1
                    diffSample=totalMean-d[P].iloc[np.random.randint(0,nrowT,nrowS)].mean()
                    positives+=int(diffSample<meanSelect)
            else:
                while (positives<threshold)&(denom<MaxReps):
                    denom+=1                
                    diffSample=totalMean-d[P].iloc[np.random.randint(0,nrowT,nrowS)].mean()
                    positives+=int(diffSample>meanSelect)
        elif tailed==2:
            print("tailed==2")
            while (positives<threshold)&(denom<MaxReps):
                denom+=1
                diffSample=d[P].iloc[np.random.randint(0,nrowT,nrowS)].abs().mean()
                positives+=int(diffSample<meanVar)#flipped this 1/12/24. PD
        else:
            print("Must use 1 or 2 tail.")
        #print(positives/denom)
        #print((meanSelect-totalMean)/totalMean)
    else:
        print("No hits. Skipping Term.")
        return(1.0,1.0,1.0,1.0)
    return([positives/denom, meanSelect, meanVar, nrowS])

#This is the main block, and needs to be edited based on your data and the parameter (score column you will be sampling in your data.)
def f(term, data):
    """
    DOCUMENTATION GOES HERE

    param term: 
    param data:

    return enrichment_df:
    """
    enrichment_df = pd.DataFrame([term,computeEnrichment(term,data,"b")],index=["term","pBS"]).transpose()
    return(enrichment_df)

def stat_results(input_df, output_path='~/IM_AN_ABANDONED_DOGGO_FILE_ADOPT_ME.csv'):
    """
    DOCUMENTATION GOES HERE

    param input_df:
    param output_path

    return stats_df
    """
    stats_df=pd.DataFrame()
    for term in MSigDB.term.unique():
        if term.startswith(PREFIX):
            result=f(term, input_df)
            stats_df=pd.concat([stats_df,
                                pd.DataFrame([term,result.pBS[0][0],
                                result.pBS[0][1],
                                result.pBS[0][2],
                                result.pBS[0][3]]).transpose()], ignore_index=True)
    stats_df.columns=["term","pBS","meanStat","meanVar","N"]   
    pcorr=mult(list(stats_df.pBS), method='fdr_bh')
    stats_df['sig']=pcorr[0]
    stats_df['adjPbs']=pcorr[1]
    stats_df.to_csv(output_path)
    
    return stats_df