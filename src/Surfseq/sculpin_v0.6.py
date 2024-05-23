# Benjamin Adam Catching
# 2023-03-02
# NIH-NIAID-LVD
# Quantitative Virology and Evolution Unit

"""
Mottled sculpin (Myoxocephalus scorpius), an Atlantic ocean fish, v0.6
From aligned, phred-score = 20 cutoff, Q20.txt files, call each consensus.
mutation in each replicate of a 8-replicate, 12-condition passaging experiment
Returns a single .csv file with information for each nucleotide-position.
"""

# US/MO/14-18947 consensus sequence
MO_aa = 'MGAQVTRQQTGTHENANIATNGSHITYNQINFYKDSYAASASKQDFSQDPSKFTEPVVEGLKAGAPVLKSPSAEACGYSDRVLQLKLGNSAIVTQEAANYCCAYGEWPNYLPDHEAVAIDKPTQPETATDRFYTLKSVKWETGSTGWWWKLPDALNNIGMFGQNVQHHYLYRSGFLIHVQCNATKFHQGALLVVAIPEHQRGAHNTNTSPGFDDIMKGEEGGTFNHPYVLDDGTSLACATIFPHQWINLRTNNSATIVLPWMNAAPMDFPLRHNQWTLAIIPVVPLGTRTTSSMVPITVSIAPMCCEFNGLRHAITQGVPTYLLPGSGQFLTTDDHSSAPALPCFNPTPEMHIPGQVRNMLEVVQVESMMEINNTESAVGMERLKVDISALTDVDQLLFNIPLDIQLDGPLRNTLVGNISRYYTHWSGSLEMTFMFCGSFMATGKLILCYTPPGGSCPTTRETAMLGTHIVWDFGLQSSVTLIIPWISGSHYRMFNNDAKSTNANVGYVTCFMQTNLIVPSESSDTCSLIGFIAAKNDFSLRLMRDSPDIGQLDHLHAAEAAYQIESIIKTATDTVKSEINAELGVVPSLNAVETGVTSNTEPEEAIQTRTVINQHGVSETLVENFLSRAALVSKRSFEYKDHTSSTARADKNFFKWTINTRSFVQLRRKLELFTYLRFDAEITILTTVAVNGSGNNTYVGLPDLTLQAMFVPTGALTPEKQDSFHWQSGSNASVFFKISDPPARITIPFMCINSAYSVFYDGFAGFEKNGLYGINPADTIGNLCVRIVNEHQPVGFTVTVRVYMKPKHIKAWAPRPPRTLPYMSIANANYKGKQRAPNALSAIIGNRDSVKTMPHNIVNTGPGFGGVFVGSFKIINYHLATTEERQSAIYVDWQSDVLVTPIAAHGRHQIARCKCNTGVYYCRHKNRSYPICFEGPGIQWIEQNEYYPARYQTNVLLAVGPAEAGDCGGLLVCPHGVIGLLTAGGGGIVAFTDIRNLLWLDTDAMEQGITDYIQNLGNAFGAGFTETISNKAKEVQDMLIGESSLLEKLLKALIKIISALVIVIRNSEDLVTVTATLALLGCHDSPWSYLKQKVCSYLGIPYVPRQGESWLKKFTEACNALRGLDWLSQKIDKFINWLKTKILPEAREKYEFVQRLKQLPVIENQVSTIEHSCPTTEQQQALFNNVQYYSHYCRKYAPLYAVEAKRVVALEKKINNYIQFKSKSRIEPVCLIIHGSPGTGKSVASNLIARAITEKLGGDIYSLPPDPKYFDGYKQQTVVLMDDLMQNPDGNDISMFCQMVSTVDFIPPMASLEEKGTLYTSPFLIATTNAGSIHAPTVSDSKALSRRFKFDVDIEVTDSYKDSNKLDMSRAVEMCKPDGCAPTNYKRCCPLICGKAIQFRDRRTNARSTIDMLVTDIIKEYRTRNSTQDKLEALFQGPPQFKEIKISVTPDTPAPDAINDLLRSVDSQEVRDYCQKKGWIVVHPSNELIVEKHISRAFITLQAIATFVSIAGVVYVIYKLFAGIQGPYTGIPNPKPKVPSLRTAKVQGPGFDFAQAIMKKNTVIARTEKGEFTMLGVYDRVAVIPTHASVGETIYINDVETKVLDACALRDLTDTNLEITIVKLDRNQKFRDIRHFLPRYEDDYNDAVLSVHTSKFPNMYIPVGQVTNYGFLNLGGTPTHRILMYNFPTRAGQCGGVVTTTGKVIGIHVGGNGAQGFAAMLLHSYFSDTQGEIVSSEKSGVCINAPAKTKLQPSVFHQVFEGSKEPAVLNPKDPRLKTDFEEAIFSKYTGNKIMLMDEYMEEAVDHYVGCLEPLDISVDPIPLESAMYGMDGLEALDLTTSAGFPYLLQGKKKRDIFNRHTRDTSEMTKMLEKYGVDLPFVTFVKDELRSREKVEKGKSRLIEASSLNDSVAMRVAFGNLYATFHNNPGTATGSAVGCDPDIFWSKIPILLDGEIFAFDYTGYDASLSPVWFACLKKVLIKLGYTHQTSFIDYLCHSVHLYKDKKYIVNGGMPSGSSGTSIFNTMINNIIIRTLLIRVYKGIDLDQFKMIAYGDDVIASYPHKIDPGLLAEAGKQYGLVMTPADKGTSFIDTNWENVTFLKRYFRADDQYPFLIHPVMPMKEIHESIRWTKDPRNTQDHVRSLCYLAWHNGEEAYNEFCRKIRSVPVGRALTLPAYSSLRRKWLDSF*'
# Fermon translation
fermon_aa = 'MGAQVTRQQTGTHENANIATNGSHITYNQINFYKDSYAASASKQDFSQDPSKFTEPVVEGLKAGAPVLKSPSAEACGYSDRVLQLKLGNSAIVTQEAANYCCAYGEWPNYLPDHEAVAIDKPTQPETSTDRFYTLRSVKWESNSTGWWWKLPDALNNIGMFGQNVQYHYLYRSGFLIHVQCNATKFHQGALLVVAIPEHQRGAHDTTTSPGFNDIMKGERGGTFNHPYVLDDGTSIACATIFPHQWINLRTNNSATIVLPWMNVAPMDFPLRHNQWTLAVIPVVPLGTRTMSSVVPITVSIAPMCCEFNGLRHAITQGVPTYLLPGSGQFLTTDDHSSAPVLPCFNPTPEMHIPGQIRNMLEMIQVESMMEINNTDGANGMERLRVDISVQADLDQLLFNIPLDIQLDGPLRNTLVGNISRYYTHWSGSLEMTFMFCGSFMATGKLILCYTPPGGSCPTTRETAMLGTHIVWDFGLQSSITLIIPWISGSHYRMFNSDAKSTNANVGYVTCFMQTNLIVPSESSDTCSLIGFIAAKDDFSLRLMRDSPDIGQSNHLHGAEAAYQVESIIKTATDTVKSEINAELGVVPSLNAVETGATSNTEPEEAIQTRTVINQHGVSETLVENFLGRAALVSKKSFEYKNHASSSAGTHKNFFKWTINTKSFVQLRRKLELFTYLRFDAEITILTTVAVNGNNDSTYMGLPDLTLQAMFVPTGALTPKEQDSFHWQSGSNASVFFKISDPPARMTIPFMCINSAYSVFYDGFAGFEKNGLYGINPADTIGNLCVRIVNEHQPVGFTVTVRVYMKPKHIKAWAPRPPRTMPYMSIANANYKGRDTAPNTLNAIIGNRASVTTMPHNIVTTGPGFGGVFVGSFKIINYHLATIEERQSAIYVDWQSDVLVTPIAAHGRHQIARCKCNTGVYYCRHRDRSYPICFEGPGIQWIEQNEYYPARYQTNVLLAAGPAEAGDCGGLLVCPHGVIGLLTAGGGGIVAFTDIRNLLWLDTDVMEQGITDYIQNLGNAFGAGFTETISNKAKEVQDMLIGESSLLEKLLKALIKIISALVIVIRNSEDLITVTATLALLGCHDSPWSYLKQKVCSYLGIPYVPRQSESWLKKFTEACNALRGLDWLSQKIDKFINWLKTKILPEAREKYEFVQRLKQLPVIEKQVSTIEHSCPTTERQQALFNNVQYYSHYCRKYAPLYAVESKRVAALEKKINNYIQFKSKSRIEPVCLIIHGSPGTGKSVASNLIARAITEKLGGDIYSLPPDPKYFDGYKQQTVVLMDDLMQNPDGNDISMFCQMVSTVDFIPPMASLEEKGTLYTSPFLIATTNAGSIHAPTVSDSKALSRRFKFDVDIEVTDSYKDSNKLDMSRAVEMCKPDNCTPTNYKRCCPLICGKAIQFRDRRTNARSTVDMLVTDIIKEYRTRNSTQDKLEALFQGPPQFKEIKISVAPDTPAPDAINDLLRSVDSQEVRDYCQKKGWIVIHPSNELVVEKHISRAFITLQAIATFVSIAGVVYVIYKLFAGIQGPYTGIPNPKPKVPSLRTAKVQGPGFDFAQAIMKKNTVIARTEKGEFTMLGVYDRVAVIPTHASVGEIIYINDVETRVLDACALRDLTDTNLEITIVKLDRNQKFRDIRHFLPRCEDDYNDAVLSVHTSKFPNMYIPVGQVTNYGFLNLGGTPTHRILMYNFPTRAGQCGGVVTTTGKVIGIHVGGNGAQGFAAMLLHSYFTDTQGEIVSNEKSGMCINAPAKTKLQPSVFHQVFEGSKEPAVLNSKDPRLKTDFEEAIFSKYTGNKIMLMDEYMEEAVDHYVGCLEPLDISVDPIPLENAMYGMEGLEALDLTTSAGFPYLLQGKKKRDIFNRQTRDTSEMTKMLEKYGVDLPFVTFVKDELRSREKVEKGKSRLIEASSLNDSVAMRVAFGNLYATFHNNPGTATGSAVGCDPDIFWSKIPILLDGEIFAFDYTGYDASLSPVWFACLKKVLIKLGYTHQTSFIDYLCHSVHLYKDRKYVINGGMPSGSSGTSIFNTMINNIIIRTLLIKVYKGIDLDQFKMIAYGDDVIASYPHKIDPGLLAEAGKHYGLVMTPADKGTSFIDTNWENVTFLKRYFRADDQYPFLIHPVMPMKEIHESIRWTKDPRNTQDHVRSLCYLAWHNGEEAYNEFCRKIRSVPVGRALTLPAYSSLRRKWLDSF*'

# Define function for extracting consensus and >50% mutants
def Q20toconsensus(q20_df):
    """
    Read through dataframe of position nucleotide counts
    """
    
    pos2nuc = {0:'A', 1:'C', 2:'G', 3:'T'}
    # Make list of max nucleotides
    consensus_list = []
    sums = []
    for i in range(len(q20_df)):
        if sum(q20_df.iloc[i]) >= 3:
            new_nuc = pos2nuc[list(q20_df.iloc[i]).index(max(list(q20_df.iloc[i])))]
        else:
            new_nuc = 'N'
        sums.append(q20_df.iloc[i].sum())
        consensus_list.append(new_nuc)
    return consensus_list, sums

# Extract locations of all Q20 files
P1_files = sorted(glob.glob('../../data/sequences/2023-02-06_passage_1/Q20/*.txt'))
P2_files = sorted(glob.glob('../../data/sequences/2023-02-27_passage_2/Q20/*.txt'))
P3_files = sorted(glob.glob('../../data/sequences/2023-03-14_passage_3/Q20/*.txt'))
P4_files = sorted(glob.glob('../../data/sequences/2023-04-06_passage_4/Q20/*.txt'))
P5_files = sorted(glob.glob('../../data/sequences/2023-05-30_passage_5/Q20/*.txt'))
P6_files = sorted(glob.glob('../../data/sequences/2023-06-18_passage_6/Q20/*.txt'))
P1_df = pd.read_csv(P1_files[1], delimiter='\t')
P1_df.columns = ['A', 'C', 'G', 'T', 'N']

# Iterate over all Q20 files, open, then find consensus, add to list
P1_consensus_list = []
for i in range(96):
    #i += 1
    if i % 8 == 0:
        print(i)
    # Define next file
    temp_file = P1_files[i]
    # Open file
    temp_df = pd.read_csv(temp_file, delimiter='\t')
    temp_df.columns = ['A', 'C', 'G', 'T', 'N']
    # Get consensus
    temp_consensus = Q20toconsensus(temp_df)
    #print(''.join(temp_consensus[0]))
    P1_consensus_list.append(''.join(temp_consensus[0]))

P2_consensus_list = []
for i in range(96):
    #i += 1
    if i % 8 == 0:
        print(i)
    # Define next file
    temp_file = P2_files[i]
    # Open file
    temp_df = pd.read_csv(temp_file, delimiter='\t')
    temp_df.columns = ['A', 'C', 'G', 'T', 'N']
    # Get consensus
    temp_consensus = Q20toconsensus(temp_df)
    #print(''.join(temp_consensus[0]))
    P2_consensus_list.append(''.join(temp_consensus[0]))

P3_consensus_list = []
for i in range(96):
    #i += 1
    if i % 8 == 0:
        print(i)
    # Define next file
    temp_file = P3_files[i]
    # Open file
    temp_df = pd.read_csv(temp_file, delimiter='\t')
    temp_df.columns = ['A', 'C', 'G', 'T', 'N']
    # Get consensus
    temp_consensus = Q20toconsensus(temp_df)
    #print(''.join(temp_consensus[0]))
    P3_consensus_list.append(''.join(temp_consensus[0]))

P4_consensus_list = []
for i in range(96):
    #i += 1
    if i % 8 == 0:
        print(i)
    # Define next file
    temp_file = P4_files[i]
    # Open file
    temp_df = pd.read_csv(temp_file, delimiter='\t')
    temp_df.columns = ['A', 'C', 'G', 'T', 'N']
    # Get consensus
    temp_consensus = Q20toconsensus(temp_df)
    #print(''.join(temp_consensus[0]))
    P4_consensus_list.append(''.join(temp_consensus[0]))

P5_consensus_list = []
for i in range(96):
    #i += 1
    if i % 8 == 0:
        print(i)
    # Define next file
    temp_file = P5_files[i]
    # Open file
    temp_df = pd.read_csv(temp_file, delimiter='\t')
    temp_df.columns = ['A', 'C', 'G', 'T', 'N']
    # Get consensus
    temp_consensus = Q20toconsensus(temp_df)
    #print(''.join(temp_consensus[0]))
    P5_consensus_list.append(''.join(temp_consensus[0]))

P6_consensus_list = []
for i in range(96):
    #i += 1
    if i % 8 == 0:
        print(i)
    # Define next file
    temp_file = P6_files[i]
    # Open file
    temp_df = pd.read_csv(temp_file, delimiter='\t')
    temp_df.columns = ['A', 'C', 'G', 'T', 'N']
    # Get consensus
    temp_consensus = Q20toconsensus(temp_df)
    #print(''.join(temp_consensus[0]))
    P6_consensus_list.append(''.join(temp_consensus[0]))


mutant_P1 = []
gene_starts = np.array([69, 317, 564, 861, 1008, 1107, 1437, 1526, 1548, 1731, 2188])
genes = ['VP4', 'VP2', 'VP3', 'VP1', '2A', '2B', '2C', '3A', '3B', '3C', '3D']

for passage in range(0, 6, 1):
    print(passage)
    current_list_of_consensus = [P1_consensus_list, P2_consensus_list, P3_consensus_list, P4_consensus_list, 
                                 P5_consensus_list, P6_consensus_list][passage]
    count = 0
    for i, cell in enumerate(['A549', 'RD', 'SH-SY5Y']):
        for j, temperature in enumerate([33, 37]):
            for k, strain in enumerate(['US/MO/14-18947', 'Fermon']):
                temporary_reference = [MO_aa, fermon_aa][k]
                for l, rep in enumerate('ABCDEFGH'):
                    #print(count)
                    # Translate depending on the reference
                    if strain == 'US/MO/14-18947':
                        temporary_sequence = str(Seq.Seq(current_list_of_consensus[count][697:7264]).translate())
                    else:
                        temporary_sequence = str(Seq.Seq(current_list_of_consensus[count][731:7301]).translate())
                    #print(len(temporary_sequence), len(temporary_reference))
                    for i in range(2187):
                        if i >= 68:
                            ith_gene = sum(i > gene_starts)
                            current_gene = genes[ith_gene]
                            gene_position = i - gene_starts[ith_gene-1] + 1
                        else:
                            current_gene = 'VP4'
                            gene_position = int(i+1)
                        
                        if temporary_reference[i] != temporary_sequence[i] and temporary_sequence[i] != 'X':
                            mutation = temporary_reference[i] + str(int(gene_position)) + temporary_sequence[i]
                            mutant_P1.append([cell, temperature, strain, rep, temporary_reference[i], int(i+1), temporary_sequence[i], current_gene, gene_position, passage+1, mutation])
                    count += 1
P1_mutant_df = pd.DataFrame(mutant_P1)
P1_mutant_df.columns = ['cell', 'temperature', 'strain', 'replicate', 'WT amino acid', 'polyprotein position', 'mutant amino acid', 'gene', 'gene position', 'passage', 'mutation']

P1_mutant_df.to_csv('../../data/consensus_mutations.csv')