
def CLINVAR():
    # SETUP
    import os
    directory = os.getcwd()
    for a in os.listdir(''.join([directory, '/results/'])):
        os.remove(os.path.join(''.join([directory, '/results/']), a))
    # refining features extracted with terminal code
    # re-organizing the list
    print('opening the file...')
    ClinVar_features = open(''.join([directory, '/input/ClinVar_features.txt']), 'r').read().splitlines()
    ClinVar_features_refined = []
    count = []
    for a in list(range(len(ClinVar_features))):
        if '<ClinVarSet ID=' in ClinVar_features[a]:
            c = a + 1
            while c < len(ClinVar_features) and '<ClinVarSet ID=' not in ClinVar_features[c]:
                c = c + 1
            new = ClinVar_features[a:a + 4]
            position = []
            for c in ClinVar_features[a + 4:c]:
                if 'SequenceLocation' in c:
                    position.append(c)
            new.append(position)
            ClinVar_features_refined.append(new)
        if a % 1000000 == 0:
            print('processed:', a, 'out of', len(ClinVar_features))
    for a in ClinVar_features_refined:
        if '  <Title>' not in a[1] or 'AND ' not in a[1]:
            print('error: Title')
        if '    <ClinVarAccession Acc=' not in a[2]:
            print('error: ClinVarAccession Acc')
        if '      <Description>' not in a[3]:
            print('error: Description')
    del count
    del new
    del position
    del ClinVar_features
    # data overview
    # finding all types of pathological conditions
    pathological_types = []
    for a in ClinVar_features_refined:
        if a[3] not in pathological_types:
            pathological_types.append(a[3])
    print('types of pathological labels:', len(pathological_types))
    del pathological_types
    # finding out how many "Pathogenic" or "Benign" variants contain "not provided" or cancer-related words
    # the variants are then eliminated
    print('filtering...')
    count = []
    no_words = ['not provided', 'cancer', 'Cancer', 'CANCER', 'carcinoma', 'Carcinoma', 'CARCINOMA', 'tumor', 'Tumor',
                'TUMOR', 'adenocarcinoma', 'Adenocarcinoma', 'ADENOCARCINOMA', 'myeloma', 'Myeloma', 'MYELOMA',
                'retinoblastoma', 'Retinoblastoma', 'RETINOBLASTOMA', 'leukemia', 'Leukemia', 'LEUKEMIA', 'lymphoma',
                'Lymphoma', 'LYMPHOMA', 'sarcoma', 'Sarcoma', 'SARCOMA', 'melanoma', 'Melanoma', 'MELANOMA', 'glioma',
                'Glioma', 'GLIOMA', 'glioblastoma', 'Glioblastoma', 'GLIOBLASTOMA']
    for a in list(range(len(ClinVar_features_refined))):
        if ('Pathogenic' in ClinVar_features_refined[a][3] or 'Benign' in ClinVar_features_refined[a][3]):
            found = 'no'
            for c in no_words:
                if c in ClinVar_features_refined[a][1].split('AND ')[1].split('<')[0]:
                    found = 'yes'
            if found == 'yes':
                count.append(ClinVar_features_refined[a][1].split('AND ')[1].split('<')[0])
                ClinVar_features_refined[a] = []
        else:
            ClinVar_features_refined[a] = []
        if a % 100000 == 0:
            print('filtering', a, 'out of', len(ClinVar_features_refined))
    print('variants with no pathology or cancer-related pathology:', len(count), '(eliminated)')
    del found
    del count
    del no_words
    # checking the quality of position annotation
    # all variants that have no indication of "referenceAllele" and "alternateAllele" are discarded
    count = []
    for a in list(range(len(ClinVar_features_refined))):
        if ClinVar_features_refined[a] != []:
            found = 'no'
            for c in list(range(len(ClinVar_features_refined[a][4]))):
                if 'referenceAllele' in ClinVar_features_refined[a][4][c] and \
                        'alternateAllele' in ClinVar_features_refined[a][4][c]:
                    found = 'yes'
                else:
                    ClinVar_features_refined[a][4][c] = []
            if found == 'no':
                count.append(ClinVar_features_refined[a][4])
                ClinVar_features_refined[a] = []
    print('variants with no position annotation:', len(count), '(eliminated)')
    del found
    del count
    # checking variants with multiple mutations and deleting them
    count = []
    for a in list(range(len(ClinVar_features_refined))):
        if ClinVar_features_refined[a] != []:
            found = []
            for c in list(range(len(ClinVar_features_refined[a][4]))):
                if 'referenceAllele' in ClinVar_features_refined[a][4][c] and \
                        'alternateAllele' in ClinVar_features_refined[a][4][c] and \
                        ClinVar_features_refined[a][4][c] not in found:
                    found.append(ClinVar_features_refined[a][4][c])
            if len(found) > 1:
                count.append(found)
                ClinVar_features_refined[a] = []
    print('variants with multiple mutations:', len(count), '(eliminated)')
    del found
    del count
    # checking the length of the mutations
    # all mutations that might affect relative position from the 5'-ss (e.g. insertions and deletions) are discarded
    single_mutations = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']
    count = []
    for a in list(range(len(ClinVar_features_refined))):
        if ClinVar_features_refined[a] != []:
            found = []
            for c in list(range(len(ClinVar_features_refined[a][4]))):
                if ClinVar_features_refined[a][4][c] != [] and \
                        (ClinVar_features_refined[a][4][c].split('referenceAllele')[1].split('=')[1].split('"')[
                             1] not in single_mutations or \
                         ClinVar_features_refined[a][4][c].split('alternateAllele')[1].split('=')[1].split('"')[
                             1] not in single_mutations):
                    found.append(
                        [ClinVar_features_refined[a][4][c].split('referenceAllele')[1].split('=')[1].split('"')[1],
                         ClinVar_features_refined[a][4][c].split('alternateAllele')[1].split('=')[1].split('"')[1]])
            if len(found) > 0:
                count.append(found)
                ClinVar_features_refined[a] = []
    print('variants with non-single nucleotide mutations:', len(count), '(eliminated)')
    del single_mutations
    del count
    del found
    # checking variants with ambiguous position annotation
    # position annotation is then condensed
    count = []
    for a in list(range(len(ClinVar_features_refined))):
        if ClinVar_features_refined[a] != []:
            found = []
            for c in list(range(len(ClinVar_features_refined[a][4]))):
                if ClinVar_features_refined[a][4][c] != []:
                    found.append(ClinVar_features_refined[a][4][c])
            if len(found) != 1 or \
                    found[0].split('start="')[1].split('"')[0] != found[0].split('stop="')[1].split('"')[0]:
                count.append(ClinVar_features_refined[a])
            else:
                ClinVar_features_refined[a][4] = found[0]
    print('variants with unclear position annotation:', len(count))
    del found
    del count
    # checking gene annotation through dictionaries created from the annotation file
    # the keys of the dictionaries are the gene name or the position
    # if the gene name or the variant position is not cleared it gets fixed
    print('preparing annotation')
    annotation = open(''.join([directory, '/input/annotation.gtf']),
                      'r').read().splitlines()
    del annotation[0:5]
    annotation_name = {}
    annotation_position = {}
    count = 0
    for a in annotation:
        if a.split('\t')[2] == 'gene':
            if a.split('\t')[8].split('"')[5] not in annotation_name:
                annotation_name[a.split('\t')[8].split('"')[5]] = ''
            if a.split('\t')[0] not in annotation_position:
                annotation_position[a.split('\t')[0]] = []
            annotation_position[a.split('\t')[0]].append([int(a.split('\t')[3]),
                                                          int(a.split('\t')[4]),
                                                          a.split('\t')[8].split('"')[5]])
        count = count + 1
        if count % 1000000 == 0:
            print('preparing annotation:', count, 'out of', len(annotation))
    del annotation
    count = []
    for a in list(range(len(ClinVar_features_refined))):
        if ClinVar_features_refined[a] != []:
            if '(' not in ClinVar_features_refined[a][1].split('AND ')[0] or \
                    ClinVar_features_refined[a][1].split('(')[1].split(')')[0] not in annotation_name:
                count.append(ClinVar_features_refined[a][1])
                position = [ClinVar_features_refined[a][4].split('Chr=')[1].split('"')[1],
                            int(ClinVar_features_refined[a][4].split('start="')[1].split('"')[0])]
                for d in annotation_position[position[0]]:
                    if position[1] >= d[0] and position[1] <= d[1]:
                        ClinVar_features_refined[a][1] = 'AND '.join([''.join(['(', d[2], ')']),
                                                                      ClinVar_features_refined[a][1].split('AND ')[1]])
                        break
    print('variants with no clear gene annotation:', len(count), '(fixed)')
    count = []
    for a in list(range(len(ClinVar_features_refined))):
        if ClinVar_features_refined[a] != []:
            if '(' not in ClinVar_features_refined[a][1].split('AND ')[0] or \
                    ClinVar_features_refined[a][1].split('(')[1].split(')')[0] not in annotation_name:
                count.append(ClinVar_features_refined[a][1])
                ClinVar_features_refined[a] = []
    print('variants with no clear gene annotation (after fixing):', len(count), '(eliminated)')
    del annotation_name
    del annotation_position
    del position
    del count
    # extracting pathogenic and benign vatiants with single nucleotide variations affecting annotated genes
    # the following features are extracted:
    # ClinVarSet ID
    # ClinVarAccession Acc
    # disease
    # gene
    # position
    # some variants are duplicate, they are condensed together
    print('condensig the dataset...')
    ClinVar_features_refined_pathogenicVSbenign = []
    duplicates = {}
    count = -1
    for a in ClinVar_features_refined:
        if a != []:
            new = [a[0].split('"')[1],
                   a[2].split('"')[1],
                   a[1].split('AND ')[1].split('<')[0],
                   a[1].split('(')[1].split(')')[0],
                   a[4].split('Chr=')[1].split('"')[1],
                   a[4].split('start=')[1].split('"')[1],
                   '>'.join([a[4].split('start=')[1].split('"')[1],
                             a[4].split('referenceAllele')[1].split('=')[1].split('"')[1],
                             a[4].split('alternateAllele')[1].split('=')[1].split('"')[1]]),
                   a[3].split('>')[1].split('<')[0]]
            # checking for duplicates
            if '>'.join([new[3], new[6]]) not in duplicates:
                ClinVar_features_refined_pathogenicVSbenign.append(new)
                count = count + 1
                duplicates['>'.join([new[3], new[6]])] = count
            else:
                ClinVar_features_refined_pathogenicVSbenign[duplicates['>'.join([new[3], new[6]])]][0] = '_'.join([
                    ClinVar_features_refined_pathogenicVSbenign[duplicates['>'.join([new[3], new[6]])]][0],
                    new[0]])
                ClinVar_features_refined_pathogenicVSbenign[duplicates['>'.join([new[3], new[6]])]][1] = '_'.join([
                    ClinVar_features_refined_pathogenicVSbenign[duplicates['>'.join([new[3], new[6]])]][1],
                    new[1]])
                ClinVar_features_refined_pathogenicVSbenign[duplicates['>'.join([new[3], new[6]])]][2] = '_'.join([
                    ClinVar_features_refined_pathogenicVSbenign[duplicates['>'.join([new[3], new[6]])]][2],
                    new[2]])
                ClinVar_features_refined_pathogenicVSbenign[duplicates['>'.join([new[3], new[6]])]][7] = '_'.join([
                    ClinVar_features_refined_pathogenicVSbenign[duplicates['>'.join([new[3], new[6]])]][7],
                    new[7]])
    del duplicates
    del count
    del new
    del ClinVar_features_refined
    # counting how many variants are Pathogenic, Benign or both
    # total number variants = length ClinVar_features_refined
    count_both = 0
    count_Pathogenic = 0
    count_Benign = 0
    count_duplicate = 0
    for a in ClinVar_features_refined_pathogenicVSbenign:
        if 'Pathogenic' in a[7] and 'Benign' in a[7]:
            count_both = count_both + 1
        elif 'Pathogenic' in a[7]:
            count_Pathogenic = count_Pathogenic + 1
        elif 'Benign' in a[7]:
            count_Benign = count_Benign + 1
        if len(a[0].split('_')) > 1:
            count_duplicate = count_duplicate + 1
    print('total variants:', len(ClinVar_features_refined_pathogenicVSbenign))
    print('Pathogenic:', count_Pathogenic)
    print('Benign:', count_Benign)
    print('both:', count_both)
    print('duplicates:', count_duplicate)
    del count_both
    del count_Pathogenic
    del count_Benign
    del count_duplicate
    # searching for splicing mutations
    # creating a dictionary of genes and exons
    # MT gene and last exons are excluded
    # genes with duplicate names are also excluded
    print('selecting splicing mutations...')
    print('preparing annotation...')
    annotation = open(''.join([directory, '/input/annotation.gtf']),
                      'r').read().splitlines()
    del annotation[0:5]
    single_genes = {}
    duplicate_genes = {}
    for a in annotation:  # creating a filter for duplicate genes
        if a.split('\t')[2] == 'gene':
            if a.split('\t')[8].split('"')[5] not in single_genes:
                single_genes[a.split('\t')[8].split('"')[5]] = ''
            else:
                duplicate_genes[a.split('\t')[8].split('"')[5]] = ''
    del single_genes
    transcripts = {}
    for a in annotation:  # creating a filter for last exons
        if a.split('\t')[2] == 'transcript' and a.split('\t')[8].split('"')[9] not in duplicate_genes:
            if a.split('\t')[8].split('"')[9] not in transcripts:
                transcripts[a.split('\t')[8].split('"')[9]] = {'direction': a.split('\t')[6],
                                                               a.split('\t')[8].split('"')[5]: []}
            else:
                transcripts[a.split('\t')[8].split('"')[9]][a.split('\t')[8].split('"')[5]] = []
    for a in annotation:
        if a.split('\t')[2] == 'exon' and a.split('\t')[8].split('"')[11] in transcripts:
            transcripts[a.split('\t')[8].split('"')[11]][a.split('\t')[8].split('"')[5]].append(
                [int(a.split('\t')[3]), int(a.split('\t')[4])])
    non_last_exons = {}
    for a in transcripts:
        non_last_exons[a] = []
        for b in transcripts[a]:
            if b != 'direction':
                transcripts[a][b] = sorted(transcripts[a][b], key=lambda x: x[0])
                if transcripts[a]['direction'] == '+':
                    for e in list(range(len(transcripts[a][b]) - 1)):
                        if [str(transcripts[a][b][e][0]), str(transcripts[a][b][e][1])] not in non_last_exons[a]:
                            non_last_exons[a].append([str(transcripts[a][b][e][0]), str(transcripts[a][b][e][1])])
                else:
                    for e in list(range(1, len(transcripts[a][b]))):
                        if [str(transcripts[a][b][e][0]), str(transcripts[a][b][e][1])] not in non_last_exons[a]:
                            non_last_exons[a].append([str(transcripts[a][b][e][0]), str(transcripts[a][b][e][1])])
    del transcripts
    dictionary = {}
    count = 0
    for a in annotation:  # collecting genes and exons
        if a.split('\t')[2] == 'gene' and a.split('\t')[8].split('"')[5] not in duplicate_genes:
            if a.split('\t')[8].split('"')[5] not in dictionary:
                dictionary[a.split('\t')[8].split('"')[5]] = [
                    [a.split('\t')[6], a.split('\t')[0], a.split('\t')[3], a.split('\t')[4],
                     a.split('\t')[8].split('"')[1]], []]
            else:
                dictionary[a.split('\t')[8].split('"')[5]][0] = [a.split('\t')[6], a.split('\t')[0], a.split('\t')[3],
                                                                 a.split('\t')[4], a.split('\t')[8].split('"')[1]]
        elif a.split('\t')[2] == 'exon' and a.split('\t')[8].split('"')[11] not in duplicate_genes:
            if a.split('\t')[8].split('"')[11] not in dictionary and [a.split('\t')[3], a.split('\t')[4]] in \
                    non_last_exons[a.split('\t')[8].split('"')[11]]:
                dictionary[a.split('\t')[8].split('"')[11]] = [[], [[a.split('\t')[3], a.split('\t')[4]]]]
            elif a.split('\t')[8].split('"')[11] in dictionary and [a.split('\t')[3], a.split('\t')[4]] in \
                    non_last_exons[a.split('\t')[8].split('"')[11]]:
                dictionary[a.split('\t')[8].split('"')[11]][1].append([a.split('\t')[3], a.split('\t')[4]])
        count = count + 1
        if count % 1000000 == 0:
            print('preparing annotation:', count, 'out of', len(annotation))
    del annotation
    del duplicate_genes
    del non_last_exons
    # scanning for SNP variants falling between -20 and +20 from 5'-splice sites
    # variants reported to be both "Pathogenic" and "Benign" are discarded
    # the positions in "range" are relative to the 3' base of the exon (which is 0)
    print('searching splicing mutations...')
    target = [-19, 21]
    ClinVar_features_refined_pathogenicVSbenign_splicing = []
    count = 0
    for a in ClinVar_features_refined_pathogenicVSbenign:
        if a[3] in dictionary and not ('Pathogenic' in a[7] and 'Benign' in a[7]):
            gene = dictionary[a[3]]
            new = [a[3], a[1], ''.join(['chr', gene[0][1]]), '', '', gene[0][0], '', '', a[6], '', '', a[2], '', [],
                   a[7]]
            # looking for positions around the 5'-ss
            splicesites = []
            if gene[0][0] == '+':
                for d in gene[1]:
                    if int(a[5]) in [z + int(d[1]) for z in list(range(target[0], target[1]))] and \
                            d not in new[13] and d[1] not in splicesites:
                        new[13].append(d)
                        splicesites.append(d[1])
                # keeping only non-ambigous exons
                if new[13] != [] and len(splicesites) == 1:
                    new[3] = new[13][0][0]
                    new[4] = new[13][0][1]
                    new[7] = str(int(new[13][0][1]) - int(new[13][0][0]) + 1)
                    new[8] = '_'.join([new[13][0][1],
                                       str(int(a[5]) - int(new[13][0][1])),
                                       '>'.join([new[8].split('>')[1], new[8].split('>')[2]])])
                    del new[13]
                    ClinVar_features_refined_pathogenicVSbenign_splicing.append(new)
            elif gene[0][0] == '-':
                for d in gene[1]:
                    if int(a[5]) in [int(d[0]) - z for z in list(range(target[0], target[1]))] and \
                            d not in new[13] and d[0] not in splicesites:
                        new[13].append(d)
                        splicesites.append(d[0])
                # keeping only non-ambigous exons
                if new[13] != [] and len(splicesites) == 1:
                    new[3] = new[13][0][0]
                    new[4] = new[13][0][1]
                    new[7] = str(int(new[13][0][1]) - int(new[13][0][0]) + 1)
                    new[8] = '_'.join([new[13][0][0],
                                       str(int(new[13][0][0]) - int(a[5])),
                                       '>'.join([new[8].split('>')[1], new[8].split('>')[2]])])
                    del new[13]
                    ClinVar_features_refined_pathogenicVSbenign_splicing.append(new)
        count = count + 1
        if count % 10000 == 0:
            print('searching splicing mutations...', count, 'out of', len(ClinVar_features_refined_pathogenicVSbenign))
    del count
    del target
    del gene
    del new
    del splicesites
    del ClinVar_features_refined_pathogenicVSbenign
    del dictionary
    # outputting the data about mutation distribution
    distribution_list = []
    distribution_dictionary = {}
    for a in list(range(-19, 21)):
        distribution_list.append(str(a))
        distribution_dictionary[str(a)] = [[], []]
    for a in ClinVar_features_refined_pathogenicVSbenign_splicing:
        if 'Pathogenic' in a[13]:
            distribution_dictionary[a[8].split('_')[1]][0].append('_'.join([a[0], a[8]]))
        elif 'Benign' in a[13]:
            distribution_dictionary[a[8].split('_')[1]][1].append('_'.join([a[0], a[8]]))
    output = open(''.join([directory, '/results/distribution.txt']), 'x')
    output.write('position')
    output.write('\t')
    output.write('# Pathogenic')
    output.write('\t')
    output.write('# Benign')
    output.write('\n')
    for a in distribution_list:
        if int(a) < 1:
            output.write(str(int(a) - 1))
        else:
            output.write(a)
        output.write('\t')
        output.write(str(len(distribution_dictionary[a][0])))
        output.write('\t')
        output.write(str(len(distribution_dictionary[a][1])))
        output.write('\n')
    del distribution_list
    del distribution_dictionary
    del output
    # choosing candidates for U1 treatment
    # making a dictionary for the assembly
    # in the initial list chromosome X and Y are in position 23 and 24 respectively
    print('choosing candidate exons for U1 treatment: preparing the assembly...')
    genome = open(''.join([directory, '/input/assembly.fna']),
                  'r').read().split('>')
    for a in list(range(len(genome))):
        if genome[a][0:2] != 'NC':
            genome[a] = []
    genome[len(genome) - 1] = []
    while [] in genome:
        genome.remove([])
    genome.insert(0, [])
    for a in list(range(1, len(genome))):
        genome[a] = genome[a].splitlines()
        genome[a].remove(genome[a][0])
        genome[a] = ''.join(genome[a])
    dictionary = {}
    for a in list(range(1, len(genome))):
        if a < 23:
            dictionary[str(a)] = genome[a]
        elif a == 23:
            dictionary['X'] = genome[a]
        elif a == 24:
            dictionary['Y'] = genome[a]
    del genome
    # selecting pathogenic variants and extracting the sequence from -31 to +110 from the 5'-ss
    # also, mitochondrial genes are eliminated
    # only variants falling in the selected positions are chosen
    # the positions are referred to the 3'-end of the exon (which is 0)
    # for some reason this part of the algorithm can be run only once (after then it gives errors)
    print('choosing candidate exons for U1 treatment: selecting the exons...')
    ClinVar_features_refined_pathogenicVSbenign_splicing_candidates = []
    positions = ['-3', '-2', '-1', '3', '4', '5', '6', '7']
    count = []
    for a in ClinVar_features_refined_pathogenicVSbenign_splicing:
        if 'Benign' not in a[13] and \
                a[2].split('chr')[1] != 'MT' and \
                a[8].split('_')[1] in positions:
            # extracting the sequence
            if a[5] == '+':
                sequence = dictionary[a[2].split('chr')[1]][int(a[4]) - 31:int(a[4]) + 110]
                sequence = list(sequence)
                for d in list(range(len(sequence))):
                    if sequence[d] == 'a':
                        sequence[d] = 'A'
                    elif sequence[d] == 't':
                        sequence[d] = 'T'
                    elif sequence[d] == 'c':
                        sequence[d] = 'C'
                    elif sequence[d] == 'g':
                        sequence[d] = 'G'
                sequence = ''.join(sequence)
            elif a[5] == '-':
                sequence = dictionary[a[2].split('chr')[1]][int(a[3]) - 111:int(a[3]) + 30]
                sequence = list(sequence)
                sequence.reverse()
                for d in list(range(len(sequence))):
                    if sequence[d] == 'A' or sequence[d] == 'a':
                        sequence[d] = 'T'
                    elif sequence[d] == 'T' or sequence[d] == 't':
                        sequence[d] = 'A'
                    elif sequence[d] == 'C' or sequence[d] == 'c':
                        sequence[d] = 'G'
                    elif sequence[d] == 'G' or sequence[d] == 'g':
                        sequence[d] = 'C'
                sequence = ''.join(sequence)
            # mutating the sequence
            # with genes on the - strands, the mutation annotation referes to the + strand anyway
            if a[5] == '+':
                if sequence[30 + int(a[8].split('_')[1])] != a[8].split('_')[2].split('>')[0]:
                    count.append(a)
                else:
                    sequence = list(sequence)
                    sequence[30 + int(a[8].split('_')[1])] = a[8].split('_')[2].split('>')[1]
                    a[9] = ''.join(sequence)
                    ClinVar_features_refined_pathogenicVSbenign_splicing_candidates.append(a)
            elif a[5] == '-':
                if a[8].split('_')[2].split('>')[0] == 'A':
                    a[8] = '_'.join([a[8].split('_')[0],
                                     a[8].split('_')[1],
                                     '>'.join(['T', a[8].split('_')[2].split('>')[1]])])
                elif a[8].split('_')[2].split('>')[0] == 'T':
                    a[8] = '_'.join([a[8].split('_')[0],
                                     a[8].split('_')[1],
                                     '>'.join(['A', a[8].split('_')[2].split('>')[1]])])
                elif a[8].split('_')[2].split('>')[0] == 'C':
                    a[8] = '_'.join([a[8].split('_')[0],
                                     a[8].split('_')[1],
                                     '>'.join(['G', a[8].split('_')[2].split('>')[1]])])
                elif a[8].split('_')[2].split('>')[0] == 'G':
                    a[8] = '_'.join([a[8].split('_')[0],
                                     a[8].split('_')[1],
                                     '>'.join(['C', a[8].split('_')[2].split('>')[1]])])
                if a[8].split('_')[2].split('>')[1] == 'A':
                    a[8] = '_'.join([a[8].split('_')[0],
                                     a[8].split('_')[1],
                                     '>'.join([a[8].split('_')[2].split('>')[0], 'T'])])
                elif a[8].split('_')[2].split('>')[1] == 'T':
                    a[8] = '_'.join([a[8].split('_')[0],
                                     a[8].split('_')[1],
                                     '>'.join([a[8].split('_')[2].split('>')[0], 'A'])])
                elif a[8].split('_')[2].split('>')[1] == 'C':
                    a[8] = '_'.join([a[8].split('_')[0],
                                     a[8].split('_')[1],
                                     '>'.join([a[8].split('_')[2].split('>')[0], 'G'])])
                elif a[8].split('_')[2].split('>')[1] == 'G':
                    a[8] = '_'.join([a[8].split('_')[0],
                                     a[8].split('_')[1],
                                     '>'.join([a[8].split('_')[2].split('>')[0], 'C'])])
                if sequence[30 + int(a[8].split('_')[1])] != a[8].split('_')[2].split('>')[0]:
                    count.append(a)
                else:
                    sequence = list(sequence)
                    sequence[30 + int(a[8].split('_')[1])] = a[8].split('_')[2].split('>')[1]
                    a[9] = ''.join(sequence)
                    ClinVar_features_refined_pathogenicVSbenign_splicing_candidates.append(a)
    for a in list(range(len(ClinVar_features_refined_pathogenicVSbenign_splicing_candidates))):
        ClinVar_features_refined_pathogenicVSbenign_splicing_candidates[a].remove(
            ClinVar_features_refined_pathogenicVSbenign_splicing_candidates[a][13])
    print('selected variants:', len(ClinVar_features_refined_pathogenicVSbenign_splicing_candidates))
    print('variants with non-corresponding sequence:', len(count), '(eliminated)')
    del count
    del sequence
    del positions
    del dictionary
    del ClinVar_features_refined_pathogenicVSbenign_splicing
    # outputting the exons
    ClinVar_features_refined_pathogenicVSbenign_splicing_candidates.insert(0,
                                                                           ['gene', 'ClinVar id', 'chr', 'start', 'end',
                                                                            'strand', 'TechNOA', 'length', 'mutations',
                                                                            'mutated sequence -31>+110', 'effect',
                                                                            'note', 'reference'])
    output = open(''.join([directory, '/results/exons.txt']), 'x')
    for a in ClinVar_features_refined_pathogenicVSbenign_splicing_candidates:
        for b in a:
            output.write(b)
            output.write('\t')
        output.write('\n')
    del output
CLINVAR()